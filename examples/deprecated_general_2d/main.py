import os
import sys
import numpy as np

base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path,
                           "abqhom",
                           "src")
sys.path.append(gmsh_path)
sys.path.append(abqhom_path)

# import abaqus homogenization tools
from abqhom.RVE import write_abq_input, finalize_model
from abqhom.abq import average_stress, apply_periodic_bc
from abqhom.utils import reuss, voigt
from abqhom.examples import hole_RVE_2d, fancy_RVE_2d

# import abaqus modules
from abaqus import Mdb, mdb, session
from abaqusConstants import (CARTESIAN, ON, PERCENTAGE, ODB, FULL)

def homogenize_stress(model_name, group_map, strain, E, nu, plane_strain=True, 
                      job_name="homogenization"):
    # write mesh into abaqus input file
    workdir = os.getcwd()
    mesh_file = os.path.join(workdir, model_name + ".inp")
    el_type = "CPE3" if plane_strain else "CPS3"
    write_abq_input(model_name, group_map, mesh_file, el_type)
    
    # import the RVE from input file
    Mdb()
    mdb.ModelFromInputFile(name=model_name, inputFileName=mesh_file)
    del mdb.models["Model-1"]
    model = mdb.models[model_name]
    RVE = model.parts["RVE"]
    all_elements = RVE.sets["ALL_ELEMENTS"]
    
    # create material
    material = model.Material(name="base_material")
    material.Elastic(table=((E, nu), ))
    
    # create section
    model.HomogeneousSolidSection(name="solid_section", 
                                  material="base_material", thickness=None)
    RVE.SectionAssignment(region=all_elements, sectionName="solid_section")
    
    # assembly
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    assembly.Instance(name="RVE-1", part=RVE, dependent=ON)

    # add static load step
    model.StaticStep(name="Load", previous="Initial", 
                     description="apply constraints")
    
    # output requests
    model.FieldOutputRequest(name="field_output", createStepName="Load", 
                              variables=("S", "E", "U", "RF", "TF", "IVOL", 
                                        "NFORC", "SENER"))
    del model.fieldOutputRequests["F-Output-1"]
    del model.historyOutputRequests["H-Output-1"]
    
    # apply periodic boundary conditions
    apply_periodic_bc(model_name, "RVE-1", "Load", group_map, strain)
    
    # job
    mdb.Job(name=job_name, model=model_name, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            nodalOutputPrecision=FULL, resultsFormat=ODB)
    mdb.jobs[job_name].submit()
    mdb.jobs[job_name].waitForCompletion()
    
    # load result file
    odb_path = os.path.join(workdir, job_name + ".odb")
    if odb_path in session.odbs.keys():
        session.odbs[odb_path].close()
    odb = session.openOdb(name=odb_path)
    
    # run averaging procedure
    sig = average_stress(odb, "RVE-1", "Load")
    
    return sig
    
# ---------------------------------------------------------------
# set working directory
workdir = os.path.join(base_path, 
                       "abqhom/examples",
                       "general_2d/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# material parameters
E = 4.0   # Young's modulus
nu = 0.3  # Possion's ratio
plane_strain = True

# RVE size
dx = dy = 10
radius = dx/4
lc = 0.4

# macroscopic strain
eps_star = 0.2

# create geometry in gmsh
model_name, group_map = hole_RVE_2d(dx=dx, dy=dy, radius=radius, lc=lc)
# model_name, group_map = fancy_RVE_2d(dx=dx, dy=dy, radius=radius, lc=lc)

# homogenization routine
C = np.empty((3,3))
for case in range(3):
    eps = np.array([0.0, 0.0, 0.0])
    if case == 2:
        eps_star *= 2
    eps[case] = eps_star

    # compute volume average of stress tensor
    print("\n------------------------------")
    print("Run case {}.".format(case))
    print("strain = " + str(eps.round(4).tolist()))
    job_name = "effective_material_2d_{}".format(case)
    sig = homogenize_stress(model_name, group_map, eps, E, nu, 
                            plane_strain=plane_strain, job_name=job_name)
    print("stress = " + str(sig.round(4).tolist()))
    
    # compute effective stiffness tensor
    C[:,case] = sig/eps_star
    
print("\n------------------------------")
print("Effective constitutive matrix:")
print(C.round(4))

# compute voigt and reuss bounds for RVE with hole
if plane_strain:
    C_bulk = lambda E, nu: E/((1 + nu)*(1 - 2*nu))*np.array([[1 - nu, nu, 0.0],
                                                              [nu, 1 - nu, 0.0],
                                                              [0.0, 0.0, 
                                                              (1 - 2*nu)/2]])
else:
    C_bulk =  lambda E, nu: E/(1-nu**2)*np.array([[1.0, nu, 0.0],
                                                  [nu, 1.0,  0.0],
                                                  [0.0, 0.0, (1 - nu)/2]])
volume_fraction = (dx*dy - np.pi*radius**2)/(dx*dy)
print("\n Voigt upper bound:")
print(voigt(volume_fraction, C_bulk(E, nu), C_bulk(1e-10, 0.0)).round(4))
print("\n Reuss lower bound:")
print(reuss(volume_fraction, C_bulk(E, nu), C_bulk(1e-10, 0.0)).round(4))

# finalize gmsh model
finalize_model()