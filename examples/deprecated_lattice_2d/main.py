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
from abqhom.utils import export_csv_file
from abqhom.examples import simple_RVE_2d

# import abaqus modules
from abaqus import Mdb, mdb, session
from abaqusConstants import (CARTESIAN, ON, PERCENTAGE, ODB, FULL,
                             DURING_ANALYSIS, MIDDLE_SURFACE, N1_COSINES)

def homogenize_stress(model_name, group_map, strain, E, nu, d,
                      job_name="homogenization"):
    # write mesh into abaqus input file
    workdir = os.getcwd()
    mesh_file = os.path.join(workdir, model_name + ".inp")
    write_abq_input(model_name, group_map, mesh_file, "B23", lattice=True)
    
    # import the RVE from input file
    Mdb()
    mdb.ModelFromInputFile(name=model_name, inputFileName=mesh_file)
    del mdb.models["Model-1"]
    model = mdb.models[model_name]
    RVE = model.parts["RVE"]
    all_elements = RVE.sets["ALL_ELEMENTS"]
    boundary_elements = RVE.sets["BOUNDARY_ELEMENTS"]
    volume_elements = RVE.sets["VOLUME_ELEMENTS"]
    
    # create material, boundary elements will be assigned half stiffness
    material1 = model.Material(name="strut_material_boundary")
    material1.Elastic(table=((0.5*E, nu), ))
    material2 = model.Material(name="strut_material_volume")
    material2.Elastic(table=((E, nu), ))
    
    # create beam section and profile
    model.CircularProfile(name="beam_profile", r=0.5*d)
    model.BeamSection(name="beam_section_boundary", 
                      integration=DURING_ANALYSIS, profile="beam_profile", 
                      material="strut_material_boundary")
    model.BeamSection(name="beam_section_volume", 
                      integration=DURING_ANALYSIS, profile="beam_profile", 
                      material="strut_material_volume")
    RVE.SectionAssignment(region=boundary_elements, 
                          sectionName="beam_section_boundary", 
                          offsetType=MIDDLE_SURFACE)
    RVE.SectionAssignment(region=volume_elements, 
                          sectionName="beam_section_volume", 
                          offsetType=MIDDLE_SURFACE)
    
    # orientation of beam elements
    # this is always (0, 0, -1) for 2d models
    RVE.assignBeamSectionOrientation(region=all_elements, method=N1_COSINES, 
                                     n1=(0.0, 0.0, -1.0))
    
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
    apply_periodic_bc(model_name, "RVE-1", "Load", group_map, strain,
                      restrict_rotation=True)
    
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
                       "Institute for Mechanics/abqhom/examples",
                       "lattice_foam_2d/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# create folder to store results
result_path = os.path.join(workdir, "..", "results")
if not os.path.isdir(result_path):
    os.mkdir(result_path)

# material parameters
E = 210.0   # Young's modulus
nu = 0.3  # Possion's ratio
aspect_ratios = np.linspace(0.01, 0.8, 101)

# RVE size
dx = dy = 10
lc = 0.4

# macroscopic strain
eps_star = 0.2

# create geometry in gmsh
model_name, group_map = simple_RVE_2d(dx=dx, dy=dy, lc=lc)

# homogenization routine
for ar in aspect_ratios: 
    d = ar*lc  # strut diameter
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
        sig = homogenize_stress(model_name, group_map, eps, E, nu, d,
                                job_name=job_name)
        print("stress = " + str(sig.round(4).tolist()))
        
        # compute effective stiffness tensor
        C[:,case] = sig/eps_star
    
    print("\n------------------------------")
    print("Effective constitutive matrix:")
    print(C.round(4))
    
    # export material tensor
    csv_path = os.path.join(result_path, "result_{}.csv".format(ar))
    export_csv_file(C, csv_path)

# finalize gmsh model
finalize_model()