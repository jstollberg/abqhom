import os
import sys
import numpy as np

# import gmsh and homogenization package
base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
foamhom_path = os.path.join(base_path, 
                            "Institute for Mechanics", 
                            "abqhom", 
                            "src")
sys.path.append(gmsh_path)
sys.path.append(foamhom_path)
from abqhom.RVE import build_hole_RVE, RVE_to_abaqus_input
from abqhom.utils import average_stress_2d, reuss, voigt

# import abaqus modules
from abaqus import Mdb, mdb, session
from abaqusConstants import (CARTESIAN, ON, PERCENTAGE, ODB, FULL)

def homogenize_stress(eps, size, radius, lc, E, nu, plane_strain=True,
                      model_name="RVE2d", job_name="effective_material_2d"):
    if np.shape(eps) != (3,):
        raise ValueError("strain must be provided in Voigt notation")
    
    # build model in gmsh
    (node_tags, node_coords, 
     el_tags, conn, node_sets) = build_hole_RVE(model_name, 2, size, 
                                                radius, lc)
    
    # write mesh into abaqus input file
    workdir = os.getcwd()
    mesh_file = os.path.join(workdir, model_name + ".inp")
    el_type = "CPE3" if plane_strain else "CPS3"
    RVE_to_abaqus_input(size, node_tags, node_coords, el_tags, conn, el_type,
                        mesh_file)
    
    # import the RVE from input file
    Mdb()
    mdb.ModelFromInputFile(name=model_name, inputFileName=mesh_file)
    del mdb.models["Model-1"]
    model = mdb.models[model_name]
    RVE = model.parts["RVE"]
    all_elements = RVE.sets["ALL_ELEMENTS"]
    
    # dummy nodes to apply periodic BCs
    dummy_1 = model.parts["DUMMY_1"]
    dummy_2 = model.parts["DUMMY_2"]
    del model.parts["DUMMY_3"]
    
    # create material
    material = model.Material(name="base_material")
    material.Elastic(table=((E, nu), ))
    
    # create section
    model.HomogeneousSolidSection(name="solid_section", 
                                  material="base_material", thickness=None)
    RVE.SectionAssignment(region=all_elements, sectionName="solid_section")
    
    # add node sets
    for set_name, set_nodes in node_sets.items():
        set_name = str(set_name)
        set_nodes = [int(i) for i in set_nodes]
        RVE.SetFromNodeLabels(name=set_name, nodeLabels=set_nodes, 
                              unsorted=True)
    
    # assembly
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    assembly.Instance(name="RVE-1", part=RVE, dependent=ON)
    assembly.Instance(name="DUMMY_1-1", part=dummy_1, dependent=ON)
    assembly.Instance(name="DUMMY_2-1", part=dummy_2, dependent=ON)
    
    # add static load step
    model.StaticStep(name="Load", previous="Initial", 
                     description="apply constraints")
    
    # get strain components for boundary conditions
    eps_xx = eps[0]
    eps_yy = eps[1]
    eps_xy = 0.5*eps[2]
    
    # apply BCs to dummy nodes
    region = assembly.instances["DUMMY_1-1"].sets["BC"]
    model.DisplacementBC(name="dummy_1_bc", createStepName="Load", 
                         region=region, u1=eps_xx*size, u2=eps_xy*size)
    region = assembly.instances["DUMMY_2-1"].sets["BC"]
    model.DisplacementBC(name="dummy_2_bc", createStepName="Load", 
                         region=region, u1=eps_xy*size, u2=eps_yy*size)
    
    # apply Dirichlet BCs to corner nodes
    region = assembly.instances["RVE-1"].sets["POINT_1"]
    model.DisplacementBC(name="point_1_bc", createStepName="Load", 
                         region=region, 
                         u1=0.0, u2=0.0)
    region = assembly.instances["RVE-1"].sets["POINT_2"]
    model.DisplacementBC(name="point_2_bc", createStepName="Load", 
                         region=region, u1=eps_xx*size, u2=eps_xy*size)
    region = assembly.instances["RVE-1"].sets["POINT_3"]
    model.DisplacementBC(name="point_3_bc", createStepName="Load", 
                         region=region, u1=eps_xy*size, u2=eps_yy*size)
    region = assembly.instances["RVE-1"].sets["POINT_4"]
    model.DisplacementBC(name="point_4_bc", createStepName="Load", 
                         region=region, u1=(eps_xx + eps_xy)*size, 
                         u2=(eps_xy + eps_yy)*size)
    
    # coupling between left and right edges
    model.Equation(name="periodic_edge_1_1", terms=((1.0, "RVE-1.EDGE_3", 1),
                                                    (-1.0, "RVE-1.EDGE_2", 1),
                                                    (-1.0, "DUMMY_1-1.BC", 1)))
    model.Equation(name="periodic_edge_1_2", terms=((1.0, "RVE-1.EDGE_3", 2),
                                                    (-1.0, "RVE-1.EDGE_2", 2),
                                                    (-1.0, "DUMMY_1-1.BC", 2)))
    
    # coupling between top and bottom edges
    model.Equation(name="periodic_edge_2_1", terms=((1.0, "RVE-1.EDGE_4", 1),
                                                    (-1.0, "RVE-1.EDGE_1", 1),
                                                    (-1.0, "DUMMY_2-1.BC", 1)))
    model.Equation(name="periodic_edge_2_2", terms=((1.0, "RVE-1.EDGE_4", 2),
                                                    (-1.0, "RVE-1.EDGE_1", 2),
                                                    (-1.0, "DUMMY_2-1.BC", 2)))
    
    # output requests
    model.FieldOutputRequest(name="field_output", createStepName="Load", 
                             variables=("S", "E", "U", "RF", "TF", "IVOL", 
                                        "NFORC", "SENER"))
    del model.fieldOutputRequests["F-Output-1"]
    del model.historyOutputRequests["H-Output-1"]
    
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
    boundary_sets = ["POINT_1", "POINT_2", "POINT_3", "POINT_4",
                     "EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4"]
    sig = average_stress_2d(odb, size, boundary_sets, "RVE-1", "Load")
    
    return sig
    
# ---------------------------------------------------------------
# set working directory
workdir = os.path.join(base_path, 
                       "Institute for Mechanics/abqhom/examples/RVE_hole/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# material parameters
E = 4.0   # Young's modulus
nu = 0.3  # Possion's ratio
plane_strain = True

# RVE size
size = 10
radius = size/4
lc = 0.4

# macroscopic strain
eps_star = 0.2

C = np.empty((3,3))
for case in range(3):
    eps = np.array([0.0, 0.0, 0.0])
    if case == 3:
        eps_star *= 2
    eps[case] = eps_star

    # compute volume average of stress tensor
    print("\n------------------------------")
    print("Run case {}.".format(case))
    print("strain = " + str(eps.tolist()))
    job_name = "effective_material_2d_{}".format(case)
    sig = homogenize_stress(eps, size, radius, lc, E, nu, 
                            plane_strain=plane_strain, job_name=job_name)
    print("stress = " + str(sig.tolist()))
    
    # compute effective stiffness tensor
    C[:,case] = sig/eps_star
    
print("\n------------------------------")
print("Effective constitutive matrix:")
print(C)

# compute voigt and reuss bounds
if plane_strain:
    C_bulk = lambda E, nu: E/((1 + nu)*(1 - 2*nu))*np.array([[1 - nu, nu, 0.0],
                                                             [nu, 1 - nu, 0.0],
                                                             [0.0, 0.0, 
                                                              (1 - 2*nu)/2]])
else:
    C_bulk =  lambda E, nu: E/(1-nu**2)*np.array([[1.0, nu, 0.0],
                                                  [nu, 1.0,  0.0],
                                                  [0.0, 0.0, (1 - nu)/2]])
volume_fraction = (size**2 - np.pi/radius**2)/size**2
print("\n Voigt upper bound:")
print(voigt(volume_fraction, C_bulk(E, nu), C_bulk(1e-10, 0.0)))
print("\n Reuss lower bound:")
print(reuss(volume_fraction, C_bulk(E, nu), C_bulk(1e-10, 0.0)))