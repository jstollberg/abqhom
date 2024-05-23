import os
import sys
import numpy as np
import random

base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path,
                           "abqhom",
                           "src")
paraqus_path = os.path.join(base_path, "paraqus", "src")
sys.path.append(gmsh_path)
sys.path.append(abqhom_path)
sys.path.append(paraqus_path)

# import abaqus homogenization tools
from abqhom.RVE import write_abq_input, finalize_model
from abqhom.abq import average_stress, apply_periodic_bc
from abqhom.utils import export_csv_file
from abqhom.examples import simple_RVE_3d

# import paraqus to export vtk files
paraqus_installed = False
try:
    from paraqus.abaqus import ODBReader
    from paraqus.writers import BinaryWriter
    paraqus_installes = True
except:
    pass

# import abaqus modules
from abaqus import Mdb, mdb, session
from abaqusConstants import (CARTESIAN, ON, PERCENTAGE, ODB, FULL,
                             DURING_ANALYSIS, MIDDLE_SURFACE, N1_COSINES)
from mesh import MeshElementArray
from regionToolset import Region

def export_vtk(odb_path, model_name, vtk_path):
    reader = ODBReader(odb_path=odb_path,
                       model_name=model_name,
                       instance_names=["RVE-1"],
                       )
    reader.add_field_export_request("U", field_position="nodes")
    reader.add_set_export_request("BOUNDARY_ELEMENTS", set_type="elements",
                                  instance_name="RVE-1")
    vtu_writer = BinaryWriter(vtk_path, clear_output_dir=False)
    instance_models = list(reader.read_instances(step_name="Load",
                                                 frame_index=-1))
    instance_model = instance_models[0]
    vtu_writer.write(instance_model)

def homogenize_stress(model_name, group_map, strain, E, nu, d, vtk_path,
                      job_name="homogenization"):
    # write mesh into abaqus input file
    workdir = os.getcwd()
    mesh_file = os.path.join(workdir, model_name + ".inp")
    write_abq_input(model_name, group_map, mesh_file, "B33", lattice=True)
    
    # import the RVE from input file
    Mdb()
    mdb.ModelFromInputFile(name=model_name, inputFileName=mesh_file)
    del mdb.models["Model-1"]
    model = mdb.models[model_name]
    RVE = model.parts["RVE"]
    volume_elements = RVE.sets["VOLUME_ELEMENTS"]
    
    # create material, boundary elements will be assigned half stiffness
    material1 = model.Material(name="strut_material_edges")
    material1.Elastic(table=((0.25*E, nu), ))
    material2 = model.Material(name="strut_material_faces")
    material2.Elastic(table=((0.5*E, nu), ))
    material3 = model.Material(name="strut_material_volume")
    material3.Elastic(table=((E, nu), ))
    
    # create beam section and profile
    model.CircularProfile(name="beam_profile", r=0.5*d)
    model.BeamSection(name="beam_section_edges", 
                      integration=DURING_ANALYSIS, profile="beam_profile", 
                      material="strut_material_edges")
    model.BeamSection(name="beam_section_faces", 
                      integration=DURING_ANALYSIS, profile="beam_profile", 
                      material="strut_material_faces")
    model.BeamSection(name="beam_section_volume", 
                      integration=DURING_ANALYSIS, profile="beam_profile", 
                      material="strut_material_volume")
    
    for i in range(1, 13):
        set_name = "ELEMENTS_EDGE_{}".format(i)
        region = RVE.sets[set_name]
        RVE.SectionAssignment(region=region, 
                              sectionName="beam_section_edges", 
                              offsetType=MIDDLE_SURFACE)
        
    for i in range(1, 7):
        set_name = "ELEMENTS_FACE_{}".format(i)
        region = RVE.sets[set_name]
        RVE.SectionAssignment(region=region, 
                              sectionName="beam_section_faces", 
                              offsetType=MIDDLE_SURFACE)
        
    RVE.SectionAssignment(region=volume_elements, 
                          sectionName="beam_section_volume", 
                          offsetType=MIDDLE_SURFACE)
    
    # orientation of beam elements
    # since we use circular beam elements we only need to pay attention that
    # the n1 direction differs from the tangent direction
    for e in RVE.elements:
        # compute tangent vector
        conn = e.connectivity
        node1, node2 = RVE.nodes[conn[0]], RVE.nodes[conn[-1]]
        coords1, coords2 = node1.coordinates, node2.coordinates
        tangent = np.array([j - i for (i, j) in zip(coords1, coords2)])
        
        # compute normal vector n1
        a = 0
        while tangent[a] == 0:
            a += 1
        b = a + 1 if a + 1 < len(tangent) else a - 1
        n1 = np.zeros_like(tangent)
        n1[b], n1[a] = tangent[a], -tangent[b]
        n1 /= np.linalg.norm(tangent)
        
        # assign orientation
        arr = MeshElementArray(elements=(e,))
        region = Region(elements=arr)
        RVE.assignBeamSectionOrientation(region=region, method=N1_COSINES, 
                                         n1=n1)
            
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
    
    # export as vtk
    if paraqus_installed:
        export_vtk(odb_path, job_name + "_{}".format(10/lc), vtk_path)
    
    return sig
    
# ---------------------------------------------------------------
# set working directory
workdir = os.path.join(base_path, 
                       "abqhom/examples/lattice_3d/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# create folder to store results
result_path = os.path.join(workdir, "..", "results")
if not os.path.isdir(result_path):
    os.mkdir(result_path)
    
# create folder to store vtk files
vtk_path = os.path.join(workdir, "..", "vtk_output")
if not os.path.isdir(vtk_path):
    os.mkdir(vtk_path)

# RVE size
dx = dy = dz = 10
lc = 1.25

# macroscopic strain
eps_star = 0.2

# create geometry in gmsh
model_name, group_map = simple_RVE_3d(dx=dx, dy=dy, dz=dz, lc=lc)

# homogenization routine
n_samples = 1  # 1000
for i in range(n_samples):
    E = random.uniform(0.01, 300.0)
    nu = random.uniform(0.0, 0.5)
    ar = random.uniform(0.01, 1.0)
    d = ar*lc  # strut diameter
    
    C = np.empty((6,6))
    for case in range(6):
        eps = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        if case == 3:
            eps_star *= 2
        eps[case] = eps_star
    
        # compute volume average of stress tensor
        print("\n------------------------------")
        print("Run case {}.".format(case + 1))
        print("strain = " + str(eps.round(4).tolist()))
        job_name = "effective_material_3d_{}".format(case)
        sig = homogenize_stress(model_name, group_map, eps, E, nu, d, vtk_path,
                                job_name=job_name)
        print("stress = " + str(sig.round(4).tolist()))
        
        # compute effective stiffness tensor
        C[:,case] = sig/eps_star
    
    print("\n------------------------------")
    print("Effective constitutive matrix:")
    print(C.round(4))
    
    # export material tensor
    csv_path = os.path.join(result_path, "result_{}_{}_{}.csv".format(ar, E, nu))
    export_csv_file(C, csv_path)

# finalize gmsh model
finalize_model()