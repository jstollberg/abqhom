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
from abqhom.RVE import make_lattice_RVE
from abqhom.utils import RVE_to_abaqus_input

# import abaqus modules
from abaqus import Mdb, mdb, session
from abaqusConstants import (DURING_ANALYSIS, CARTESIAN, ON, MIDDLE_SURFACE,
                              N1_COSINES, PERCENTAGE, ODB, FULL, CENTROID)

def average_stress_2d(eps, size, lc, E, nu, d, model_name="RVE2d", 
                      job_name="effective_material_2d"):
    """
    Set up the homogenization BVP and compute the averaged stress tensor.
    
    Periodic boundary conditions are applied depending on the macroscopic
    strain tensor.

    Parameters
    ----------
    eps : numpy.ndarray
        The macroscopic strain tensor in Voigt notation.
    size : float
        Edge length of the RVE.
    lc : float
        Characteristic cell size.
    E : float
        Young's modulus of the strut material.
    nu : float
        Poisson's ratio of the strut material.
    d : float
        Diameter of the struts.
    model_name : str, optional
        Name of the Gmsh and Abaqus model. The default is "RVE2d".
    job_name : str, optional
        Name of the Abaqus job. The default is "effective_material_2d".

    Returns
    -------
    sig : numpy.ndarray
        Volume averaged stress tensor in Voigt notation.

    """
    if np.shape(eps) != (3,):
        raise ValueError("strain must be provided in Voigt notation")
    
    # build model in gmsh
    (node_tags, node_coords, 
     el_tags, conn, node_sets) = make_lattice_RVE(model_name, 2, size, lc)
    
    # write mesh into abaqus input file
    workdir = os.getcwd()
    mesh_file = os.path.join(workdir, model_name + ".inp")
    RVE_to_abaqus_input(size, node_tags, node_coords, el_tags, conn, "B23",
                        mesh_file)
    
    # import the RVE from input file
    Mdb()
    mdb.ModelFromInputFile(name=model_name, inputFileName=mesh_file)
    del mdb.models["Model-1"]
    model = mdb.models[model_name]
    RVE = model.parts["RVE"]
    all_elements = RVE.sets["ALL_ELEMENTS"]
    # all_nodes = RVE.sets["ALL_NODES"]
    
    # dummy nodes to apply periodic BCs
    dummy_1 = model.parts["DUMMY_1"]
    dummy_2 = model.parts["DUMMY_2"]
    del model.parts["DUMMY_3"]
    
    # create material
    material = model.Material(name="strut_material")
    material.Elastic(table=((E, nu), ))
    
    # create beam section and profile
    model.CircularProfile(name="beam_profile", r=0.5*d)
    model.BeamSection(name="beam_section", integration=DURING_ANALYSIS, 
                      profile="beam_profile", material="strut_material")
    RVE.SectionAssignment(region=all_elements, sectionName="beam_section", 
                          offsetType=MIDDLE_SURFACE)
    
    # orientation of beam elements
    # since we use circular beams, the n1-orientation is not important and we 
    # just use some default parameters here
    RVE.assignBeamSectionOrientation(region=all_elements, method=N1_COSINES, 
                                     n1=(0.0, 0.0, -1.0))
    
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
                         u1=0.0, u2=0.0, ur3=0.0)
    region = assembly.instances["RVE-1"].sets["POINT_2"]
    model.DisplacementBC(name="point_2_bc", createStepName="Load", 
                         region=region, u1=eps_xx*size, u2=eps_xy*size, 
                         ur3=0.0)
    region = assembly.instances["RVE-1"].sets["POINT_4"]
    model.DisplacementBC(name="point_4_bc", createStepName="Load", 
                         region=region, u1=eps_xy*size, u2=eps_yy*size, 
                         ur3=0.0)
    region = assembly.instances["RVE-1"].sets["POINT_3"]
    model.DisplacementBC(name="point_3_bc", createStepName="Load", 
                         region=region, u1=(eps_xx + eps_xy)*size, 
                         u2=(eps_xy + eps_yy)*size, ur3=0.0)
    
    # coupling between left and right edges
    model.Equation(name="periodic_edge_1_1", terms=((1.0, "RVE-1.EDGE_2", 1),
                                                    (-1.0, "RVE-1.EDGE_4", 1),
                                                    (-1.0, "DUMMY_1-1.BC", 1)))
    model.Equation(name="periodic_edge_1_2", terms=((1.0, "RVE-1.EDGE_2", 2),
                                                    (-1.0, "RVE-1.EDGE_4", 2),
                                                    (-1.0, "DUMMY_1-1.BC", 2)))
    model.Equation(name="periodic_edge_1_6", terms=((1.0, "RVE-1.EDGE_2", 6),
                                                    (-1.0, "RVE-1.EDGE_4", 6)))
    
    
    # coupling between top and bottom edges
    model.Equation(name="periodic_edge_2_1", terms=((1.0, "RVE-1.EDGE_3", 1),
                                                    (-1.0, "RVE-1.EDGE_1", 1),
                                                    (-1.0, "DUMMY_2-1.BC", 1)))
    model.Equation(name="periodic_edge_2_2", terms=((1.0, "RVE-1.EDGE_3", 2),
                                                    (-1.0, "RVE-1.EDGE_1", 2),
                                                    (-1.0, "DUMMY_2-1.BC", 2)))
    model.Equation(name="periodic_edge_2_6", terms=((1.0, "RVE-1.EDGE_3", 6),
                                                    (-1.0, "RVE-1.EDGE_1", 6)))
    
    # output requests
    model.FieldOutputRequest(name="field_output", createStepName="Load", 
                             variables=("S", "E", "SE", "U", "RF", "SF", 
                                        "TF", "NFORC", "ESF1", "SENER",
                                        "IVOL", "CF", "NFORCSO"))
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
    # session.viewports["Viewport: 1"].setValues(displayedObject=odb)
    # session.viewports["Viewport: 1"].makeCurrent()
    
    # compute volume average of stress
    forces = (odb.steps["Load"].frames[-1].fieldOutputs["SF"]
              .getSubset(CENTROID))
    sig_xx = sig_yy = sig_xy = 0.0
    for f in forces.values:
        assert all(f.data[1::] == 0.0)
        axial_force = f.data[0]
        
        # compute strut length and direction
        el = f.elementLabel
        n1, n2 = conn[int(np.where(el_tags == el)[0])]
        x1, x2 = node_coords[int(n1 - 1)], node_coords[int(n2 - 1)]
        l = np.linalg.norm(x1 - x2)
        rx, ry, _ = np.abs(x1 - x2)/l
        
        # update stress components
        sig_xx += axial_force*rx*rx*l
        sig_yy += axial_force*ry*ry*l
        sig_xy += axial_force*rx*ry*l
        
    # divide by RVE volume
    sig_xx /= size**2
    sig_yy /= size**2
    sig_xy /= size**2
    sig = np.array([sig_xx, sig_yy, sig_xy])
    print("sig =", sig)
    
    return sig
    
# ---------------------------------------------------------------
# set working directory
workdir = os.path.join(base_path, 
                       "Institute for Mechanics/abqhom/examples/abq")
os.chdir(workdir)

# material parameters
E = 210   # Young's modulus
nu = 0.3  # Possion's ratio

# RVE size
size = 10

# parameter space
diameters = [0.1]
cell_sizes = [5]#[0.39]
# cell_sizes = np.linspace(0.4, 2.0, 30)
#diameters = np.linspace(0.05, 0.5, 50)

# macroscopic strain
eps_star = 0.1

# run homogenization routine
for lc in cell_sizes:
    for d in diameters:
        C = np.empty((3, 3))
        alpha = d/lc  # aspect ratio
        for case in range(3):
            eps = np.array([0.0, 0.0, 0.0])
            if case == 3:
                eps_star *= 2
            eps[case] = eps_star
        
            # compute volume average of stress tensor
            print("Run case {}.".format(case))
            job_name = "effective_material_2d_{}".format(case)
            sig = average_stress_2d(eps, size, lc, E, nu, d, job_name=job_name)
            
            # compute effective stiffness tensor
            C[:,case] = sig/eps_star
            
        path = os.path.join(os.path.dirname(workdir), "C",
                            #"C_{}_{}.csv".format(lc, d))
                            "C_{}.csv".format(alpha))
        with open(path, "w") as file:
            file.write("{},{},{}\n".format(C[0,0], C[0,1], C[0,2]))
            file.write("{},{},{}\n".format(C[1,0], C[1,1], C[1,2]))
            file.write("{},{},{}".format(C[2,0], C[2,1], C[2,2]))
    
        print("Effective constitutive matrix:")
        print(C)