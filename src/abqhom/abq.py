import numpy as np
import gmsh
from abqhom.RVE import create_boundary_sets
from abqhom.utils import voigt_to_strain_tensor
from abaqus import mdb
from abaqusConstants import THREE_D, TWO_D_PLANAR, DEFORMABLE_BODY, ON, UNSET

BOUNDARY_SETS_2D = ["POINT_1", "POINT_2", "POINT_3", "POINT_4",
                    "EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4"]

BOUNDARY_SETS_3D = ["POINT_1", "POINT_2", "POINT_3", "POINT_4", "POINT_5",
                    "POINT_6", "POINT_7", "POINT_8", "EDGE_1", "EDGE_2", 
                    "EDGE_3", "EDGE_4", "EDGE_5", "EDGE_6", "EDGE_7", "EDGE_8",
                    "EDGE_9", "EDGE_10", "EDGE_11", "EDGE_12", "FACE_1", 
                    "FACE_2", "FACE_3", "FACE_4", "FACE_5", "FACE_6"]

def average_stress(odb, instance_name, step_name, frame=-1):
    """
    Compute the volume averaged stress tensor.

    Parameters
    ----------
    odb : Odb Object
        Abaqus output database object resulting from the microscale simulation.
    instance_name : str
        Name of the RVE instance.
    step_name : step
        Name of the simulation step.
    frame : int, optional
        Identifier of the frame of interest. The default is -1.

    Returns
    -------
    sig : numpy.ndarray
        Volume averaged stress in Voigt notation.

    """
    # get model dimension
    instance = odb.rootAssembly.instances[instance_name]
    dim = 3 if instance.embeddedSpace == THREE_D else 2
    
    # compute RVE size
    if dim == 2:
        dx = np.linalg.norm(np.array(instance.nodeSets["POINT_2"].nodes[0]
                                     .coordinates)
                            - np.array(instance.nodeSets["POINT_1"].nodes[0]
                                       .coordinates))
        dy = np.linalg.norm(np.array(instance.nodeSets["POINT_4"].nodes[0]
                                     .coordinates)
                            - np.array(instance.nodeSets["POINT_1"].nodes[0]
                                       .coordinates))
    elif dim == 3:
        dx = np.linalg.norm(np.array(instance.nodeSets["POINT_5"].nodes[0]
                                     .coordinates)
                            - np.array(instance.nodeSets["POINT_1"].nodes[0]
                                       .coordinates))
        dy = np.linalg.norm(np.array(instance.nodeSets["POINT_4"].nodes[0]
                                     .coordinates)
                            - np.array(instance.nodeSets["POINT_1"].nodes[0]
                                       .coordinates))
        dz = np.linalg.norm(np.array(instance.nodeSets["POINT_1"].nodes[0]
                                     .coordinates)
                            - np.array(instance.nodeSets["POINT_2"].nodes[0]
                                       .coordinates))
    
    # initialize averaged stress tensor
    sig = np.zeros((dim, dim), dtype=float)
    
    # get internal forces
    f1_output = odb.steps[step_name].frames[frame].fieldOutputs["NFORC1"]
    f2_output = odb.steps[step_name].frames[frame].fieldOutputs["NFORC2"]
    if dim == 3:
        f3_output = odb.steps[step_name].frames[frame].fieldOutputs["NFORC3"]
    else:
        f3_output = f1_output
    
    # loop over all boundary nodes
    boundary_sets = BOUNDARY_SETS_3D if dim == 3 else BOUNDARY_SETS_2D
    for set_name in boundary_sets:
        region = instance.nodeSets[set_name]
        f1_sub = f1_output.getSubset(region=region).values
        f2_sub = f2_output.getSubset(region=region).values
        f3_sub = f3_output.getSubset(region=region).values
        
        # collect force values
        f1_dict = {}
        f2_dict = {}
        f3_dict = {}
        for f1, f2, f3 in zip(f1_sub, f2_sub, f3_sub):
            f1_dict[f1.nodeLabel] = f1_dict.get(f1.nodeLabel, 0.0) - f1.data
            f2_dict[f2.nodeLabel] = f2_dict.get(f2.nodeLabel, 0.0) - f2.data
            f3_dict[f3.nodeLabel] = f3_dict.get(f3.nodeLabel, 0.0) - f3.data
            
        # update stress tensor
        for node in region.nodes:
            x = node.coordinates[0:dim]
            f_inner = np.array([f1_dict[node.label], f2_dict[node.label], 
                                f3_dict[node.label]])[0:dim]
            sig += np.outer(f_inner, x)
            
    # divide by RVE volume
    if dim == 2:
        sig = np.array([sig[0,0], sig[1,1], sig[0,1]])/(dx*dy)
        return sig
    
    sig = np.array([sig[0,0], sig[1,1], sig[2,2], 
                    sig[1,2], sig[0,2], sig[0,1]])/(dx*dy*dz)
    return sig

def _create_ref_points(model_name):
    """
    Create reference points for periodic boundary conditions.

    Parameters
    ----------
    model : str
        Name of the model. Should be the same for Abaqus and Gmsh.

    Returns
    -------
    instance_sets : list
        List containing the part names and set names for the reference points.

    """
    model = mdb.models[model_name]
    dim = gmsh.model.getDimension()
    part_dim = THREE_D if dim == 3 else TWO_D_PLANAR
    assembly = model.rootAssembly
    
    # create reference points
    model.Part(dimensionality=part_dim, name="REF_1", type=DEFORMABLE_BODY)
    ref_1 = model.parts["REF_1"]
    ref_1.ReferencePoint(point=(0.0, 0.0, 0.0))
    ref_1.Set(name="BC", referencePoints=(ref_1.referencePoints[1],))
    assembly.Instance(name="REF_1-1", part=ref_1, dependent=ON)
    
    model.Part(dimensionality=part_dim, name="REF_2", type=DEFORMABLE_BODY)
    ref_2 = model.parts["REF_2"]
    ref_2.ReferencePoint(point=(0.0, 0.0, 0.0))
    ref_2.Set(name="BC", referencePoints=(ref_2.referencePoints[1],))
    assembly.Instance(name="REF_2-1", part=ref_2, dependent=ON)
    
    instance_sets = [("REF_1-1", "BC"), ("REF_2-1", "BC")]
    if dim == 2:
        return instance_sets
    
    model.Part(dimensionality=part_dim, name="REF_3", type=DEFORMABLE_BODY)
    ref_3 = model.parts["REF_3"]
    ref_3.ReferencePoint(point=(0.0, 0.0, 0.0))
    ref_3.Set(name="BC", referencePoints=(ref_3.referencePoints[1],))
    assembly.Instance(name="REF_3-1", part=ref_3, dependent=ON)
    
    instance_sets.append(("REF_3-1", "BC"))
    return instance_sets

def _add_equation_constraints(periodic_pairs, reference_points, model_name,
                              instance_name, restrict_rotation):
    """
    Add equation constraints to the Abaqus model.

    Parameters
    ----------
    periodic_pairs : list
        TODO.
    reference_points : list
        TODO.
    model_name : str
        Name of the model. Should be the same for Abaqus and Gmsh.
    instance_name : str
        Name if the instance to apply the constraints on.
    restrict_rotation : bool
        If true, periodic constraints will be applied on rotational degrees
        of freedom.

    Returns
    -------
    None.

    """
    model = mdb.models[model_name]
    dim = gmsh.model.getDimension()
    
    eqn_tag = 1
    for pair, ref in zip(periodic_pairs, reference_points):
        
        # apply constraints on translational degrees of freedom
        for dof in range(1, dim + 1):
            terms = []
            for i in pair:
                terms.append((i[1], instance_name + "." + i[0], dof))
            for j in ref:
                terms.append((j[1], j[0][0] + "." + j[0][1], dof))
            model.Equation(name="periodic_bc_{}_{}".format(eqn_tag, dof), 
                           terms=terms)
            
        # apply constraints on rotational degrees of freedom
        if restrict_rotation:
            rotational_dofs = [4, 5, 6] if dim == 3 else [6]
            for dof in rotational_dofs:
                terms = []
                for i in pair:
                    terms.append((i[1], instance_name + "." + i[0], dof))
                model.Equation(name="periodic_bc_{}_{}".format(eqn_tag, dof), 
                               terms=terms)

        eqn_tag += 1
    
def _boundary_conditions_2d(model_name, instance_name, step_name, eps, 
                            restrict_rotation):
    """
    TODO

    Parameters
    ----------
    model_name : TYPE
        DESCRIPTION.
    instance_name : TYPE
        DESCRIPTION.
    step_name : TYPE
        DESCRIPTION.
    eps : TYPE
        DESCRIPTION.
    restrict_rotation : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    model = mdb.models[model_name]
    assembly = model.rootAssembly
    instance = assembly.instances[instance_name]
    
    ur3 = 0.0 if restrict_rotation else UNSET
    
    # add reference points
    ref_points = _create_ref_points(model_name)
    ref_1 = ref_points[0]
    ref_2 = ref_points[1]
    
    # apply Dirichlet BCs to corner points
    for p in ("POINT_1", "POINT_2", "POINT_3", "POINT_4"):
        region = instance.sets[p]
        x = region.nodes[0].coordinates
        u = np.dot(eps, x)
        model.DisplacementBC(name=p.lower() + "_bc", createStepName=step_name, 
                             region=region, u1=u[0], u2=u[1], ur3=ur3)
        
    # compute RVE size
    dx = np.linalg.norm(np.array(instance.sets["POINT_2"]
                                 .nodes[0].coordinates) - 
                        np.array(instance.sets["POINT_1"]
                                 .nodes[0].coordinates))
    dy = np.linalg.norm(np.array(instance.sets["POINT_4"]
                                 .nodes[0].coordinates) - 
                        np.array(instance.sets["POINT_1"]
                                 .nodes[0].coordinates))
    
    # apply Dirichlet BCs to reference points
    region = assembly.instances[ref_1[0]].sets[ref_1[1]]
    x = np.array([dx, 0.0, 0.0])
    u = np.dot(eps, x)
    model.DisplacementBC(name="ref_1_bc", createStepName=step_name, 
                          region=region, u1=u[0], u2=u[1])
    
    region = assembly.instances[ref_2[0]].sets[ref_2[1]]
    x = np.array([0.0, dy, 0.0])
    u = np.dot(eps, x)
    model.DisplacementBC(name="ref_2_bc", createStepName=step_name, 
                         region=region, u1=u[0], u2=u[1])
    
    # define coupling between edges
    pairs = [[["EDGE_2", 1.0], ["EDGE_4", -1.0]], 
             [["EDGE_3", 1.0], ["EDGE_1", -1.0]]]
    references = [[[ref_1, -1.0]],
                  [[ref_2, -1.0]]]
    
    # add equation constraints
    _add_equation_constraints(pairs, references, model_name, instance_name,
                              restrict_rotation)
    
def _boundary_conditions_3d(model_name, instance_name, step_name, eps, 
                            restrict_rotation):
    """
    TODO

    Parameters
    ----------
    model_name : TYPE
        DESCRIPTION.
    instance_name : TYPE
        DESCRIPTION.
    step_name : TYPE
        DESCRIPTION.
    eps : TYPE
        DESCRIPTION.
    restrict_rotation : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    model = mdb.models[model_name]
    assembly = model.rootAssembly
    instance = assembly.instances[instance_name]
    
    ur1 = ur2 = ur3 = 0.0 if restrict_rotation else UNSET
    
    # add reference points
    ref_points = _create_ref_points(model_name)
    ref_1 = ref_points[0]
    ref_2 = ref_points[1]
    ref_3 = ref_points[2]
    
    # apply Dirichlet BCs to corner points
    for p in ("POINT_1", "POINT_2", "POINT_3", "POINT_4", "POINT_5", "POINT_6",
              "POINT_7", "POINT_8"):
        region = instance.sets[p]
        x = region.nodes[0].coordinates
        u = np.dot(eps, x)
        model.DisplacementBC(name=p.lower() + "_bc", createStepName=step_name, 
                             region=region, u1=u[0], u2=u[1], u3=u[2], ur1=ur1,
                             ur2=ur2, ur3=ur3)
        
    # compute RVE size
    dx = np.linalg.norm(np.array(instance.sets["POINT_5"]
                                 .nodes[0].coordinates) - 
                        np.array(instance.sets["POINT_1"]
                                 .nodes[0].coordinates))
    dy = np.linalg.norm(np.array(instance.sets["POINT_4"]
                                 .nodes[0].coordinates) - 
                        np.array(instance.sets["POINT_1"]
                                 .nodes[0].coordinates))
    dz = np.linalg.norm(np.array(instance.sets["POINT_1"]
                                 .nodes[0].coordinates) - 
                        np.array(instance.sets["POINT_2"]
                                 .nodes[0].coordinates))
    
    # apply Dirichlet BCs to reference points
    region = assembly.instances[ref_1[0]].sets[ref_1[1]]
    x = np.array([dx, 0.0, 0.0])
    u = np.dot(eps, x)
    model.DisplacementBC(name="ref_1_bc", createStepName=step_name, 
                          region=region, u1=u[0], u2=u[1], u3=u[2])
    
    region = assembly.instances[ref_2[0]].sets[ref_2[1]]
    x = np.array([0.0, dy, 0.0])
    u = np.dot(eps, x)
    model.DisplacementBC(name="ref_2_bc", createStepName=step_name, 
                         region=region, u1=u[0], u2=u[1], u3=u[2])
    
    region = assembly.instances[ref_3[0]].sets[ref_3[1]]
    x = np.array([0.0, 0.0, dz])
    u = np.dot(eps, x)
    model.DisplacementBC(name="ref_3_bc", createStepName=step_name, 
                         region=region, u1=u[0], u2=u[1], u3=u[2])
    
    # define coupling between edges
    pairs = [[["FACE_2", 1.0], ["FACE_1", -1.0]], 
             [["FACE_4", 1.0], ["FACE_3", -1.0]],
             [["FACE_6", 1.0], ["FACE_5", -1.0]],
             [["EDGE_5", 1.0], ["EDGE_3", -1.0]],
             [["EDGE_7", 1.0], ["EDGE_1", -1.0]],
             [["EDGE_12", 1.0], ["EDGE_9", -1.0]],
             [["EDGE_10", 1.0], ["EDGE_11", -1.0]],
             [["EDGE_4", 1.0], ["EDGE_6", -1.0]],
             [["EDGE_8", 1.0], ["EDGE_2", -1.0]]]
    references = [[[ref_1, -1.0]],
                  [[ref_2, -1.0]],
                  [[ref_3, -1.0]],
                  [[ref_2, 1.0], [ref_1, -1.0]],
                  [[ref_1, -1.0], [ref_2, -1.0]],
                  [[ref_2, -1.0], [ref_3, -1.0]],
                  [[ref_2, 1.0], [ref_3, -1.0]],
                  [[ref_1, 1.0], [ref_3, -1.0]],
                  [[ref_1, -1.0], [ref_3, -1.0]]]
    
    # add equation constraints
    _add_equation_constraints(pairs, references, model_name, instance_name,
                              restrict_rotation)
    
def apply_periodic_bc(model_name, instance_name, step_name, group_map, strain,
                      restrict_rotation=False):
    """
    TODO

    Parameters
    ----------
    model_name : TYPE
        DESCRIPTION.
    instance_name : TYPE
        DESCRIPTION.
    step_name : TYPE
        DESCRIPTION.
    group_map : TYPE
        DESCRIPTION.
    strain : TYPE
        DESCRIPTION.
    restrict_rotation : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    # read model dimension
    dim = gmsh.model.getDimension()
    
    # convert Voigt to tensor
    if strain.shape == (3,):
        strain = np.array([strain[0], strain[1], 0.0, 0.0, 0.0, strain[2]])
    eps = voigt_to_strain_tensor(strain)
     
    # some shortcuts
    model = mdb.models[model_name]
    assembly = model.rootAssembly
    instance = assembly.instances[instance_name]
    
    # add boundary node sets
    boundary_sets = create_boundary_sets(model_name, group_map)
    for set_name, set_nodes in boundary_sets.items():
        set_name = str(set_name)
        set_nodes = [int(i) for i in set_nodes]
        instance.part.SetFromNodeLabels(name=set_name, nodeLabels=set_nodes, 
                                        unsorted=True)
    
    if dim == 2:
        _boundary_conditions_2d(model_name, instance_name, step_name, eps,
                                restrict_rotation)
        
    else:
        _boundary_conditions_3d(model_name, instance_name, step_name, eps,
                                restrict_rotation)
        
    gmsh.model.setCurrent(current)