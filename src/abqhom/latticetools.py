import numpy as np

from abqhom.RVE import find_boundary_elements

def _point_in_cylinder(p1, p2, radius, q):
    axis = p2 - p1
    
    between_circ_faces = (np.dot(q - p1, axis) >= 0 
                          and np.dot(q - p2, axis) <= 0)
    inside_shell = (np.linalg.norm(np.cross(q - p1, axis)) 
                    <= radius*np.linalg.norm(axis))
    inside_cylinder = between_circ_faces and inside_shell
    
    return inside_cylinder

def approx_rel_density(model_name, group_map, node_tags, node_coords, el_tags, 
                       conn, radius, dx, dy, dz):
    # get boundary elements
    boundary_elements = find_boundary_elements(model_name, group_map, el_tags, 
                                               conn)
    face_elements = set(np.hstack([j for i, j in boundary_elements.items() 
                               if "FACE" in i]))
    edge_elements = set(np.hstack([j for i, j in boundary_elements.items() 
                               if "EDGE" in i]))
    
    volume_struts = 0.0
    for e, el_conn in zip(el_tags, conn):
        n1, n2 = node_tags[int(el_conn[0] - 1)], node_tags[int(el_conn[1] - 1)]
        p1, p2 = node_coords[int(n1 - 1)], node_coords[int(n2 - 1)]
        length = np.linalg.norm(p2 - p1)
        
        # get weight factor
        weight = 1.0
        if set([e]) <= face_elements:
            weight = 0.5
        elif set([e]) <= edge_elements:
            weight = 0.25
        
        # update strut volume
        volume_struts += weight*length*np.pi*radius*radius
        
    rel_density = volume_struts/(dx*dy*dz)
    
    return rel_density

def monte_carlo_density(node_tags, node_coords, el_tags, conn, radius, x, y, z, 
                        dx, dy, dz, tol=1e-3, min_points=500, max_points=2000):
    # start monte carlo density integration
    rel_density_old = 0.0
    n_strut_points = 0
    n_total = 0
    change = 2*tol
    while change > tol or n_total < min_points:
        # generate random sample point
        sample_x = np.random.uniform(low=x - dx/2, high=x + dx/2, size=(1,))
        sample_y = np.random.uniform(low=y - dy/2, high=y + dy/2, size=(1,))
        sample_z = np.random.uniform(low=z - dz/2, high=z + dz/2, size=(1,))
        sample = np.concatenate((sample_x, sample_y, sample_z))
        
        # check if the sample point lies in one of the struts
        for e, el_conn in zip(el_tags, conn):
            n1 = node_tags[int(el_conn[0] - 1)]
            n2 = node_tags[int(el_conn[1] - 1)]
            p1, p2 = node_coords[int(n1 - 1)], node_coords[int(n2 - 1)]
            
            if _point_in_cylinder(p1, p2, radius, sample):
                n_strut_points += 1
                break

        # update total point counter
        n_total += 1
            
        # update relative density
        rel_density = n_strut_points/n_total
        if rel_density_old != 0:
            change = abs(rel_density - rel_density_old)/rel_density_old
        rel_density_old = rel_density
        
        if n_total == max_points:
            break
        
    return rel_density, change, n_total