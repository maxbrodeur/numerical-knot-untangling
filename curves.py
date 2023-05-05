from matplotlib import pyplot as plt
import numpy as np
from typing import Union

# UTILITY FUNCTIONS FOR MANIPULATING OBJ KNOTS -------------------------------------------------

# Read obj file and return array of vertices, dictionary of edges and dictionary of backward edges
def read_obj(filename):
    '''# read obj file and return array of vertices, dictionary of edges and dictionary of backward edges
    - vertices: array of vertices [(x,y,z)]
    - edges: dictionary of edges {u:v} where u is the index of the vertex and v is the index of the vertex it connects to
    - backward_edges: dictionary of backward edges {v:u} where v is the index of the vertex and u is the index of the vertex it connects to
    
    Note: the edge dictionnaries contain indices starting from 0, not 1 as in the obj file, the reason is 
    to mach the indices of the vertices array
    '''
    with open(filename) as f:
        lines = f.readlines()
        vertices = [line.split()[1:] for line in lines if line.startswith('v ')]
        vertices = [[float(x) for x in vertex] for vertex in vertices]
        vertices = np.array(vertices)
        edges = [line.split()[1:] for line in lines if line.startswith('l ')]
        # subtract 1 to start index from 0
        new_edges = {int(edge[0])-1:int(edge[1])-1 for edge in edges} 
        backward_edges = {v:k for k,v in new_edges.items()}

        # check for duplicates in edge values (incoming and outgoing)
        if len(new_edges.keys()) < len(edges) or len(backward_edges.keys()) < len(edges): # if keys were hit more than twice, then there are duplicates
            # attempt to fix by inverting duplicate edges, will raise an exception if it failed to fix
            new_edges = fix_obj(edges)

    # sort vertices in the order of the edges
    # check for multiple components (>=2 distinct loops)
    components = [] # list of components of the knot
    component = []
    new_vertices = np.zeros((len(vertices), 3))
    cur = list(new_edges.keys())[0]
    visited = set()
    for i in range(len(vertices)):
        new_vertices[i] = vertices[cur]
        component.append(i)
        visited.add(cur)
        if cur in new_edges:
            next = new_edges[cur]
        # dead end TODO: fix this
            
        # if next is already visited, then we have a new component
        if next in visited:
            components.append(component.copy())
            component = []
            # remove visited edges
            new_edges = {k:v for k,v in new_edges.items() if k not in visited}
            if len(new_edges) > 0:
                cur = list(new_edges.keys())[0]
        else:
            cur = next
    
    # update edges dictionary
    shifted_components = [component[1:] + component[:1] for component in components] # close component loops
    new_edges = {}
    for component, shifted_component in zip(components, shifted_components):
        for u,v in zip(component, shifted_component):
            new_edges[u] = v
    
    backward_edges = {v:k for k,v in new_edges.items()}

    return new_vertices, new_edges, backward_edges, components

def create_edges(gamma):
    '''# create edges from gamma
    - gamma: array of vertices [(x,y,z)]
    
    returns a dictionary of edges {u:v} where u is the index of the vertex and v is the index of the vertex it connects to
    '''
    edges = {}

    # chose most isolated edge as start point (less chance of picking wrong neighbor)
    start = -1
    max_sum_dist = 0
    for u in range(len(gamma)):
        sum_dist = 0
        vertex = gamma[u]
        for v in gamma:
            sum_dist += np.linalg.norm(vertex - v)
        
        if sum_dist > max_sum_dist:
            max_sum_dist = sum_dist
            start = u
    
    # create edges
    # find first neighbor
    u = start
    neighbor = None
    while neighbor != start:
        min_dist = np.inf
        neighbor = None
        for v in range(len(gamma)):
            dist = np.linalg.norm(gamma[u] - gamma[v])
            # if smallest distance and v not already an outgoing neighbor
            if u!=v and dist < min_dist and v not in edges.values():
                min_dist = dist
                neighbor = v
        edges[u] = neighbor
        u = edges[u]

    # sort the vertices and update the edges dictionary
    new_vertices = np.zeros((len(vertices), 3))
    cur = 0
    next = edges[cur]
    for i in range(len(vertices)):
        new_vertices[i] = vertices[cur]
        cur = next
        next = edges[cur]

    # update edges dictionary
    edges = {i:(i+1)%len(vertices) for i in range(len(vertices))}

    return vertices, edges, {v:k for k,v in edges.items()}

# Attempting to fix obj by inverting duplicate edges
def fix_obj(edges):
    new_edges = {}
    dup = []
    for u,v in edges:
        u = int(u) - 1 # subtract 1 to start index from 0
        v = int(v) - 1 
        if u not in new_edges:
            new_edges[u] = v
        elif v not in new_edges:
            new_edges[v] = u
            dup.append(u+1) # add 1 to start index from 1 to match obj file
        # if both u and v are already in new_edges, then there is a duplicate outgoing edge aswell
        else :
            dup.append(u+1)
            dup.append(v+1)
            raise Exception(f'Found duplicate vertices: {set(dup)} \nIn a loop, for each vertex v, there should exist a unique outgoing uv-edge, and for each vertex u there should exists a unique incoming uv-edge.'\
                            f'\nAttempted to fix by inverting duplicate edges, but failed.')

# Create a torus knot
def torus_knot(p=2, q=3, n=200, r1 = 4, r2 = 2):
    ''' # Create a torus knot
    - r1 is the radius of the circle from the origin (z axis)
    - r2 is the radius of the circle around the r1 circle axis
    - p is number of times the knot wraps around the r1 circle
    - q is number of times the knot wraps around the r2 circle
    
    returns vertices and edges

    note: if p and q are not relatively prime, then the torus knot has more than one component
    '''
    phi = np.linspace(0, 2 * np.pi, n+1)
    # remove last element to avoid duplicate at 0 and 2pi
    phi = phi[:-1]
    r = r2 * np.cos(q*phi) + r1
    X = r * np.cos(p*phi)
    Y = r * np.sin(p*phi)
    Z = -r2* np.sin(q*phi)

    edges = {i:(i+1)%n for i in range(n)}

    return np.array([X, Y, Z]).T, edges, {v:k for k,v in edges.items()}

# Closes the loop by connecting the discontinuity
def check_loop(edges):
    ''' # close the loop by connecting discontinuity
    logic: 
    - go forwards until hitting a discontinuity, 
    - go backwards until hitting a discontinuity
    - connect the two discontinuities 

    returns a new edges dictionary with the loop closed'''


    backward_edges = {v:k for k,v in edges.items()}    
        
    first = list(edges.keys())[0]
    v = edges[first]
    # go forwards until hitting a discontinuity
    while v in edges.keys() and v != first:
        v = edges[v]
    # if not full loop, go backwards until hitting a discontinuity
    if v != first:
        w = first
        while w in backward_edges.keys():
            w = backward_edges[w]
        edges[v] = w
        print('Closed loop')
    else:
        print('Already closed loop')
    return edges

def write_obj(filename, vertices):
    with open(filename, 'w') as f:
        X, Y, Z = vertices
        for x,y,z in zip(X, Y, Z):
            f.write('v {} {} {}\n'.format(x,y,z))
        for i in range(1, len(X)):
            f.write('l {} {}\n'.format(i, (i+1)))
        f.write('l {} {}\n'.format(len(X), 1))

def invert_YZ(vertices):
    return np.array([[x, z, y] for x,y,z in vertices])

def set_Z(vertices:np.array, Z:float):
    return np.array([[x, y, Z] for x,y,_ in vertices])

def decimate(vertices:np.array, edges:dict, percentage:float=0.5) -> Union[np.array, dict]:
    '''# Decimate a knot by a percentage
            - percentage = 0.5 means that every other vertex is removed
            - percentage converted to step (0.5 = 2)
            - creates new vertices and edges by stepping through the vertices'''
    # randomly sample new vertices
    new_size = int(len(vertices)*percentage)
    # new_index -> old_index
    old_indices = np.random.choice(len(vertices), size=new_size, replace=False)
    old_indices.sort()
    # old_index -> new_index
    new_indices = {old_index:new_index for new_index, old_index in enumerate(old_indices)}
    
    new_vertices = vertices[old_indices]

    # fix new edges
    new_edges = {}
    for new_index in range(new_size): 
        old_index = old_indices[new_index]
        old_neighbor = edges[old_index]
        # find closest sampled neighbor
        new_neighbor = old_neighbor
        while new_neighbor not in old_indices:
            new_neighbor = edges[new_neighbor]
        # convert old index to new index
        new_neighbor = new_indices[new_neighbor]
        new_edges[new_index] = new_neighbor
        
    return new_vertices, new_edges

if __name__ == '__main__':

    vertices, _, _ = torus_knot(3,2,500)
    X, Y, Z = vertices[:,0], vertices[:,1], vertices[:,2]
    ax = plt.axes(projection='3d')

    plt.plot(X, Y, Z, c='black')
    plt.scatter(X, Y, Z, c=range(len(X)), cmap='twilight_shifted')
    plt.title('Trefoil Knot')


