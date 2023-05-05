# config files for the knot untangling
# options:
#   - geometry: 'file' or 'torus'
#   - filename: path to the .obj file
#   - frame folder: path to the folder where the frames are saved
#   - visuals: 'clean', 'arrows', 'gradient computation'
#   - figure setup: '2d 2x1', '2d 1x1', '3d' 
#   - scattered: True or False 
#   - flip: True or False (switch Y and Z coordinates (blender uses Y as up, but we want Z as up))
#   - tolerance: float (tolerance for the "unknottedness"/energy)
#   - max iterations: int (number of iterations)
#   - n: int (number of vertices in case of generated geometry)
#   - p: int (torus knot)
#   - q: int (torus knot)
#   - r1: float (inner radius)
#   - r2: float (outer radius)
#   - varying: True or False (dynamic time step optimization)
#   - clean plotting: True or False 
#   - angle: tuple (only for 3d plots)
#   - headless: True or False 
#   - title: 'n', 'energy', 'dt', TODO: 'dt/energy', 'dt/n', 'dt/energy/n
#   - name: string (only for 3d plots)


default_config = {
    'geometry': 'file', 
    'filename': './obj/reefknot2.obj',
    'frame folder': './frames/reefknot',
    'figure setup': '2d 2x1',
    'visuals': 'clean',
    'scattered': True,
    'flip': True,       # switch Y and Z coordinates (blender uses Y as up, but we want Z as up)
    'tolerance': 0.1,    # tolerance for the "unknottedness"/energy
    'max iterations': 5_000,
}

# trefoil
trefoil = {
    'geometry': 'torus',
    'frame folder': './frames/torus',
    'figure setup': '2d 2x1',
    'visuals': 'clean',
    'n': 100,             # number of vertices
    'p': 2,              # p/q is the torus knot
    'q': 3,
    'r1': 2,             # inner radius
    'r2': 1,             # outer radius
    'tolerance': 0.08,    # tolerance for the "unknottedness"/energy
    'max iterations': 400,
}

# torus 4 1
torus_4_1 = {
    'geometry': 'torus',
    'frame folder': './frames/torus',
    'figure setup': '3d',
    'clean plotting': True,
    'scattered': True,
    'n': 50,             # number of vertices
    'p': 4,              # p/q is the torus knot
    'q': 1,
    'r1': 8,             # inner radius
    'r2': 1,             # outer radius
    'varying': True,
    'tolerance': 0.1,    # tolerance for the "unknottedness"/energy
    'max iterations': 3_500
}


# overhand config
overhand = {
    'geometry': 'file',
    'headless': True,
    'filename': './obj/overhand.obj',
    'frame folder': './frames/overhand',
    'figure setup': '2d 2x1',
    'angle': (-25, 35, 180),
    'visuals': 'clean',
    'scattered': True,
    'flip': True,        # switch Y and Z coordinates (blender uses Y as up, but we want Z as up)
    'tolerance': 0.1,    # tolerance for the "unknottedness"/energy
    'max iterations': 601,
}

# arrows hour glass
hourglass_gradient = {
    'geometry': 'file',
    'filename': './obj/hourglass.obj',
    'headless': True,   # live plotting reduces speed
    'frame folder': './hourglass',
    'visuals': 'arrows',
    'figure setup': '2d 1x1',
    'tolerance': 0.99,    # tolerance for the "unknottedness"/energy
    'max iterations': 2001,
}

# energy visualization hourglass
hourglass_energy = {
    'geometry': 'file',
    'filename': './obj/hourglass.obj',
    'title': 'energy',
    'headless': True,   # live plotting reduces speed
    'frame folder': './energy_frames/hourglass/computation_with_energy',
    'visuals': 'gradient computation',
    'figure setup': '2d 1x1',
}

# reef knot 1
reefknot1 = {
    'geometry': 'file',
    'filename': './obj/reefknot1.obj',
    'headless': True,   # live plotting reduces speed
    'frame folder': './frames/reefknot',
    'visuals': 'clean',
    'scattered': True,
    'flip': True,       # switch Y and Z coordinates (blender uses Y as up, but we want Z as up)
    'figure setup': '2d 2x1',
    'angle': (20, 80, 0),
    'tolerance': 0.1,    # tolerance for the "unknottedness"/energy
    'max iterations': 1501,
}

# reef knot 2 
reefknot2 = {
    'title': 'n',   
    'name': 'Reef Knot 2',
    'geometry': 'file',
    'headless': True,   # live plotting reduces speed
    'filename': './obj/reefknot2.obj',
    'frame folder': './frames/reefknot/reefknot2_final_01',
    'visuals': 'clean',
    'scattered': True,
    'flip': True,       # switch Y and Z coordinates (blender uses Y as up, but we want Z as up)
    'figure setup': '2d 2x1',
    'tolerance': 0.01,    # tolerance for the "unknottedness"/energy
    'max iterations': 2001,
}
