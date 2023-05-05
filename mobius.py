import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp
import time

# import from other files
from curves import * # geometry utilities
from configs import * # configuration files

# GLOBAL VARIABLES
gamma = edges = backward_edges = config = max_edge_length = fig = axes = components = components_inv = None

# Description of variables:
#   gamma:  the curve itself, a list of N points in R^3, 
#           each point is a 3D numpy array
#
#   edges:  a dictionary representing the edges of the curve, 
#           where edges[u]=v represents the edge uv 
#           where u and v are indices of points in gamma 
#
#   backward_edges: 
#           a dictionary representing the edges of the curve, 
#           where backward_edges[u]=v represents the edge vu 
#           where u and v are indices of points in gamma
#
#   config: a dictionary of configuration parameters - see load_config() method
#
#   max_edge_length: 
#           the maximum length of an edge in the curve
#
#   fig:    the matplotlib figure
#
#   axes:   the matplotlib axes, array
# 
#   components: 
#           a list of lists, 
#           where each list is a component of the curve 
#           Note: two components are disjoint – they do not share any points and aren't linked by edges
# 
#   components_inv: 
#           a dictionary of indices to integers TODO: should just be a list,
#           where components_inv[i] = j if the ith component is linked to the jth component by an edge, 
#           to make it easier to determine if two points are in the same component

# find shortest path from u to v
def curve_distance(u ,v):
    global gamma, edges, backward_edges

    # find the path from u to v
    rpath = []
    lpath = []
    rlength = llength = 0
    rcur = lcur = u
    rnext = lnext = u
    # go around the curve to the right until we reach v
    while rnext != v and lnext != v:
        # step right
        rnext = edges[rcur]
        rlength += np.linalg.norm(gamma[rnext] - gamma[rcur])
        rpath.append(rcur)
        rcur = rnext
        # step left
        lnext = backward_edges[lcur]
        llength += np.linalg.norm(gamma[lnext] - gamma[lcur])
        lpath.append(lcur)
        lcur = lnext

    length = 0
    path = []

    if rnext == v:
        path = rpath
        length = rlength
    elif lnext == v:
        path = lpath
        length = llength

    path.append(v)
    
    return length, path

# compute the mobius energy and gradient contribution of a pair of points
def Mobius_gradient(pair):
    global gamma, edges, backward_edges, config
    
    u,v = pair

    # initialize
    energy = 0
    gradient = np.zeros((len(gamma), 3))
    
    
    # distance term (term A) --------------------------------------------------
    
    # energy term
    A = 1 / np.linalg.norm(gamma[u] - gamma[v]) ** 2 # inverse squared distance
    # gradient term
    distance = np.linalg.norm(gamma[u] - gamma[v])
    gradient[v] +=   2 * (gamma[u] - gamma[v]) / distance ** 4
    gradient[u] +=  -2 * (gamma[u] - gamma[v]) / distance ** 4

    # check components
    if components_inv[u] != components_inv[v]:
        energy += A
        return energy, gradient, distance
    
    # curve distance term (term B) --------------------------------------------
    length, path = curve_distance(u, v)
    B = 1 / length ** 2
    energy += A - B
        
    # gradient terms
    for i in path[1:-1]:
        if i != u and i != v:
            # lI - lI+1 case
            dist_I = np.linalg.norm(gamma[i] - gamma[backward_edges[i]])
            dist_I1 = np.linalg.norm(gamma[edges[i]] - gamma[i])
            der = 2 * (gamma[i] - gamma[backward_edges[i]]) / (length ** 3 * dist_I) \
                - 2 * (gamma[edges[i]] - gamma[i]) / (length ** 3 * dist_I1)
            gradient[i] += der
            
    # u is the first vertex and v is the last vertex in the shortest path
    if path[0] == u:
        # lI+1 case since u is the first vertex and term is |u+1 - u|
        dist = np.linalg.norm(gamma[edges[u]]-gamma[u])
        der = -2 * (gamma[edges[u]]-gamma[u])/(length**3 * dist)
        gradient[u] += der

        # and lI case for v since v is the last vertex and term is |v - v-1|
        dist = np.linalg.norm(gamma[v]-gamma[backward_edges[v]])
        der = 2 * (gamma[v]-gamma[backward_edges[v]])/(length**3 * dist)
        gradient[v] += der
    # v is the first vertex and u is the last vertex
    else:
        # lI case since v is the first vertex and term is |v+1 - v|
        dist = np.linalg.norm(gamma[edges[v]]-gamma[v])
        der = -2 * (gamma[edges[v]]-gamma[v])/(length**3 * dist)
        gradient[v] += der

        # and lI+1 case for u since u is the last vertex and term is |u - u-1|
        dist = np.linalg.norm(gamma[u]-gamma[backward_edges[u]])
        der = 2 * (gamma[u]-gamma[backward_edges[u]])/(length**3 * dist)
        gradient[u] += der

    return energy, gradient, distance      

# gradient descent on the mobius energy
def Mobius_gradient_descent():
    global gamma, config, max_edge_length, edges

    clear_frames()
    compute_max_edge_length()

    tolerance = config['tolerance']
    max_iters = config['max iterations']

    count = 0
    txt = fig.texts[0]
    mp.set_start_method('fork')

    uv_pairs = [(u,v) for u in range(len(gamma)) for v in range(len(gamma)) if u!=v]

    t=0
    true_time = []
    dts = []
    energies = []
    start = time.time()
    dt = config['dt']
    while count < max_iters:
        # launch multiprocessing over all uv pairs
        pool = mp.Pool(processes=mp.cpu_count())
        results = pool.map(Mobius_gradient, uv_pairs)

        # add up the energies and gradients
        energy = 0
        gradient = np.zeros((len(gamma), 3))
        min_dist = np.inf
        for e, g, dist in results:
            energy += e
            gradient += g
            min_dist = min(min_dist, dist)

        # prevent overshooting
        if config['varying']:
            m = gradient.max()
            dt = min_dist*tolerance/m
            # while (dt*gradient).max() > min_dist*tolerance:
            #     dt /= 10
        t+=dt 

        # visualization
        if config['visuals'] == 'clean':
            clean_plot(count)
        elif config['visuals'] == 'arrows':
            plot_arrows(gradient, energy, count, min_dist)
        else:
            fill_var = 30
            norms = np.linalg.norm(gradient, axis=1)
            txt.set_text(f'{"time: "}{t:.8f}'.ljust(fill_var)+
                         f'{"dt: "}{dt:.2e}'.ljust(fill_var)+
                         f'{"energy: "}{energy:.5f}'.ljust(fill_var)+
                         f'{"curve length: "}{curve_length():.3f}'.ljust(fill_var)+
                        f'{"max grad: "}{norms.max():.2f}'.ljust(fill_var))
            draw_angles(norms,gradient)
        
        plt.savefig( f'{config["frame folder"]}/frame' + str(count).zfill(5) + '.pdf' )
        plt.savefig( f'{config["frame folder"]}/frame' + str(count).zfill(5) + '.png' )
        count += 1

        # gradient descent
        gamma -= dt*gradient
        print(f'energy: {energy}', end='\r')
        true_time.append(time.time()-start)
        dts.append(dt)
        energies.append(energy)
    plt.show()
    plot_energy(true_time, energies, dts)



# UTILITY -----------------------------------------------------------------------------------

def curve_length():
    length = 0
    cur = 0
    next = -1
    while next != 0:
        next = edges[cur]
        length += np.linalg.norm(gamma[cur] - gamma[next])
        cur = next
    return length

def compute_max_edge_length():
    global max_edge_length, gamma, edges, backward_edges
    max_edge_length = 0
    for u,v in edges.items():
        max_edge_length = max(max_edge_length, np.linalg.norm(gamma[u]-gamma[v]))
    max_edge_length = max_edge_length * config['max_edge_length_ratio']
    print("max length: ", max_edge_length, end='\r')

# PLOTTING ----------------------------------------------------------------------------------

# visualize the gradient computation
def plot_gradient_computation():
    global gamma, edges, backward_edges, config, fig
    clear_frames()

    plt.figure(fig) # activate the figure

    compute_max_edge_length()

    # energy gradient visualization
    global count
    count = 0

    # initialize the energy
    energy = 0
    # initialize the gradient
    gradient = np.zeros((len(gamma), 3))
    for u in range(len(gamma)):
        v = edges[u]
        while v != u:
            # energy terms, -----------------------------------------------------------
            # only path and length are usefule for actual gradient computation
            length, path = curve_distance(u, v)
            
            pair = (u,v)
            e, g, _ = Mobius_gradient(pair)
            energy += e
            gradient += g
            
            draw_vertices(u, v, gradient,length, path, energy, count)
            count+=1
            v = edges[v]
    
def clean_plot(iters):
    global fig, axes, gamma, config, components 
    X, Y, Z = gamma[:,0], gamma[:,1], gamma[:,2]
    X, Y, Z = np.array([*X,X[0]]), np.array([*Y,Y[0]]), np.array([*Z,Z[0]])
    for ax in axes.flatten():
        # check projection of axis
        ax.cla()
        if config['figure setup'] == '3d 1x1':
            ax.axis('off')
        # TODO: make this more general
        colors = ['black', 'green', 'blue', 'red', 'purple', 'yellow', 'pink', 'brown', 'gray', 'orange']
        if ax.name == '3d':
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])

            # TODO: always plot components
            if len(components) > 1:
                for component in components:
                    # plot each component and close them
                    color = colors.pop(0)
                    ax.plot([*X[component], X[component][0]], [*Y[component], Y[component][0]], [*Z[component], Z[component][0]], color='black', linewidth=1)
                    if 'scattered' in config and config['scattered' ]:
                        ax.scatter(X[component], Y[component], Z[component], color=color, s=20, alpha=1.0)
            else:
                ax.plot(X, Y, Z, color='black', linewidth=1)
                if 'scattered' in config and config['scattered' ]:
                    ax.scatter(X, Y, Z, cmap='twilight_shifted', c=range(len(X)), s=10)
            if 'zlims' in config:
                ax.set_zlim(config['zlims'])
        else:
            ax.set_xticks([])
            ax.set_yticks([])
            if len(components) > 1:
                for component in components:
                    # plot each component and close them
                    color = colors.pop(0)
                    ax.plot([*X[component], X[component][0]], [*Y[component], Y[component][0]], color='black', linewidth=1)
                    if 'scattered' in config and config['scattered']:
                        ax.scatter(X[component], Y[component], color=color, s=20, alpha=1.0)
            else:
                ax.plot(X, Y, color='black', linewidth=1)
                if 'scattered' in config and config['scattered']:
                    ax.scatter(X, Y, cmap='twilight_shifted', c=range(len(X)), s=10)

    title = ''
    if config['title'] == 'dt':
        title = r'\(dt = {}\)'.format(config['dt'])
    elif config['title'] == 'energy':
        title = r'\(\\mathcal{{E}} = {}\)'.format(config['energy'])
    elif config['title'] == 'n':
        title = r'\(n = {}\)'.format(iters)
    
    if config['figure setup'] == '3d 1x1' or config['figure setup'] == '2d 1x1':
        ax.set_title(title)
    else:
        fig.suptitle(title)

    if not config['headless']:
        plt.draw()
        plt.pause(0.001)

def plot_energy(time,energy,dts):
    dts = np.array(dts)
    fig = plt.figure(figsize=(10,10))
    ax = plt.gca()
    ax.plot(time, energy, color='black', alpha=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # log scale xticks
    scat = plt.scatter(time, energy, 
                c=dts, cmap='rainbow', norm=colors.LogNorm(min(dts), max(dts)), # log scale for dt values
                marker='o', s=15, alpha=1)
    xticks = np.linspace(min(time), max(time), 10)
    yticks = np.linspace(min(energy), max(energy), 5)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Energy')
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    # log scale colorbar
    # ts = np.logspace(np.log10(min(dts)), np.log10(max(dts)), 10, endpoint=True)
    if config['varying']:
        cax = fig.add_axes([0.12, 0.92, 0.80, 0.02])
        plt.colorbar(label=r'\(dt\)', orientation='horizontal', pad=0.2, cax=cax)
    else:
        s = "{0:.2g}".format(dts[0])
        if "e" in s:
            base, exponent = s.split("e")
            s = r"\(dt={0} \times 10^{{{1}}}\)".format(base, int(exponent))
        ax.set_title(s)
    name = config['name'].split(' ')
    name = '_'.join(name)
    folder = config['frame folder']
    plt.savefig(f'{folder}/{name}_{config["geometry"]}_energy_vs_time_{time[-1]}.pdf')
    plt.savefig(f'{folder}/{name}_{config["geometry"]}_energy_vs_time_{time[-1]}.png')
    # plt.show()
    
def draw_angles(norms=None, gradients=np.zeros(0)): 
    global gamma, max_edge_length

    X, Y, Z = gamma[:,0], gamma[:,1], gamma[:,2]
    X, Y, Z = np.array([*X, X[0]]), np.array([*Y, Y[0]]), np.array([*Z, Z[0]]) # close the loop
    gx, gy, gz = gradients[:,0], gradients[:,1], gradients[:,2]
    gx, gy, gz = np.array([*gx, gx[0]]), np.array([*gy, gy[0]]), np.array([*gz, gz[0]]) # close the loop
    g = [(gx, gy, gz), (gx, gy, gz), (gx, gy, gz), (gx, gy, gz)]
    data = [(X,Y,Z), (X,Y,Z), (X,Y,Z), (X,Y,Z)]
    titles = ['Normal view', 'Top view', 'Side view', 'Front view']

    # compute_max_edge_length()

    # scatter 3D vertices for all 4 subplots from different angles (normal, top, side, front)
    for ax in axes.flatten():
        ax.cla()
        dat = data.pop(0)
        grad = g.pop(0)
        if len(dat) == 3:
            if 'zlims' in config:
                ax.set_zlim(config['zlims'])
            ax.scatter3D(*dat, c = range(len(X)))
            ax.quiver(*dat, *grad, length=max_edge_length/2, color='black', arrow_length_ratio=0.1, normalize=True)
            ax.plot(*dat, c = 'black')
            # ax.set_xticks([])
            # ax.set_yticks([])
            # ax.set_zticks([])
        elif len(dat) == 2:
            X, Y = dat
            grad = np.array(grad).T
            norms = np.linalg.norm(grad, axis=1).reshape(-1)
            ax.scatter(*dat, c = norms)
            
            grad = grad / np.stack([norms,norms],axis=1)
            grad *= max_edge_length/2

            grad = list(grad.T)
            ax.quiver(*dat, *grad, color='black', width=0.005, headwidth=2, headlength=2, headaxislength=1.5)
            
            cmap = cm.get_cmap('rainbow')
            c = cmap(norms)
            ax.plot(*dat, c = 'black')

        # scatter plot
        ax.set_title(titles.pop(0))

    if not config['headless']:
        plt.draw()
        plt.pause(0.0001)

# delete all frames in ./frames
def clear_frames():
    folder = config['frame folder']
    import os
    for filename in os.listdir(folder):
        if filename.endswith(".pdf") or filename.endswith(".png"):
            os.remove(os.path.join(folder, filename))

# draw the vertices for visualing the energy calculation
# - u,v: the two vertices to draw
# - gradients: the gradient at each vertex
# - length: distance between u and v on the curve
# - path: the shortest path between u and v
def draw_vertices(u,v, gradients, length, path, energy, count):
    global gamma, edges, max_edge_length, fig, axes, components

    X, Y, Z = gamma[:,0], gamma[:,1], gamma[:,2]
    gx, gy, gz = gradients[:,0], gradients[:,1], gradients[:,2]

    norms = np.linalg.norm(gradients, axis=1)

    # plot the path between u and v
    path = np.array([[gamma[u,0],gamma[u,1],gamma[u,2]] for u in path])
    
    # distance between u and v
    distance = np.linalg.norm(gamma[u]-gamma[v])
    uv = np.array([gamma[u],gamma[v]])

    # txt.set_text(f'length: {length:.7f}\n'f'distance: {distance:.7f}\n'f'max edge length: {max_edge_length:.7f}\n'f'energy: {energy:.7f}\n')
    c = np.zeros(len(gamma)) 
    c[u]=1
    c[v]=2
    
    for ax in axes.flatten():
            ax.clear()
            # turn off grid
            ax.grid(False)
            # turn off ticks
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            # turn off axis
            ax.axis('off')

            if len(components)>1:
                for component in components:
                    ax.plot([*X[component], X[component[0]]], [*Y[component], Y[component[0]]], [*Z[component], Z[component[0]]], color='black')
            else:
                ax.plot([*X, X[0]], [*Y, Y[0]], [*Z, Z[0]], color='black')
            
            ax.scatter(X[u], Y[u], Z[u], c='green')
            ax.scatter(X[v], Y[v], Z[v], c='red')
            d = max_edge_length/3
            ax.text(gamma[u,0]+d,gamma[u,1]+d,gamma[u,2],'u')
            ax.text(gamma[v,0]+d,gamma[v,1]+d,gamma[v,2],'v')
            
            ax.plot(path[:,0], path[:,1], path[:,2], color='pink', label=r'\(d(u,v)= {}\)'+f'{length:.3f}')
            ax.plot(uv[:,0], uv[:,1], uv[:,2], color='blue', label=r'\(\|\gamma^{(u)}-\gamma^{(v)}\|=\) '+f'{distance:.3f}')
            if config['title'] == 'energy':
                ax.set_title(r'\(\mathcal{E}=\)'+f'{energy:.3f}')
            plt.legend(loc='upper center', ncol=2)

    if not config['headless']:
        plt.draw()
        plt.pause(0.0001)        

    # save figure as frame with 5 digits using count
    folder = config['frame folder']
    plt.savefig(f'{folder}'f'/frame{count:05d}.pdf')
    plt.savefig(f'{folder}'f'/frame{count:05d}.png')

def plot_arrows(gradients, energy, count, min_dist):
    global gamma, edges, max_edge_length, fig, axes

    X, Y, Z = gamma[:,0], gamma[:,1], gamma[:,2]
    gx, gy, gz = gradients[:,0], gradients[:,1], gradients[:,2]

    plt.rcParams["axes.titley"] = 0.95
    # compute_max_edge_length()

    for ax in axes.flatten():
            ax.clear()
            ax.set_title(r'\(\mathcal{{E}}={energy:.2f}\)'.format(energy=energy))
            # turn off grid
            ax.grid(False)
            # turn off ticks
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            # turn off axis
            ax.axis('off')

            if len(components)>1:
                for component in components:
                    ax.plot([*X[component], X[component[0]]], [*Y[component], Y[component[0]]], [*Z[component], Z[component[0]]], color='black')
            else:
                ax.plot([*X, X[0]], [*Y, Y[0]], [*Z, Z[0]], color='black')

            ax.scatter(X, Y, Z, color='black', s=8)
            ax.quiver(X, Y, Z, gx, gy, gz, length=min_dist , color='blue', normalize=True)

    if not config['headless']:
        plt.draw()
        plt.pause(0.005)

    folder = config['frame folder']
    plt.savefig(f'{folder}'f'/frame{count:05d}.pdf')
    plt.savefig(f'{folder}'f'/frame{count:05d}.png')
    
def load_config(config_):
    global config, gamma, edges, backward_edges, fig, axes, components, components_inv
    config = config_

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica",
        "font.size": 20
    })
    

    # default values
    if 'headless' not in config:
        config['headless'] = False
    if 'dt' not in config:
        config['dt'] = 0.01
    if 'frame folder' not in config:
        config['frame folder'] = './frames'
    if 'tolerance' not in config:
        config['tolerance'] = 1e-2
    if 'varying' not in config:
        config['varying'] = True
    if 'figure setup' not in config:
        config['figure setup'] = '3d'
    if 'clean plotting' not in config:
        config['clean plotting'] = False
        config['ignore edges'] = False
    if 'max iterations' not in config:
        config['max iterations'] = 99_999
    if 'title' not in config:
        config['title'] = ''

    # load geometry
    geom = ''
    if 'geometry' not in config:
        geom = 'circle'
        gamma, edges, backward_edges = torus_knot(1, 0, config['n'], 3, 0)
    elif config['geometry'] == 'torus':
        p = config['p']
        q = config['q']
        n = config['n']
        r1 = config['r1']
        r2 = config['r2']
        geom = 'torus knot' + f' p={p}, q={q}, n={n}, r1={r1}, r2={r2}'
        gamma, edges, backward_edges = torus_knot(p, q, n, r1, r2)
    elif config['geometry'] == 'file':
        filename = config['filename']
        gamma, edges, backward_edges, components = read_obj(filename)
        geom = filename.split('/')[-1].split('.')[0]
    elif config['geometry'] == 'equilateral triangle':
        geom = 'equilateral triangle'
        gamma, edges, backward_edges = simple_equilateral_triangle()
    elif config['geometry'] == 'circle':
        geom = 'circle'
        gamma, edges, backward_edges = torus_knot(1, 0, config['n'], 3, 0)

    if components is None:
        components = [np.arange(len(gamma))]

    # dictionary indicating which component a vertex belongs to
    components_inv = {v:k for k, component in enumerate(components) for v in component}

    if 'name' not in config:
        config['name'] = geom

    if 'ignore edges' in config and config['ignore edges']:
        gamma, edges, backward_edges = create_edges(gamma)

    # check if loop is closed, if not close it
    edges = check_loop(edges)
    
    # transformations
    if 'flip' in config and config['flip']:
        gamma = invert_YZ(gamma) # blender uses Y as up, but we want Z as up
    if 'flatten' in config and config['flatten']:
        gamma = set_Z(gamma, 0)  # flatten the curve
    if 'decimate' in config and config['decimate'] < 1:
        gamma, edges = decimate(gamma, edges, config['decimate'])
        backward_edges = {v:k for k,v in edges.items()}
    
    # Plotting configurations
    if config['figure setup'] == '3d':
        # 4 3d subplots in a 2x2 grid
        angle = (30,30)
        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(221, projection='3d', xlabel='x', ylabel='y', zlabel='z')
        if 'angle' in config:
            angle = config['angle']
        ax1.view_init(*angle)
        ax2 = fig.add_subplot(222, projection='3d', xlabel='x', ylabel='y')
        ax2.view_init(90, 0)
        ax3 = fig.add_subplot(223, projection='3d', xlabel='x', ylabel='z')
        ax3.view_init(0, 90)
        ax4 = fig.add_subplot(224, projection='3d', xlabel='x', ylabel='y', zlabel='z')
        ax4.view_init(0, 0)
        axes = np.array([ax1,ax2,ax3,ax4])
        txt = fig.text(0.07, 0.97, '', fontsize=10, color='black')
    elif config['figure setup'] == '2d 2x1':
        # 2 subplots in a 1x2 grid, one in 3D and one in 2D
        angle = (35, 45, 0)
        if 'angle' in config:
            angle = config['angle']
        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(121, projection='3d', xlabel='x', ylabel='y', zlabel='z')
        ax1.view_init(*angle)
        ax2 = fig.add_subplot(122, xlabel='x', ylabel='y')
        axes = np.array([ax1,ax2])
        txt = fig.text(0.07, 0.97, '', fontsize=10, color='black')
    elif config['figure setup'] == '2d 1x1':
        # 1 subplot in 2D
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, xlabel='x', ylabel='y', projection='3d')
        ax.view_init(0, 90)
        ax.view_init(0,90,-90) # for hourglass
        axes = np.array([ax])
        txt = fig.text(0.07, 0.97, '', fontsize=10, color='black')
    elif config['figure setup'] == '3d 1x1':
        # 1 subplot in 3D
        angle = (30, 25)
        if 'angle' in config:
            angle = config['angle']
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, xlabel='x', ylabel='y', zlabel='z', projection='3d')
        ax.view_init(*angle)
        axes = np.array([ax])
        txt = fig.text(0.07, 0.97, '', fontsize=10, color='black')

    # plt.rcParams['figure.constrained_layout.use'] = True

    if 'max_edge_length_ratio' not in config:
        config['max_edge_length_ratio'] = 1

    print(f'Loaded {geom} – with {len(gamma)} vertices')

    if 'visuals' not in config:
        config['visuals'] = 'clean'

    if config['visuals'] == 'gradient computation':
        plot_gradient_computation()
        exit(0)

# main ------------------------------------------------------------------------------------
if __name__ == "__main__":

    # configs imported from config.py
    config = default_config
    # config = torus_4_1

    load_config(config)
    Mobius_gradient_descent()