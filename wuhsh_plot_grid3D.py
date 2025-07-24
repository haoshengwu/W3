import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pi = 3.14159265358979323846

def read_grid3D(filename):
    with open(filename, 'r') as file:
        # Read first line
        first_line = file.readline().strip()
        if ' ' in first_line:
            dims = list(map(int, first_line.split()))
        else:
            dims = [int(first_line[i:i+10]) for i in range(0, len(first_line), 10)]
        
        nr, n_pol, n_tor = dims[0], dims[1], dims[2]  # 重命名变量
        print(f"Grid dimensions: nr={nr}, np={n_pol}, nt={n_tor}")
        
        # Create array
        grid3d = np.zeros((n_tor, n_pol, nr, 3))  # 使用新变量名
        print(f"Array shape: {grid3d.shape}")
        
        for k in range(n_tor):  # 使用新变量名
            # Read phi value
            line = file.readline().strip()
            phi = float(line)
            grid3d[k, :, :, 2] = phi
            
            # Read R coordinates
            for i in range(n_pol):  # 使用新变量名
                for j in range(nr):
                    line = file.readline().strip()
                    value = float(line)
                    grid3d[k, i, j, 0] = value
            
            # Read Z coordinates  
            for i in range(n_pol):  # 使用新变量名
                for j in range(nr):
                    line = file.readline().strip()
                    value = float(line)
                    grid3d[k, i, j, 1] = value
    return grid3d

def plot_grid3d_tor_cut_section(grid_slice, ax, dim1_range=None, dim2_range=None):
    """
    Plot a 2D grid in 3D XYZ coordinate system
    
    Parameters:
    -----------
    grid_slice : numpy.ndarray
        3D array with shape (dim1, dim2, 3) where last dimension is [R, Z, phi]
    ax : matplotlib 3D axis
        3D axis to plot on
    dim1_range : tuple or None
        Range for first dimension (start, end). If None, use full range
    dim2_range : tuple or None  
        Range for second dimension (start, end). If None, use full range
    """
    
    # Validate input
    if len(grid_slice.shape) != 3 or grid_slice.shape[2] != 3:
        raise ValueError("grid_slice must be 3D array with shape (dim1, dim2, 3)")
    
    dim1_size, dim2_size = grid_slice.shape[0], grid_slice.shape[1]
    
    # Set default ranges if not provided
    if dim1_range is None:
        dim1_range = (0, dim1_size)
    if dim2_range is None:
        dim2_range = (0, dim2_size)
    
    # Validate and apply ranges
    dim1_start = max(0, dim1_range[0])
    dim1_end = min(dim1_size, dim1_range[1])
    dim2_start = max(0, dim2_range[0])
    dim2_end = min(dim2_size, dim2_range[1])
    
    # Extract the specified range
    grid_region = grid_slice[dim1_start:dim1_end, dim2_start:dim2_end, :]
    
    # Extract R, Z, phi coordinates
    R = grid_region[:, :, 0]
    Z = grid_region[:, :, 1] 
    phi = grid_region[:, :, 2]
    
    # Convert cylindrical to Cartesian coordinates
    X = R * np.cos(np.radians(phi))
    Y = R * np.sin(np.radians(phi))
    
    # Plot grid lines along first dimension (blue)
    for i in range(grid_region.shape[0]):
        ax.plot(X[i, :], Y[i, :], Z[i, :], 'b-', alpha=0.7, linewidth=1)
    
    # Plot grid lines along second dimension (red)
    for j in range(grid_region.shape[1]):
        ax.plot(X[:, j], Y[:, j], Z[:, j], 'b-', alpha=0.7, linewidth=1)

def plot_line(line, ax, color='b', alpha=0.7, linewidth=1):
    """
    Plot a single line in 3D space
    
    Parameters:
    -----------
    line : numpy.ndarray
        2D array with shape (n_points, 3) where columns are [R, Z, phi]
    ax : matplotlib 3D axis
    color : str
        Line color
    """
    R = line[:, 0]
    Z = line[:, 1] 
    phi = line[:, 2]
    X = R * np.cos(np.radians(phi))
    Y = R * np.sin(np.radians(phi))
    
    # Plot single line
    ax.plot(X, Y, Z, color=color, alpha=alpha, linewidth=linewidth)

def plot_grid3d_rad_cut_section(grid3d, ir, ax, close=False, title=None):
    grid2d = grid3d[:, :, ir, :]  #POL-TOR
    ntor=grid2d.shape[0]
    npol=grid2d.shape[1]
    print("ntor, npol",ntor, npol)
    for i in range(1,ntor):
        line = grid2d[0:i,i,:]
        plot_line(line, ax)
    for i in range(ntor, npol-ntor+1):
        line = grid2d[:,i,:]
        plot_line(line, ax)
    
    for i in range(npol-ntor+1,npol):
        j=i-(npol-ntor)
        line = grid2d[j:ntor,i,:]
        plot_line(line, ax)

def plot_grid3d_pol_line(grid3d, it1, it2, ax, close=False, title=None):
    ntol,npol,nrad,_ = grid3d.shape
    for i in range(it1,it2+1):
        line=grid3d[i,i:npol-ntor+i,nrad-1,:]
        plot_line(line, ax)

    return


def plot_rz_surface(filename, ax, phi_start, phi_end, phi_step):
    """
    Read RZ coordinates from file and plot as 3D surface
    
    Parameters:
    -----------
    filename : str
        Path to file containing RZ coordinates (two columns)
    ax : matplotlib 3D axis
        3D axis to plot on
    phi_start : float
        Starting phi angle in degrees
    phi_end : float
        Ending phi angle in degrees
    phi_step : float
        Phi angle step size in degrees
    """
    
    # Read RZ coordinates
    data = np.loadtxt(filename)
    R_line = data[:, 0]
    Z_line = data[:, 1]
    n_points = len(R_line)
    
    # Generate phi angles
    phi_angles = np.arange(phi_start, phi_end + phi_step, phi_step)
    n_phi = len(phi_angles)
    
    # Create meshgrid
    R_mesh = np.tile(R_line[:, np.newaxis], (1, n_phi))
    Z_mesh = np.tile(Z_line[:, np.newaxis], (1, n_phi))
    phi_mesh = np.tile(phi_angles[np.newaxis, :], (n_points, 1))
    
    # Convert to Cartesian coordinates
    X_mesh = R_mesh * np.cos(np.radians(phi_mesh))
    Y_mesh = R_mesh * np.sin(np.radians(phi_mesh))
    
    # Plot surface
    ax.plot_surface(X_mesh, Y_mesh, Z_mesh, 
                   color='lightgray', alpha=0.3, 
                   edgecolor='none')



fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

plot_rz_surface('./DEBUG_gzinfo_inner',ax, 0, 60, 5)
plot_rz_surface('./DEBUG_gzinfo_outer',ax, 0, 60, 5)

solgrid3d = read_grid3D('./grid3D.dat')
ntor, npol, nrad, _ = solgrid3d.shape

plane_start = solgrid3d[0, 0:(npol-ntor), :, :] 
plot_grid3d_tor_cut_section(plane_start, ax)

plane_end = solgrid3d[ntor-1, (ntor-1):npol-1, :, :] 
plot_grid3d_tor_cut_section(plane_end, ax)

plot_grid3d_rad_cut_section(solgrid3d, nrad-1, ax)

plot_grid3d_pol_line(solgrid3d, 1, ntor-2,ax=ax)


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Grid')

plt.show()