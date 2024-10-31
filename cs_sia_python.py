#!/usr/bin/env python
# coding: utf-8



import numpy as np
import os
from scipy.sparse.linalg import cgs
from scipy.sparse.linalg import cgs, LinearOperator
os.chdir('/')


# Dimension of a cubed sphere (cs) face
cx = 32


def cs_to_vec(cs_dof, vec_cs):
    # Flatten both arrays into 1D arrays
    vec_unordered = vec_cs.flatten()
    dof_unordered = cs_dof.flatten()
    
    # Sort the dof_unordered and get the indices
    I = np.argsort(dof_unordered)
    
    # Reorder vec_unordered based on the sorted indices
    vec = vec_unordered[I]
    
    return vec

def vec_to_cs(vec):
    cs = vec.reshape((32, 6, 32))
    return cs


def read_cs_bin(fnam, kx=None, prec='real*4', cx=510):
    """
    Function to read in cube sphere binary output

    Parameters:
    fnam : str
        Input path and file name.
    kx : list or int, optional
        Vertical indices to read (default is [1]).
    prec : str, optional
        Numeric precision (default is 'real*4').
    cx : int, optional
        Dimension of the cubed sphere face (default is 510).

    Returns:
    fld : np.array
        Output array of specified dimensions.
    """
    # Default values
    if kx is None:
        kx = [1]  # Default to [1] if not provided

    # Determine the number of bytes per element based on precision
    prec_dict = {
        'int8': 1,
        'int16': 2,
        'uint16': 2,
        'int32': 4,
        'uint32': 4,
        'single': 4,
        'real*4': 4,
        'float32': 4,
        'int64': 8,
        'uint64': 8,
        'double': 8,
        'real*8': 8,
        'float64': 8,
    }

    if prec not in prec_dict:
        raise ValueError(f"Unsupported precision: {prec}")

    preclength = prec_dict[prec]
    dtype = '>f8' if preclength == 8 else '>f4'  # Big-endian float64 or float32
    fld = np.zeros((cx, 6, cx, len(kx)), dtype=np.float64 if preclength == 8 else np.float32)

    with open(fnam, 'rb') as fid:
        for k in range(len(kx)):
            # Adjust kx for 0-based indexing in Python
            skip = (kx[k] - 1) * cx  * 6 * cx
            fid.seek(skip * preclength)  # Move to the correct position in the file
            
            # Read data and reshape
            data = fid.read(cx * 6 * cx * preclength)
            if len(data) < cx * 6 * cx * preclength:
                raise ValueError('Reached end of file or insufficient data read.')

            # Unpack the data according to precision, specifying column-major order
            fld[:, :, :, k] = np.frombuffer(data, dtype=dtype).reshape((cx, 6, cx), order='F')

    return fld

# Set folder path and read the geometric arrays
subfold = 'GridMask_Laure'
land_i = read_cs_bin(f'{subfold}/mask_32x6x32.bin', [1], 'float64', cx=32)
topo_i = read_cs_bin(f'{subfold}/topo_malai_32x6x32.bin', [1], 'float64', cx=32)
RAC_i = read_cs_bin(f'{subfold}/RAC.data', [1], 'float64', cx=32)
DXC_i = read_cs_bin(f'{subfold}/DXC.data', [1], 'float64', cx=32)
DYC_i = read_cs_bin(f'{subfold}/DYC.data', [1], 'float64', cx=32)
DXG_i = read_cs_bin(f'{subfold}/DXG.data', [1], 'float64', cx=32)
DYG_i = read_cs_bin(f'{subfold}/DYG.data', [1], 'float64', cx=32)
XC_i = read_cs_bin(f'{subfold}/XC.data', [1], 'float64', cx=32)
YC_i = read_cs_bin(f'{subfold}/YC.data', [1], 'float64', cx=32)

#Convert it into the right format

land = np.zeros((32, 32, 6), dtype=int)
for k in range(len(land_i[0, 0, :])):
    for j in range(len(land_i[0, :, 0])):
        for i in range(len(land_i[:, 0, 0])):
            land[k, i, j] = land_i[i, j, k]

topo = np.zeros((32, 32, 6), dtype=int)
for k in range(len(topo_i[0, 0, :])):
    for j in range(len(topo_i[0, :, 0])):
        for i in range(len(topo_i[:, 0, 0])):
            topo[k, i, j] = topo_i[i, j, k]

RAC = np.zeros((32, 32, 6), dtype=int)
for k in range(len(RAC_i[0, 0, :])):
    for j in range(len(RAC_i[0, :, 0])):
        for i in range(len(RAC_i[:, 0, 0])):
            RAC[k, i, j] = RAC_i[i, j, k]

DXC = np.zeros((32, 32, 6), dtype=int)
for k in range(len(DXC_i[0, 0, :])):
    for j in range(len(DXC_i[0, :, 0])):
        for i in range(len(DXC_i[:, 0, 0])):
            DXC[k, i, j] = DXC_i[i, j, k]

DYC = np.zeros((32, 32, 6), dtype=int)
for k in range(len(DYC_i[0, 0, :])):
    for j in range(len(DYC_i[0, :, 0])):
        for i in range(len(DYC_i[:, 0, 0])):
            DYC[k, i, j] = DYC_i[i, j, k]

DXG = np.zeros((32, 32, 6), dtype=int)
for k in range(len(DXG_i[0, 0, :])):
    for j in range(len(DXG_i[0, :, 0])):
        for i in range(len(DXG_i[:, 0, 0])):
            DXG[k, i, j] = DXG_i[i, j, k]

DYG = np.zeros((32, 32, 6), dtype=int)
for k in range(len(DYG_i[0, 0, :])):
    for j in range(len(DYG_i[0, :, 0])):
        for i in range(len(DYG_i[:, 0, 0])):
            DYG[k, i, j] = DYG_i[i, j, k]

XC = np.zeros((32, 32, 6), dtype=int)
for k in range(len(XC_i[0, 0, :])):
    for j in range(len(XC_i[0, :, 0])):
        for i in range(len(XC_i[:, 0, 0])):
            XC[k, i, j] = XC_i[i, j, k]

YC = np.zeros((32, 32, 6), dtype=int)
for k in range(len(YC_i[0, 0, :])):
    for j in range(len(YC_i[0, :, 0])):
        for i in range(len(YC_i[:, 0, 0])):
            YC[k, i, j] = YC_i[i, j, k]

# Physical parameters

# Glen's law parameter: Pa^3 per year
# Corresponds to soft ice
Aglen = 5e-15

# Exponent in Glen's Law
nglen = 3

# Density x gravity
rhog = 917 * 9.8

# Compute Csia
Csia = 2 * Aglen / (nglen + 2) * (rhog) ** nglen

# Model time step
dt = 50

# Topology file path
topol_file = 'topol.txt'

# Open the file for reading and parse its content
obj_array = []

class DataRow:
    def __init__(self, field1, field2, field3, field4):
        self.Field1 = field1
        self.Field2 = field2
        self.Field3 = field3
        self.Field4 = field4

# Open the file for reading
with open(topol_file, 'r') as fid:
    for line in fid:
        tokens = line.split()
        if len(tokens) >= 7:
            # Create a new DataRow object and assign fields
            new_row = DataRow(tokens[0], tokens[2], tokens[4], tokens[6])
            # Append the new object to the array
            obj_array.append(new_row)

# Dictionary mapping edge names to numbers
direction_dict = {
    'N.Edge': 1,
    'S.Edge': 2,
    'E.Edge': 3,
    'W.Edge': 4
}

facets = np.zeros((2, 6, 4), dtype=int)

# Populate the facets array based on obj_array
for i in range(24):
    # Convert the string fields from obj_array to integers
    facet = int(obj_array[i].Field2)
    neighbor_facet = int(obj_array[i].Field4)
    edge = direction_dict[obj_array[i].Field1]
    neighbor_edge = direction_dict[obj_array[i].Field3]
    
    # Update the facets array with neighbor information
    facets[0, facet - 1, edge - 1] = neighbor_facet
    facets[1, facet - 1, edge - 1] = neighbor_edge

# Generate list of unique cell IDs (Degrees of Freedom)
dofs_i = np.arange(1, cx * cx * 6 + 1).reshape((32, 6, 32))

dofs = []

for i in range(len(dofs_i)):
    dofs_t = np.transpose(dofs_i[i], (1, 0))
    dofs.append(dofs_t)

dofs = np.array(dofs)

num_dofs = cx * cx * 6
# dof_matrix stores the relationships: dof ID, north neighbor, south neighbor, east neighbor, west neighbor
dof_matrix = np.zeros((num_dofs, 5), dtype=int)
dof_matrix[:, 0] = np.arange(1, num_dofs + 1)

# Initialize the doffaces array
doffaces_i = np.zeros((6, 32, 32), dtype=int)
for i in range(6):
    doffaces_i[i, :, :] = dofs[:, :, i]

doffaces = []

for i in range(len(doffaces_i)):
    d = np.transpose(doffaces_i[i], (1, 0))
    doffaces.append(d)

doffaces = np.array(doffaces)

# Loop through each face and compute the neighbors
for i in range(6):
    east_dofs = np.zeros_like(doffaces[i, :, :])
    west_dofs = np.zeros_like(doffaces[i, :, :])
    north_dofs = np.zeros_like(doffaces[i, :, :])
    south_dofs = np.zeros_like(doffaces[i, :, :])

    # Assign internal edges
    east_dofs[:, :-1] = doffaces[i, :, 1:]
    west_dofs[:, 1:] = doffaces[i, :, :-1]
    north_dofs[:-1, :] = doffaces[i, 1:, :]
    south_dofs[1:, :] = doffaces[i, :-1, :]

    # Handle boundary edges using neighbor facet information
    # North edge
    neighbor = facets[0, i, 0] - 1
    edge = facets[1, i, 0]
    if neighbor >= 0:  # Ensure neighbor index is valid
        if edge == 1:
            north_dofs[-1, :] = doffaces[neighbor, -1, ::-1]
        elif edge == 2:
            north_dofs[-1, :] = doffaces[neighbor, 0, :]
        elif edge == 3:
            north_dofs[-1, :] = doffaces[neighbor, :, -1]
        elif edge == 4:
            north_dofs[-1, :] = doffaces[neighbor, ::-1, 0]

    # South edge
    neighbor = facets[0, i, 1] - 1
    edge = facets[1, i, 1]
    if neighbor >= 0:  # Ensure neighbor index is valid
        if edge == 1:
            south_dofs[0, :] = doffaces[neighbor, -1, :]
        elif edge == 2:
            south_dofs[0, :] = doffaces[neighbor, 0, ::-1]
        elif edge == 3:
            south_dofs[0, :] = doffaces[neighbor, ::-1, -1]
        elif edge == 4:
            south_dofs[0, :] = doffaces[neighbor, :, 0]

    # East edge
    neighbor = facets[0, i, 2] - 1
    edge = facets[1, i, 2]
    if neighbor >= 0:  # Ensure neighbor index is valid
        if edge == 1:
            east_dofs[:, -1] = doffaces[neighbor, -1, :]
        elif edge == 2:
            east_dofs[:, -1] = doffaces[neighbor, 0, ::-1]
        elif edge == 3:
            east_dofs[:, -1] = doffaces[neighbor, ::-1, -1]
        elif edge == 4:
            east_dofs[:, -1] = doffaces[neighbor, :, 0]

    # West edge
    neighbor = facets[0, i, 3] - 1
    edge = facets[1, i, 3]
    if neighbor >= 0:  # Ensure neighbor index is valid
        if edge == 1:
            west_dofs[:, 0] = doffaces[neighbor, -1, ::-1]
        elif edge == 2:
            west_dofs[:, 0] = doffaces[neighbor, 0, :]
        elif edge == 3:
            west_dofs[:, 0] = doffaces[neighbor, :, -1]
        elif edge == 4:
            west_dofs[:, 0] = doffaces[neighbor, ::-1, 0]

    # Update dof_matrix with neighbor information
    current_dofs = doffaces[i, :, :]
    dof_matrix[current_dofs.flatten() - 1, 1] = north_dofs.flatten()
    dof_matrix[current_dofs.flatten() - 1, 2] = south_dofs.flatten()
    dof_matrix[current_dofs.flatten() - 1, 3] = east_dofs.flatten()
    dof_matrix[current_dofs.flatten() - 1, 4] = west_dofs.flatten()

# Convert cube sphere fields to vectors using cs_to_vec
dxc = cs_to_vec(dofs, DXC)
dyc = cs_to_vec(dofs, DYC)
dxg = cs_to_vec(dofs, DXG)
dyg = cs_to_vec(dofs, DYG)
rac = cs_to_vec(dofs, RAC)
land_vec = cs_to_vec(dofs, land)
xc = cs_to_vec(dofs, XC)
yc = cs_to_vec(dofs, YC)
topo_vec = cs_to_vec(dofs, topo)

# Calculate flow to neighboring cells based on land cells
flow2north = (land_vec == 1) & (land_vec[dof_matrix[:, 1] - 1] == 1)
flow2south = (land_vec == 1) & (land_vec[dof_matrix[:, 2] - 1] == 1)
flow2east = (land_vec == 1) & (land_vec[dof_matrix[:, 3] - 1] == 1)
flow2west = (land_vec == 1) & (land_vec[dof_matrix[:, 4] - 1] == 1)

# Initialize thickness (h_vec)
h_vec = 10 * np.ones(len(dof_matrix))
h_vec[land_vec == 0] = 0

# Save initial thickness
h_vec0 = h_vec.copy()

def cgfunc(x, dofs, A_columns):
    """
    Implements matrix action for a subset of the degrees of freedom.

    Parameters:
    x : array_like
        Input vector.
    dofs : array_like
        Degrees of freedom.
    A_columns : array_like
        Matrix action columns.

    Returns:
    y : array_like
        Result of the matrix-vector product.
    """
    # Ensure indices are within bounds of x
    valid_indices = (dofs >= 0) & (dofs <= len(x))
    if not np.all(valid_indices):
        raise IndexError("Some indices in 'dofs' are out of bounds for the given vector.")
    
    x_vals = x[dofs-1]  # Get corresponding values from x based on dofs
    y = np.sum(A_columns * x_vals, axis=1)
    return y


def smb(surf, lat):
    """
    Calculates surface mass balance and equilibrium line altitude.

    Parameters:
    surf : ndarray
        Surface elevation.
    lat : ndarray
        Latitude.

    Returns:
    a : ndarray
        Surface mass balance.
    ela : ndarray
        Equilibrium line altitude.
    """
    ela = np.zeros_like(lat)
    a = np.zeros_like(lat)
    
    ela[lat < -30] = 5500 + 100 * lat[lat < -30]
    ela[lat > 30] = 5500 - 100 * lat[lat > 30]
    ela[(lat >= -30) & (lat <= 30)] = 2500
    
    elev_diff = surf - ela
    a[elev_diff > 0] = 1 * elev_diff[elev_diff > 0]
    a[elev_diff < 0] = 2 * elev_diff[elev_diff < 0]
    
    a[a > 5] = 5
    a[a < -20] = -20
    
    return a, ela

def matvec(x):
    return cgfunc(x, dof_matrix, Acols)

# Time iteration loop
for it in range(1000):
    # Begin building linear system
    surf_vec = topo_vec + h_vec

    # East neighbor minus west
    dsdx_c = (surf_vec[dof_matrix[:, 3] - 1] - surf_vec[dof_matrix[:, 4] - 1]) / (
        dxc[dof_matrix[:, 0] - 1] + dxc[dof_matrix[:, 4] - 1])

    # North neighbor minus south
    dsdy_c = (surf_vec[dof_matrix[:, 1] - 1] - surf_vec[dof_matrix[:, 2] - 1]) / (
        dyc[dof_matrix[:, 0] - 1] + dyc[dof_matrix[:, 2] - 1])

    # Diffusivity calculation
    D_vec = Csia * (dsdx_c**2 + dsdy_c**2 + 1e-8)**((nglen - 1) / 2) * h_vec**(nglen + 2)

    # Diffusivity at south and west faces
    Dx_vec = 0.5 * (D_vec[dof_matrix[:, 0] - 1] + D_vec[dof_matrix[:, 4] - 1])
    Dy_vec = 0.5 * (D_vec[dof_matrix[:, 0] - 1] + D_vec[dof_matrix[:, 2] - 1])

    Dx_west = Dx_vec
    Dx_east = Dx_vec[dof_matrix[:, 3] - 1]
    Dy_south = Dy_vec
    Dy_north = Dy_vec[dof_matrix[:, 1] - 1]

    Brhs = h_vec.copy()

    # Construct matrix components
    Iland = land_vec == 1
    I = Iland
    Amat00 = np.ones_like(h_vec)
    Amat00[Iland] = Amat00[Iland] + dt / rac[I] * (
        +dyg[I] / dxc[dof_matrix[I, 4] - 1] * Dx_west[I]
        + dyg[dof_matrix[I, 3] - 1] / dxc[I] * Dx_east[I]
        + dxg[I] / dyc[dof_matrix[I, 2] - 1] * Dy_south[I]
        + dxg[dof_matrix[I, 1] - 1] / dyc[I] * Dy_north[I]
    )

    Brhs[I] += dt / rac[I] * (
        -dyg[I] / dxc[dof_matrix[I, 4] - 1] * Dx_west[I] * (topo_vec[I] - topo_vec[dof_matrix[I, 4] - 1])
        + dyg[dof_matrix[I, 3] - 1] / dxc[I] * Dx_east[I] * (topo_vec[dof_matrix[I, 3] - 1] - topo_vec[I])
        - dxg[I] / dyc[dof_matrix[I, 2] - 1] * Dy_south[I] * (topo_vec[I] - topo_vec[dof_matrix[I, 2] - 1])
        + dxg[dof_matrix[I, 1] - 1] / dyc[I] * Dy_north[I] * (topo_vec[dof_matrix[I, 1] - 1] - topo_vec[I])
    )

    # Surface mass balance
    a_smb, ela = smb(surf_vec[I], yc[I])
    Brhs[I] += dt * a_smb

    # Construct matrix action components
    AmatNorth = np.zeros_like(h_vec)
    I = flow2north
    AmatNorth[I] -= dt / rac[I] * dxg[dof_matrix[I, 1] - 1] / dyc[I] * Dy_north[I]

    AmatSouth = np.zeros_like(h_vec)
    I = flow2south
    AmatSouth[I] -= dt / rac[I] * dxg[I] / dyc[dof_matrix[I, 2] - 1] * Dy_south[I]

    AmatEast = np.zeros_like(h_vec)
    I = flow2east
    AmatEast[I] -= dt / rac[I] * dyg[dof_matrix[I, 3] - 1] / dxc[I] * Dx_east[I]

    AmatWest = np.zeros_like(h_vec)
    I = flow2west
    AmatWest[I] -= dt / rac[I] * dyg[I] / dxc[dof_matrix[I, 4] - 1] * Dx_west[I]

    # Combine matrix components
    Acols_i = [Amat00, AmatNorth, AmatSouth, AmatEast, AmatWest]
    Acols = np.transpose(Acols_i, (1,0))
    n = len(h_vec)

    A = LinearOperator((n, n), matvec=matvec)
    # Conjugate gradient solver
    h_vecnew, cginfo = cgs(A, Brhs)
    h_vecnew[h_vecnew < 0] = 0

    # Check for convergence information
    if cginfo < 0:
        print(f"CGS failed to converge, flag: {cginfo}")
    elif cginfo == 0:
        print("CGS converged successfully.")
    else:
        print(f"CGS converged in {cginfo} iterations.")

    # Update thickness
    h_vec = h_vecnew

    
H = vec_to_cs(h_vec)
H0 = vec_to_cs(h_vec0)  

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def crossmap(TS, cx=None):
    """
    Plot a cube-sphere field as a cross after rotating and reflecting each face.
    
    Parameters:
    TS (numpy.ndarray): Input cube sphere field of shape (n, 6, n).
    cx (list or tuple): Two-element list or tuple [cmin, cmax] to set scaling for color axis.
    """
    # Default color axis scaling if not provided
    if cx is None:
        cx = [np.min(TS), np.max(TS)]
    
    # Extract each face and apply rotation and reflection
    # Extract each face and apply rotation and reflection
    a1 = np.fliplr(np.rot90(TS[:, 0, :], 3))
    a2 = np.fliplr(np.rot90(TS[:, 1, :], 3))
    a3 = np.fliplr(np.rot90(TS[:, 2, :], 3))  # Rotate a3 by 180°
    a4 = np.fliplr(np.rot90(TS[:, 3, :], 1))
    a5 = np.fliplr(np.rot90(TS[:, 4, :], 3))
    a6 = np.fliplr(np.rot90(TS[:, 5, :], 1))
    
    tmp = np.vstack((a2, np.rot90(a4, -1), np.rot90(a5, -1), a1))

    # Plot the main cross section
    fig, ax_main = plt.subplots(figsize=(8, 6))
    ax_main.set_position([0.05, 0.35, 0.8, 0.3])  # Adjusted position for colorbar space
    h1 = ax_main.pcolormesh(tmp.T, shading='flat', cmap='bone')
    ax_main.set_xticks([])
    ax_main.set_yticks([])

    # Color limits for all plots
    h1.set_clim(cx)

    # Add top and bottom faces
    pos = ax_main.get_position()
    
    # Top face (a3 rotated)
    ax_top = fig.add_axes([pos.x0 + pos.width / 4 * 2, pos.y0 + pos.height, pos.width / 4, pos.height])
    h_top = ax_top.pcolormesh(np.rot90(a3, 2).T, shading='flat', cmap='bone')
    h_top.set_clim(cx)
    ax_top.set_xticks([])
    ax_top.set_yticks([])

    # Bottom face (a6 rotated)
    ax_bottom = fig.add_axes([pos.x0 + pos.width / 4 * 2, pos.y0 - pos.height, pos.width / 4, pos.height])
    h_bottom = ax_bottom.pcolormesh(np.rot90(a6).T, shading='flat', cmap='bone')
    h_bottom.set_clim(cx)
    ax_bottom.set_xticks([])
    ax_bottom.set_yticks([])
    
crossmap(H)
