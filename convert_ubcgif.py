# ## Imports
import argparse
import os
import pandas as pd
import numpy as np
import xarray as xr

# Convert data from UBCGIF to netCDF format
# Written by Joshua Simmons (DARE) 08-2022

# File format described here:
# https://gif.eos.ubc.ca/software/FAQ
# https://grav3d.readthedocs.io/en/latest/content/files/meshfile.html#
# https://grav3d.readthedocs.io/en/latest/content/files/modelfile.html

# Call 
# python convert_ubcgif.py mesh.msh model.den --output output.nc
# mesh (.msh file) and model (.den or .sus file) are required
# output (.nc file) is optional - by default output will be named according to the mesh file

def parse_delta_cell(x0, cell_num, txt):
    '''
    Parse the cell width specification into vertices to store.
    We will store the i,j vertices per (https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html)
    Input:
        x0: x coordinate of the first vertex
        cell_num: number of cells
        txt: text to parse
    Output:
        x_out: x coordinates of the vertices (i,j)
        cell_size: cell width per cell 
    '''
    cell_size = np.zeros((cell_num,))
    cnt = 0
    for cell_bit in txt.split():
        cell_rep = cell_bit.split('*')
        if cell_rep.__len__() == 1:
            cell_rep.insert(0, '1')
        cell_size[cnt:cnt+int(cell_rep[0])] = float(cell_rep[1])
        cnt += int(cell_rep[0])
    x_out = x0 + cell_size.cumsum()
    # now get i,j vertices
    x_out = np.insert(x_out[:-1],0,x0)
    return x_out, cell_size

if __name__ == "__main__":
    # setup arguement parsing
    parser = argparse.ArgumentParser(description='Convert UBCGIF mesh and model files to netCDF format')
    parser.add_argument('msh_file', metavar='msh_file', type=str,
                        help='Mesh file location')
    parser.add_argument('mod_file', metavar='mod_file', type=str,
                        help='Model file location')
    parser.add_argument('--output',  type=str, nargs='*', dest='out_file',
                        help='Output file location')

    # First get the args in
    args = parser.parse_args()

    # get out_file if not provided
    if args.out_file is None:
        out_file = args.msh_file.split('.')[0] + '.nc'
    else:
        out_file = args.out_file

    # parse the mesh file
    print('Parsing mesh file {}...'.format(args.msh_file))
    with open(args.msh_file, 'r') as f:
        lines = f.readlines()
        metadata = {}
        # Number of cells in the East/North/vertical direction.
        cell_count = lines[0]
        for elem, name in zip(cell_count.split(), ['NE','NN','NZ']):
            metadata[name] = int(elem)
        # Coordinates, in meters, of the southwest top corner
        coords_0 = lines[1]
        for elem, name in zip(coords_0.split(), ['E0','N0','Z0']):
            metadata[name] = float(elem)
        # nth cell width in the easting direction (ordered W to E), northing direction (ordered S to N), thickness (ordered top to bottom
        cells_E = lines[2]
        metadata['E'], metadata['Ewidth'] = parse_delta_cell(metadata['E0'], metadata['NE'], cells_E)
        metadata['Espec'] = cells_E
        cells_N = lines[3]
        metadata['N'], metadata['Nwidth'] = parse_delta_cell(metadata['N0'], metadata['NN'], cells_N)
        metadata['Nspec'] = cells_N
        cells_Z = lines[4]
        metadata['Z'], metadata['Zwidth'] = parse_delta_cell(metadata['Z0'], metadata['NZ'], cells_Z)
        metadata['Zspec'] = cells_Z

    print('Parsing model file {} (this may take a while)...'.format(args.mod_file))
    # now read in (pandas read_csv is the quickest) and convert to correct shape
    den = pd.read_csv(args.mod_file,header=None)
    den = den.values.reshape((metadata['NN'],metadata['NE'],metadata['NZ']),order='C')
    print('Overall shape: {}'.format(den.shape))
    print('Saving to {}'.format(out_file))
    # construct and save as netcdf
    xr_out = xr.Dataset(
        {'data': (['N', 'E', 'Z'], den)},
        coords={'N': metadata['N'], 'E': metadata['E'], 'Z': metadata['Z']},
        attrs={_: metadata[_] for _ in metadata if _ not in ['N','E','Z']}
    )
    xr_out.to_netcdf(out_file,mode='w')
    print('Done!')
