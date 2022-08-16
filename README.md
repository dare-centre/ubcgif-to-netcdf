# Convert data from UBCGIF to netCDF format
Written by Joshua Simmons (DARE) 08-2022

A python script for converting UBCGIF files to netCDF format.

## Example call 
`python convert_ubcgif.py mesh.msh model.den --output output.nc`

Inputs:
    - `mesh.msh`: mesh file
    - `model.den`: density file
    - `output`: output netCDF file (e.g., output.nc)
Note:
mesh (.msh file) and model (.den or .sus file) are required. output (.nc file) is optional - by default output will be named according to the mesh file.

There is a provided notebook with example for plotting:

- `example_load_netcdf.ipynb`: An example of loading the converted netCDF file with xarray and plotting for a slice.

Output netcdf variables:
- `data`: model data (NN,NE,NZ) per [the GRAV3D manual](https://grav3d.readthedocs.io/en/latest/content/files/modelfile.html)

Output netcdf coordinates - each refers to a location of bottom-left corner of mesh cell per [matplotlib plotting spec](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolor.html):
- `N`: cell vertex in North direction
- `E`: cell vertex in East direction
- `Z`: cell vertex in vertical direction

Output netcdf metadata contains:
- `NE`: number of cells in East direction
- `NN`: number of cells in North direction
- `NZ`: number of cells in vertical direction
- `EO`: Coordinates, in meters, of the southwest top corner in Easting
- `NO`: Coordinates, in meters, of the southwest top corner in Northing
- `ZO`: Coordinates, in meters, of the southwest top corner in vertical
- `Espec`: Easting spacing of mesh cells as per mesh file specification
- `Nspec`: Northing spacing of mesh cells as per mesh file specification
- `Zspec`: Vertical spacing of mesh cells as per mesh file specification
  
## Requirements
- xarray

## File format
File format described here:

- https://gif.eos.ubc.ca/software/FAQ
- https://grav3d.readthedocs.io/en/latest/content/files/meshfile.html#
- https://grav3d.readthedocs.io/en/latest/content/files/modelfile.html

## To Do
- [ ] Not tested currently for complex grid specification in the mesh file - proceed with caution and please raise an issue if you encounter problems.

## License
Distributed under the GNU General Public License v3.