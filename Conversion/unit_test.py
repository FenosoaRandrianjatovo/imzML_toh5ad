# Source https://pyimzml.readthedocs.io/en/latest/index.html

from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import anndata

# Specify the path of  imzML file with the corresponding ibd (they should be in the same folder)
imzml_file = '/home/fenosoa/scratch/Maya_Project/meta_base_MSI_data/pm.imzML'


# Create the parser object (this automatically handles the ibd file)
p = ImzMLParser(imzml_file)

# Initialize lists to store the intensity data and coordinates
spectra_list = []
coords_list = []

# Loop over all pixel coordinates using the index
for idx, coords in enumerate(p.coordinates):
    # Get m/z values and intensities for this pixel using its index
    mzs, intensities = p.getspectrum(idx)
    spectra_list.append(intensities)
    coords_list.append(coords)

# Assume the first pixel's m/z axis is common for all pixels
common_mzs, _ = p.getspectrum(0)

# Convert the list of spectra to a NumPy array
# Each row corresponds to a pixel and each column to an m/z value.
data_matrix = np.array(spectra_list)

# Create a DataFrame for spatial coordinates.
# Adjust the column names based on whether you have 2D or 3D data.
if len(coords_list[0]) == 3:
    coords_df = pd.DataFrame(coords_list, columns=['x', 'y', 'z'])
else:
    coords_df = pd.DataFrame(coords_list, columns=['x', 'y'])

# Create an AnnData object:
# - X holds the intensity data (pixels x m/z values)
# - obs contains spatial coordinates for each pixel.
adata = anndata.AnnData(X=data_matrix, obs=coords_df)

# Annotate the variable (feature) axis with the common m/z values.
adata.var['m/z'] = common_mzs

# Save the AnnData object to an h5ad file.
output_file = 'pm_output.h5ad'
adata.write(output_file)

print(f"Conversion complete. Data has been saved to {output_file}")

# This code is working
