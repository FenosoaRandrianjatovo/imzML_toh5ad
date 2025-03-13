from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import anndata

# Specify the path to your imzML file (the corresponding ibd should be in the same folder)
imzml_file = '/home/fenosoa/scratch/Maya_Project/meta_base_MSI_data/pm.imzML'

# Create the parser object (this automatically handles the ibd file)
parser = ImzMLParser(imzml_file)

# Retrieve a common m/z axis from the first pixel
first_coords = parser.coordinates[0]
common_mzs, _ = parser.getspectrum(first_coords)

# Initialize lists to store intensity arrays and spatial coordinates
all_intensities = []
all_coords = []

# Loop over all pixel coordinates in the imzML file
for coords in parser.coordinates:
    mzs, intensities = parser.getspectrum(coords)
    # If the m/z axis for each spectrum might differ,
    # you would need to re-bin intensities to a common m/z axis.
    all_intensities.append(intensities)
    all_coords.append(coords)

# Convert the list of intensity arrays to a NumPy array.
# The resulting matrix will have shape (n_pixels, n_mz_values)
data_matrix = np.array(all_intensities)

# Create a DataFrame for spatial coordinates.
# imzML coordinates can be (x, y) or (x, y, z). Adjust accordingly.
if len(all_coords[0]) == 2:
    coords_df = pd.DataFrame(all_coords, columns=['x', 'y'])
elif len(all_coords[0]) == 3:
    coords_df = pd.DataFrame(all_coords, columns=['x', 'y', 'z'])
else:
    coords_df = pd.DataFrame(all_coords, columns=[f'coord_{i}' for i in range(len(all_coords[0]))])

# Create an AnnData object.
# X holds the intensity data, and the obs DataFrame holds spatial coordinates.
adata = anndata.AnnData(X=data_matrix, obs=coords_df)

# Add the common m/z axis as a variable annotation
adata.var['m/z'] = common_mzs

# Save the AnnData object to a single h5ad file.
output_file = 'output.h5ad'
adata.write(output_file)

print(f"Conversion complete. The data has been saved to {output_file}")
