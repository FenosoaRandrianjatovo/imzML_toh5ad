# Install necessary packages if not already installed
if (!requireNamespace("Cardinal", quietly = TRUE)) {
  BiocManager::install("Cardinal")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter")
}

library(Cardinal)
library(SingleCellExperiment)
library(zellkonverter)

# Specify the path to your imzML file (ensure the corresponding ibd is in the same folder)
imzml_file <- '/home/fenosoa/scratch/Maya_Project/meta_base_MSI_data/pm.imzML'

# Read the imzML data; the ibd file is automatically located if it is in the same directory
ms_data <- readImzML(imzml_file)

# Extract the intensity matrix:
# - Rows correspond to pixels
# - Columns correspond to m/z values
intensity_matrix <- intensity(ms_data)

# Extract spatial coordinates for each pixel
coords <- coordinates(ms_data)
coords_df <- as.data.frame(coords)

# Extract the m/z values from the imzML data
mz_values <- mz(ms_data)

# For the SingleCellExperiment, rows should represent features (m/z values) and columns cells (pixels).
# Thus, transpose the intensity matrix.
intensity_matrix_t <- t(intensity_matrix)

# Create a SingleCellExperiment object:
# - The assay 'counts' holds the transposed intensity data.
# - colData contains the spatial coordinates.
sce <- SingleCellExperiment(
  assays = list(counts = intensity_matrix_t),
  colData = coords_df
)

# Annotate the feature data (rowData) with the m/z values.
rowData(sce)$mz <- mz_values

# Specify the output file name
output_file <- "output.h5ad"

# Write the SingleCellExperiment to a single h5ad file
writeH5AD(sce, output_file)

cat("Conversion complete. Data saved to", output_file, "\n")
