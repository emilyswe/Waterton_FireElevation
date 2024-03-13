
# You don't need to install vegan if it's already installed, just load the library
library(vegan)

# Set seed for reproducibility
set.seed(2)
# Import your community matrix from the CSV file
community_matrix_raw <- read.csv("Input/04_Cleaned_Waterton.csv")

# Need to extract only species data columns 
# Assuming species data is from column 6 to 116. 
community_matrix <- community_matrix_raw[, 6:116] # Adjust X to the index for column DL

# Retain the site information for plotting. Assuming it's in the fourth column, adjust if it's different.
site_info <- community_matrix_raw[, 4]

# Convert any non-numeric columns in the community matrix to numeric 
community_matrix <- data.frame(lapply(community_matrix, as.numeric))

#save matrix
write_csv(community_matrix, "Input/05_Cleaned_Matrix.csv")

# Import your community matrix from the CSV file you exported
community_matrix <- read.csv("Input/05_Cleaned_Matrix.csv")

# Run NMDS with your community matrix, specifying Bray-Curtis dissimilarity
# and setting the number of attempts to find a stable solution to 100
example_NMDS <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# Look at the Shepard plot of the NMDS
stressplot(example_NMDS)

plot(example_NMDS)