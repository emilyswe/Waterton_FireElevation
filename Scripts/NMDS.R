
# You don't need to install vegan if it's already installed, just load the library
library(vegan)

# Set seed for reproducibility
set.seed(2)

# Import your community matrix from the CSV file you exported
community_matrix <- read.csv("~/Documents/local-git/Waterton/Output/Wisconsin_Scaled_Matrix.csv", row.names = 1)

# Run NMDS with your community matrix, specifying Bray-Curtis dissimilarity
# and setting the number of attempts to find a stable solution to 100
example_NMDS <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# Look at the Shepard plot of the NMDS
stressplot(example_NMDS)

plot(example_NMDS)