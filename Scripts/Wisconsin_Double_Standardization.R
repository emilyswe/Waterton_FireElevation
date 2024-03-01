library(vegan)

# Read your cleaned community matrix from the CSV file
community_matrix <- read.csv("~/Documents/local-git/Waterton/Input/Cleaned_Matrix_V2.csv")

# Apply max scaling - divide each species count by its maximum value across all sites
community_matrix_max_scaled <- decostand(community_matrix[, -1], method = "max")

# Apply total scaling - divide each site's scaled species counts by the total scaled abundance for that site
community_matrix_wisconsin_scaled <- decostand(community_matrix_max_scaled, method = "total")

# Add the 'Location' column back to the scaled dataframe
community_matrix_wisconsin_scaled <- cbind(Location = community_matrix$Location, community_matrix_wisconsin_scaled)

# Save the Wisconsin scaled data to a new CSV file
write.csv(community_matrix_wisconsin_scaled, "~/Documents/local-git/Waterton/Output/Wisconsin_Scaled_Matrix.csv", row.names = FALSE)
