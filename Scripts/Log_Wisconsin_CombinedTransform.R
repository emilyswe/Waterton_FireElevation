library(vegan)

# Read your cleaned community matrix from the CSV file
community_matrix <- read.csv("~/Documents/local-git/Waterton/Input/Cleaned_Matrix_V2.csv")

# Apply log transformation to the community data (adding 1 to avoid log(0))
community_matrix_log_transformed <- log1p(community_matrix[, -1])

# Perform max scaling on the log-transformed data
community_matrix_max_scaled <- decostand(community_matrix_log_transformed, method = "max")

# Perform total scaling on the max-scaled data
community_matrix_wisconsin_scaled <- decostand(community_matrix_max_scaled, method = "total")

# Add the 'Location' column back to the scaled dataframe
community_matrix_wisconsin_scaled <- cbind(Location = community_matrix$Location, community_matrix_wisconsin_scaled)

# Save the combined transformed data to a new CSV file
write.csv(community_matrix_wisconsin_scaled, "~/Documents/local-git/Waterton/Input/Combined_Log_Wisconsin_Scaled_Matrix.csv", row.names = FALSE)
