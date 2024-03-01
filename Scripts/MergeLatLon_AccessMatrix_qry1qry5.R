
# Load necessary libraries
library(dplyr)

# Read the CSV files into dataframes
df_8 <- read.csv("~/Documents/local-git/Waterton/Input/8_qryWatertonCROSSTAB_firstvisitperlocation.csv")
df_1 <- read.csv("~/Documents/local-git/Waterton/Input/1_qryWatertonBASE.csv")

# Select only the necessary columns from df_1 (location, latitude, longitude)
df_1_selected <- df_1 %>% select(location, latitude, longitude)

# Remove duplicate rows based on location in df_1
df_1_unique <- df_1_selected %>% distinct(location, .keep_all = TRUE)

# Filter df_1_unique to only have locations present in df_8
df_1_filtered <- df_1_unique[df_1_unique$location %in% df_8$location, ]

# Merge the dataframes based on 'location' using a left join
merged_df <- left_join(df_8, df_1_filtered, by = "location")

# Rearrange columns to place latitude and longitude next to location
merged_df <- merged_df %>% select(names(df_8)[1:which(names(df_8) == "location")], latitude, longitude, everything())

# Write the merged dataframe to a new CSV file
write.csv(merged_df, "~/Documents/local-git/Waterton/Output/firstvisitperlocation.csv", row.names = FALSE)