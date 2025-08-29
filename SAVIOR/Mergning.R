library(tidyverse)
library(igraph)

df <- read_csv("overlap_result_hit_125.csv")


# Part 1: Creating unique IDs for overlapping structures ------------------

# Create a two-column dataframe of connections: from 'Source_Identifier' to 'Overlap_hits'
edge_list <- df %>% 
  # only rows with overlaps
  filter(`Overlap count` > 0) %>%
  # Split the comma-separated Overlap_hits into separate rows
  separate_rows(`Overlap hits`, sep = ",") %>%
  # Keep only the two columns that define the connection
  select(from = `Source Identifier`, to = `Overlap hits`) %>%
  # Ensure there are no empty 'to' values
  filter(to != "")

# This graph represents all the overlap relationships
graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Find the connected components (i.e., the groups of identical structures)
components <- components(graph)
group_membership <- components$membership

# Creating a clean lookup table that maps each Source_Identifier to a numeric GroupID
group_lookup <- tibble(
  `Source Identifier` = names(group_membership),
  GroupID = group_membership
)


# Join the GroupIDs back to the original dataframe
df_grouped <- df %>%
  left_join(group_lookup, by =  'Source Identifier')

# Assign a unique GroupID to compounds that had no overlaps
max_group_id <- max(df_grouped$GroupID, na.rm = TRUE)
na_rows <- is.na(df_grouped$GroupID)
new_group_ids <- max_group_id + seq_len(sum(na_rows))

df_grouped$GroupID[na_rows] <- new_group_ids


# Split the dataframe into a list of dataframes based on the GroupID
# where each element contains all entries for that structure.
list_of_groups <- split(df_grouped, df_grouped$GroupID)
saveRDS(list_of_groups, file = "list_of_structural_groups.rds") # saving list



# saving in df
final_df <- bind_rows(list_of_groups)

final_df <- final_df %>%
  arrange(GroupID, `Source Identifier`)

write.csv(final_df, "unique_structural_groups.csv", row.names = FALSE)



# Part 2: Visualisation and creating overlap matrix -----------------------

library(UpSetR)
library(ggplot2)

# working with final_df from previous part
final_df <- read.csv('unique_structural_groups.csv')

# 1. Get the list of sources for each unique structural group (GroupID)
source_per_group <- final_df %>%
  group_by(GroupID) %>%
  summarise(Sources = list(unique(Source)), .groups = 'drop')


# 2. Prepare the data for ComplexUpset
upset_data <- source_per_group %>%
  unnest(Sources) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = GroupID,
    names_from = Sources,
    values_from = present,
    values_fill = 0
  )


# main upset plot commmand
upset(as.data.frame(upset_data), 
      sets = c("NCATS", "ChEMBL", "VTIs", "Manual_curated", "HTIs", "POSTera"),
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.7, 0.3),
      number.angles = 0,  # Horizontal text
      point.size = 4,
      line.size = 2.5,
      mainbar.y.label = "Number of Structures",
      sets.x.label = "Database Size",
      # text.scale: intersection size, intersection degree, sets labels, numbers above bars, main bar y, set size x
      text.scale = c(2, 1.8, 1.5, 2.5, 2, 1.2),  # Increased 4th element for bar text
      sets.bar.color = "#2E86AB",
      main.bar.color = "#A23B72",
      matrix.color = "#F18F01",
      shade.color = "grey90")


# extra: checking numbers in each column
upset_data %>% 
  filter( POSTera == TRUE & rowSums(.[,c("NCATS", "ChEMBL", "VTIs", "Manual_curated", "HTIs")]) == 0) %>% 
  nrow()
