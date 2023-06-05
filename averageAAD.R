# load required packages
library(dplyr)
library(readxl)

# read in the Excel file
df <- read_excel("C:/Table_2022Q1_CORRIDORS_W_AADT.xlsx", sheet = "2021Q1_corridors")

# group the data by LINK_ID and calculate mean of RASTERVALU
grouped <- df %>%
  group_by(LINK_ID) %>%
  summarize(AADT = mean(RASTERVALU))

# add AADT_YEAR column
grouped$AADT_YEAR <- "2022Q1"

# rename AADT column to RASTERVALU
result <- grouped %>%
  select(LINK_ID, AADT, AADT_YEAR) %>%
  rename(RASTERVALU = AADT)

# write result to Excel file
write.xlsx(result, "/path/Table_2022Q1_CORRIDORS_W_AADT_final1.xlsx", rownames = FALSE)
