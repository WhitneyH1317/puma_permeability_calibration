# requirements.R

# List all required packages here
required_packages <- c(
  "amt",
  "tidyverse",
  "survival",
  "lubridate",
  "corrplot",
  "sf",
  "mapview",
  "raster",
  "lwgeom",
  "terra"
)

# Install any packages that are missing
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}

# Load all required packages
lapply(required_packages, library, character.only = TRUE)

message("All packages loaded.")