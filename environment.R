###### Environment Setup 
### Set working directory
# Dynamic user name
username <- Sys.getenv("USERNAME")

# Define path
input_path <- file.path("C:/Users",  username, "add file path")

# Set wd
setwd(input_path)
test 