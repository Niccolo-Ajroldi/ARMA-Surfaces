
# required packages
needed_packages <- c(
  "fda", 
  "rgl", 
  "animation", 
  "latex2exp", 
  "progress", 
  "freqdom", 
  "multiwave"
)

# new packages
new_packages  <- needed_packages[!(needed_packages %in%installed.packages()[,"Package"])] 

# install required packages
if(length(new_packages))
{
  install.packages(new_packages)
}