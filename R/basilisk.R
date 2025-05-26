#TODO write this function correctly for creating the Cellbender env maybe use that function in cellbender call and if its not there create it
cellbender_env <- basilisk::BasiliskEnvironment(
  envname = "cellbender",
  pkgname = "DOtools",
  packages = c("python=3.7.16",
               "cellbender=0.3.2")
  )

#TODO write this function correctly for creating the remaining python function env maybe use that function in cellbender call and if its not there create it
dotools_env <- basilisk::BasiliskEnvironment(
  envname = "dotools",
  pkgname = "DOtools",
  packages = c("python=3.11.11",
               "celltypist=1.6.3",
               "scanpro=0.3.2")
  )
