# PyBLJ
PyBLJ (Python Beta for Lennard Jones) is a program designed to parse .FChk files and produce electronic betas, particularly useful for deriving Lennard-Jones parameters in the SDLJ force field paradigm (citation forthcoming!)

# Generating Test.FChk
In order to run PyBLJ, a readable Test.FChk output file from Gaussian09 must be present in the working directory. Example Gaussian09 input file are included in this Github repository (methanol.com, ethanol.com, etc.)

# Running PyLBJ
Once the Test.FChk file has been produced, PyBLJ can be run. PyBLJ requires the included 'lebded' file in the working directory to build the radial grids.
