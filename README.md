# PyBLJ
PyBLJ (Python Beta for Lennard Jones) is a program designed to parse .FChk files and produce electronic betas, particularly useful for deriving Lennard-Jones parameters in the SDLJ force field paradigm.

# Generating Test.FChk
In order to run PyBLJ, a readable Test.FChk output file from Gaussian09 must be present in the working directory. Example Gaussian09 input file are included in this Github repository (methanol.com, ethanol.com, etc.)

# Running PyLBJ
Once the Test.FChk file has been produced, PyBLJ can be run. PyBLJ requires the included 'lebded' file in the working directory to build the radial grids. PyBLJ also requires Gaussian09 installed in your working environment, as it calls to Gaussian to generate density files specified by radial grids built using the 'lebded' file. Once your working directory has been confirmed to contain Test.FChk, lebded, and can run g09, the command to run PyLBJ is simply executing it in Python. 

# Output
After running PyLBJ, a list of the valence beta, ionization energy (kcal/mol), and partial charges will be generated for corresponding to each atom in the initial .com Gaussian09 file. Betas and ionization energies can be converted to LJ parameters via the SDLJ method found at this link:
