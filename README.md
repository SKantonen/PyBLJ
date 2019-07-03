# PyBLJ
PyBLJ (Python Beta for Lennard Jones) Python program designed to parse .FChk files and produce electronic betas

# Running PyBLJ
The example methanol input, input.com, can be run using g09 to produce Test.FChk, which is the most important file for PyBLJ to parse. Once the calculation is complete, PyBLJ can be run inside the directory with Test.FChk + lebded (included in this repository) to begin producing radial grids to be used with Test.FChk to produce radial densities. PyBLJ then uses these to fit betas and ionization energies for each atom-in-molecule.
