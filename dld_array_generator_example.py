from dld_array_generator import DLD_Array

filename = "DLD example"

g = DLD_Array()

# Update channel desired channel parameters
g.C_D = 10e3 # channel length
g.C_L = 200 # channel width
g.C_N = 2 # bend pairs (producing 4 bends in total)

g.generate(filename)
