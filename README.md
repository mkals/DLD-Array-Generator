# DLD Array Generator

A package to generate DLD arrays in DXF format. All distances in micrometers, angles in degrees.

The array is parameterized by attributes of the `DLD_Array` class. Initializing an instance sets up an array with default values. Each of the attrbutes (outlined below) can then be set as desired. When the configurtation is complete, call `generate(filename)` to output the array as in a `.dxf` file format.

## Installation
`pip install dld-array-generator`

## Examples Usage
```
from dld_array_generator import DLD_Array

filename = "DLD example"

g = DLD_Array()

# Update channel desired channel parameters
g.C_D = 10e3 # channel length
g.C_L = 200 # channel width
g.C_N = 2 # bend pairs (producing 4 bends in total)

g.generate(filename)
```

## Attributes
### Channel attributes:
    C_D = total channel length
    C_L = channel width
    C_N = number of channel bend pairs
    C_separation = distance between parralell channel segments

### Turn attributes:
    T_N = number of walls in turns
    T_W = thckness of walls in turns

### DLD-array unit cell attributes:
    D_D = center-to-center downstram pilar distance
    D_L = center-to-center lateral pilar distance
    D = pilar dimater
    theta = array angle

### IO port attributes:
    IO_D = diameter of inlets/outlets
    IO_F = filet radius for inlets/oulets
    IO_L = IO stem length

    O_LP = percentage width of left outlet
    O_S = separation of outlets
    O_BR = bend radius
    O_D = downstream distance for outlet before bend
