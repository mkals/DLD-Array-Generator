# DLD-Array Generator

A package to generate DLD arrays in DXF format. All distances in micrometers, angles in degrees.

The array is parameterized by attributes of the `DLD_Array` class. Initializing an instance sets up an array with default values. Each of the attrbutes (outlined below) can then be set as desired. When the configurtation is complete, call `generate(filename)` to output the array as in a `.dxf` file format.

See example usage in [here](blob/master/dld_array_generator_example.py). 

## Channel attributes:
    C_D = total channel length
    C_L = channel width
    C_N = number of channel bend pairs
    C_separation = distance between parralell channel segments

## Turn attributes:
    T_N = number of walls in turns
    T_W = thckness of walls in turns

## DLD-array unit cell attributes:
    D_D = center-to-center downstram pilar distance
    D_L = center-to-center lateral pilar distance
    D = pilar dimater
    theta = array angle

## IO port attributes:
    IO_D = diameter of inlets/outlets
    IO_F = filet radius for inlets/oulets
    IO_L = IO stem length

    O_LP = percentage width of left outlet
    O_S = separation of outlets
    O_BR = bend radius
    O_D = downstream distance for outlet before bend
