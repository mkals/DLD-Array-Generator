#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 08:37:10 2020

@author: mkals
"""

import numpy as np
import ezdxf

# Paper with variable name convention
# https://link.springer.com/article/10.1007/s40820-019-0308-7
# - D: Downstream
# - L: Lateral


# Sim termination
SIM = False
PRINT = False

SCALING_FACTOR = 100 if PRINT else 1

sim_rows = 10

# Channel definition

C_D = 6e3 # um, total channel length
C_L = 100 # um, channel width

C_N = 1 # number of channel bend pairs

C_d = C_D/(C_N*2) if C_N != 0 else C_D # um, length of each channel segment

C_separation = C_L # um, channel separation

T_N = 4 # number of walls in turns
T_W = 4 # um, thickness of walls in turn


# Unit cell definiton

D_D = 30 # um, center-to-center downstream pilar distance
D_L = 20 # um, center-to-center lateral pilar distance
D = 10 # um, pilar diameter

theta = 2 # deg, array period


# IO definition

IO_D = 800 # um, diameter of inlets/outlets
IO_F = 100 # um, filet radius for inlets/oulets
IO_L = 200 # um, IO stem length

O_LP = 0.3 # percentage width of left outlet
O_S = 4e3 # um, separation of outlets
O_BR = 400 # um, bend radius
O_D = 100 # um, downstream distance for outlet before bend

# SCALE

C_D = SCALING_FACTOR * C_D
C_L = SCALING_FACTOR * C_L
C_d = SCALING_FACTOR * C_d
C_separation = SCALING_FACTOR * C_separation
T_W = SCALING_FACTOR * T_W
D_D = SCALING_FACTOR * D_D
D_L = SCALING_FACTOR * D_L
D = SCALING_FACTOR * D
IO_D = SCALING_FACTOR * IO_D
IO_F = SCALING_FACTOR * IO_F
IO_L = SCALING_FACTOR * IO_L
O_S = SCALING_FACTOR * O_S
O_BR = SCALING_FACTOR * O_BR
O_D = SCALING_FACTOR * O_D


if PRINT:
    IO_D = C_L # um, diameter of inlets and outlets
    O_S = 5e4 # um, separation of outlets
    O_D = 5e3
    O_BR = 2e4 # um, bend radius


# L = x, D = y

posts_in_L = int(C_L / D_L)
posts_in_D = int(C_d / D_D)
displacement_per_D = np.tan(theta*np.pi/180) * D_D

x0, y0 = 0, D_D

# Create a new DXF document.
doc = ezdxf.new(dxfversion='R2010')

doc.header['$MEASUREMENT'] = 1
doc.header['$INSUNITS'] = 13

# Create new table entries (layers, linetypes, text styles, ...).
doc.layers.new('TEXTLAYER', dxfattribs={'color': 2})

# DXF entities (LINE, TEXT, ...) reside in a layout (modelspace,
# paperspace layout or block definition).
msp = doc.modelspace()

#msp.add_text(
#    'Test',
#    dxfattribs={
#        'layer': 'TEXTLAYER'
#    }).set_pos((0, 0.2), align='CENTER')

class Line:
    def __init__(self, x0, y0, x1, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1

    def translate(self, d_x, d_y):
        self.x0 = self.x0+d_x
        self.y0 = self.y0+d_y
        self.x1 = self.x1+d_x
        self.y1 = self.y1+d_y

    def to_modelspace(self, msp):
        msp.add_lwpolyline([(self.x0, self.y0), (self.x1, self.y1)])


    def min_y(self):
        return min([self.y0, self.y1])

    def max_y(self):
        return max([self.y0, self.y1])

    def min_x(self):
        return min([self.x0, self.x1])

    def max_x(self):
        return max([self.x0, self.x1])

    def invert_y(self):
        self.y0 = - self.y0
        self.y1 = - self.y1


    def rotate(self, ox, oy, theta):

        dx0 = self.x0 - ox
        dy0 = self.y0 - oy

        self.x0 = dx0*np.cos(theta) - dy0*np.sin(theta)
        self.y0 = dx0*np.sin(theta) + dy0*np.cos(theta)

        dx1 = self.x1 - ox
        dy1 = self.y1 - oy

        self.x1 = dx1*np.cos(theta) - dy1*np.sin(theta)
        self.y1 = dx1*np.sin(theta) + dy1*np.cos(theta)


class Circle:
    def __init__(self, x, y, r):
        self.x = x
        self.y = y
        self.r = r

    def translate(self, d_x, d_y):
        self.x = self.x+d_x
        self.y = self.y+d_y

    def to_modelspace(self, msp):
        msp.add_circle((self.x, self.y), radius=self.r)


    def min_y(self):
        return self.y - self.r

    def max_y(self):
        return self.y + self.r

    def min_x(self):
        return self.x - self.r

    def max_x(self):
        return self.x + self.r

    def invert_y(self):
        self.y = - self.y

    def rotate(self, ox, oy, theta):

        dx = self.x - ox
        dy = self.y - oy

        self.x = dx*np.cos(theta) - dy*np.sin(theta)
        self.y = dx*np.sin(theta) + dy*np.cos(theta)


class Arc:

    def __init__(self, x, y, r, start_angle, end_angle):

        self.x = x
        self.y = y
        self.r = r
        self.start_angle = start_angle
        self.end_angle = end_angle

    def translate(self, d_x, d_y):
        self.x = self.x+d_x
        self.y = self.y+d_y

    def to_modelspace(self, msp):
        msp.add_arc(
            center=(self.x, self.y),
            radius=self.r,
            start_angle=self.start_angle,
            end_angle=self.end_angle
        )

    def start(self):
        return (self.start_angle % 360) * np.pi/180

    def end(self):
        return (self.end_angle % 360) * np.pi/180

    def min_y(self):
        if np.cos(self.start()) > 0 and np.cos(self.end()) < 0:
            return self.y + self.r
        return self.y + self.r * min(np.sin(self.start()), np.sin(self.end()))
    def max_y(self):
        if np.cos(self.start()) < 0 and np.cos(self.end()) > 0:
            return self.y + self.r
        return self.y + self.r * max(np.sin(self.start()), np.sin(self.end()))

    def min_x(self):
        if np.sin(self.start()) > 0 and np.sin(self.end()) < 0:
            return self.x + self.r
        return self.x + self.r * min(np.cos(self.start()), np.cos(self.end()))


    def max_x(self):
        if np.sin(self.start()) < 0 and np.sin(self.end()) > 0:
            return self.x + self.r
        return self.x + self.r * max(np.cos(self.start()), np.cos(self.end()))


    def invert_y(self):
        self.y = - self.y
        old_start_angle = self.start_angle
        self.start_angle = -self.end_angle
        self.end_angle = -old_start_angle

    def rotate(self, ox, oy, theta):

        dx = self.x - ox
        dy = self.y - oy

        self.x = dx*np.cos(theta) - dy*np.sin(theta)
        self.y = dx*np.sin(theta) + dy*np.cos(theta)

        self.start_angle = self.start_angle + theta * 180/np.pi
        self.end_angle = self.end_angle + theta * 180/np.pi


class Block:

    def __init__(self):
        self.shapes = []
        self.origin = (0,0)
        self.rotation = 0

    def append(self, e):
        self.shapes.append(e)

    def translate(self, d_x, d_y):
        self.origin += (d_x, d_y)
        for e in self.shapes:
            e.translate(d_x, d_y)

    def mirror_y(self, translating = True):

        if translating:
            max_y = min([e.min_y() for e in self.shapes])
            min_y = max([e.max_y() for e in self.shapes])
            mean_y = (max_y+min_y)/2

        for e in self.shapes:
            e.invert_y()

        if translating:
            self.translate(0, mean_y*2)


    def to_modelspace(self, msp):
        for e in self.shapes:
            e.to_modelspace(msp)

    def merge(self, other):
        for shape in other.shapes:
            self.append(shape)

    def rotate(self, theta, translating = False):

        if translating:
            max_y = min([e.min_y() for e in self.shapes])
            min_y = max([e.max_y() for e in self.shapes])
            mean_y = (max_y+min_y)/2

            max_x = min([e.min_x() for e in self.shapes])
            min_x = max([e.max_x() for e in self.shapes])
            mean_x = (max_x+min_x)/2

        for s in self.shapes:
            s.rotate(self.origin[0], self.origin[1], theta*np.pi/180)

        if translating:
            self.translate(mean_x*2, mean_y*2)


def make_dld_block():

    b = Block()

    min_edge_points = [0]
    max_edge_points = min_edge_points.copy()

    row_count = sim_rows if SIM else posts_in_D

    for d in range(row_count):
        for l in range(posts_in_L): # +2 for complete array

            displacement = d*displacement_per_D % D_L
            center = (l*D_L + displacement - D_L, (d+1/2)*D_D)
            radius = D/2

            x_min = 0
            x_max = C_L

            if center[0] - radius < x_min:
                # only space for partial circle

                if center[0] + radius > x_min:

                    phi = np.arccos(( center[0] - x_min ) / radius)
                    the = 180 - phi * 180/np.pi

                    b.append(Arc(center[0], center[1], radius, -the, the))

                    d_y = np.sin(phi)*radius

                    min_edge_points.append(center[1]-d_y)
                    min_edge_points.append(center[1]+d_y)

            elif center[0] + radius > x_max:
                if center[0] - radius < x_max:

                    phi = np.arccos(( center[0] - x_max ) / radius)
                    the = 180 - phi * 180/np.pi

                    b.append(Arc(center[0], center[1], radius, the, -the))

                    d_y = np.sin(phi)*radius

                    max_edge_points.append(center[1]-d_y)
                    max_edge_points.append(center[1]+d_y)

            else:
                b.append(Circle(center[0], center[1], radius))


    min_edge_points.append(C_d)
    max_edge_points.append(C_d)

    # make lines on left side
    for i in range(int(len(min_edge_points)/2)):
        b.append(Line(0, min_edge_points[2*i], 0, min_edge_points[2*i+1]))

    # make lines on right side
    for i in range(int(len(max_edge_points)/2)):
        b.append(Line(C_L, max_edge_points[2*i], C_L, max_edge_points[2*i+1]))

    return b


def make_dld_corner():
    b = Block()

    r0 = C_separation/2
    r1 = r0 + C_L

    x = r1
    y = 0

    t0 = 0
    t1 = 180

    b.append(Arc(x, y, r0, t0, t1))
    b.append(Arc(x, y, r1, t0, t1))

    Rs = np.linspace(r0, r1, T_N+1, endpoint=False)[1:] # center radii
    Rl = np.concatenate((Rs+T_W/2, Rs-T_W/2)) # wall radii

    for r in Rl:
        b.append(Arc(x, y, r, t0, t1))

    for r in Rs:
        b.append(Arc(x+r, y, T_W/2, t1, t0))
        b.append(Arc(x-r, y, T_W/2, t1, t0))

    return b


def make_io_port(W):

    b = Block()

    F = IO_F
    R = IO_D/2
    L = IO_L

    x = W/2
    y = IO_L + np.sqrt((F+R)**2-(F+W/2)**2)

    theta = np.arccos((F + W/2)/(F + R)) * 180/np.pi

    y_line = IO_L

    b.append(Line(0, 0, 0, y_line))
    b.append(Line(W, 0, W, y_line))

    if SIM:
        b.append(Line(0, y_line, W, y_line))
        return b

    b.append(Arc(x, y, R, -theta, 180+theta))

    b.append(Arc(-F, L, F, 0, theta))
    b.append(Arc(F+W, L, F, 180-theta, 180))

    return b

def make_outlet(W):

    if SIM:

        angle = 10
        outlet = Block()

        outlet.append(Line(0, 0, W, 0))
        outlet.translate(O_BR, 0)
        outlet.rotate(angle)
        outlet.translate(-O_BR, -O_BR)

        outlet.append(Arc(-O_BR, -O_BR, O_BR, 0, angle))
        outlet.append(Arc(-O_BR, -O_BR, O_BR+W, 0, angle))

    else:
        outlet = make_io_port(W)

        outlet.rotate(90)

        outlet.translate(-O_S/2, 0)

        outlet.append(Line(-O_S/2, 0, -O_BR, 0))
        outlet.append(Line(-O_S/2, W, -O_BR, W))

        outlet.append(Arc(-O_BR, -O_BR, O_BR, 0, 90))
        outlet.append(Arc(-O_BR, -O_BR, O_BR+W, 0, 90))

    outlet.translate(0, O_BR)

    return outlet

def make_outlets():

    outlet_small = make_outlet(C_L * O_LP)
    outlet_large = make_outlet(C_L * (1-O_LP))

    outlet_small.mirror_y(translating = False)
    outlet_small.rotate(180)
    outlet_small.translate(C_L, 0)

    b = Block()
    b.merge(outlet_small)
    b.merge(outlet_large)

    b.translate(0, O_D)
    b.append(Line(0, 0, 0, O_D))
    b.append(Line(C_L, 0, C_L, O_D))

    return b




# Generate patterns

dld_array = make_dld_block()
dld_array.to_modelspace(msp)

if SIM:

    doc.saveas('sim.dxf')

if C_N != 0:
    corner_top = make_dld_corner()
    corner_top.translate(0,C_d)
    corner_top.to_modelspace(msp)

    corner_bottom = make_dld_corner()
    corner_bottom.mirror_y(translating = False)
    corner_bottom.translate(C_L+C_separation, 0)
    corner_bottom.to_modelspace(msp)


# Repeat patterns

for i in range(C_N*2):

    dld_array.rotate(180, True)
    dld_array.translate(C_L+C_separation, 0)
    dld_array.to_modelspace(msp)

    if i>1 and i%2 == 0:
        corner_top.translate(2*(C_L+C_separation), 0)
        corner_top.to_modelspace(msp)

    if i>1 and i%2 == 1:
        corner_bottom.translate(2*(C_L+C_separation), 0)
        corner_bottom.to_modelspace(msp)


# Generate inlet

inlet = make_io_port(C_L)
inlet.mirror_y(translating = False)
inlet.to_modelspace(msp)


# Generate outlet

outlets = make_outlets()
outlets.translate(C_N*2*(C_L+C_separation), C_d)
outlets.to_modelspace(msp)


# Save DXF document.
doc.saveas('../Python Output/DLD mask.dxf')
