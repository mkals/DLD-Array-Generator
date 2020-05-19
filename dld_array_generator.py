#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 09:40:07 2020

@author: Morten Kals
"""

import ezdxf
import numpy as np


class DLD_Array:
    """ A package to generate DLD arrays in DXF format. All distances in micrometers, angles in degrees.

    Channel attributes:
        C_D = total channel length
        C_L = channel width
        C_N = number of channel bend pairs
        C_separation = distance between parralell channel segments

    Turn attributes:
        T_N = number of walls in turns
        T_W = thckness of walls in turns

    DLD-array unit cell attributes:
        D_D = center-to-center downstram pilar distance
        D_L = center-to-center lateral pilar distance
        D = pilar dimater
        theta = array angle

    IO port attributes:
        IO_D = diameter of inlets/outlets
        IO_F = filet radius for inlets/oulets
        IO_L = IO stem length

        O_LP = percentage width of left outlet
        O_S = separation of outlets
        O_BR = bend radius
        O_D = downstream distance for outlet before bend

    """
    ### Geometry definition

    C_D = 6e3 # um, total channel length
    C_L = 100 # um, channel width

    C_N = 1 # number of channel bend pairs

    C_d = C_D/(C_N*2) if C_N != 0 else C_D # um, length of each channel segment

    C_separation = 100 # um, channel separation

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


    for_simulation = False

    ###


    def make_modelspace():
        # Create a new DXF document.
        doc = ezdxf.new(dxfversion='R2010')

        doc.header['$MEASUREMENT'] = 1
        doc.header['$INSUNITS'] = 13

        # Create new table entries (layers, linetypes, text styles, ...).
        doc.layers.new('TEXTLAYER', dxfattribs={'color': 2})

        # DXF entities (LINE, TEXT, ...) reside in a layout (modelspace,
        # paperspace layout or block definition).
        msp = doc.modelspace()

        return doc, msp


    def make_dld_block(self):

        b = Block()

        posts_in_L = int(self.C_L / self.D_L)
        posts_in_D = int(self.C_d / self.D_D)
        displacement_per_D = np.tan(self.theta*np.pi/180) * self.D_D

        min_edge_points = [0]
        max_edge_points = min_edge_points.copy()

        for d in range(posts_in_D):
            for l in range(posts_in_L): # +2 for complete array

                displacement = d*displacement_per_D % self.D_L
                center = (l*self.D_L + displacement - self.D_L, (d+1/2)*self.D_D)
                radius = self.D/2

                x_min = 0
                x_max = self.C_L

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


        min_edge_points.append(self.C_d)
        max_edge_points.append(self.C_d)

        # make lines on left side
        for i in range(int(len(min_edge_points)/2)):
            b.append(Line(0, min_edge_points[2*i], 0, min_edge_points[2*i+1]))

        # make lines on right side
        for i in range(int(len(max_edge_points)/2)):
            b.append(Line(self.C_L, max_edge_points[2*i], self.C_L, max_edge_points[2*i+1]))

        return b


    def make_dld_corner(self):
        b = Block()

        r0 = self.C_separation/2
        r1 = r0 + self.C_L

        x = r1
        y = 0

        t0 = 0
        t1 = 180

        b.append(Arc(x, y, r0, t0, t1))
        b.append(Arc(x, y, r1, t0, t1))

        Rs = np.linspace(r0, r1, self.T_N+1, endpoint=False)[1:] # center radii
        Rl = np.concatenate((Rs+self.T_W/2, Rs-self.T_W/2)) # wall radii

        for r in Rl:
            b.append(Arc(x, y, r, t0, t1))

        for r in Rs:
            b.append(Arc(x+r, y, self.T_W/2, t1, t0))
            b.append(Arc(x-r, y, self.T_W/2, t1, t0))

        return b


    def make_io_port(self, W):

        b = Block()

        F = self.IO_F
        R = self.IO_D/2
        L = self.IO_L

        x = W/2
        y = self.IO_L + np.sqrt((F+R)**2-(F+W/2)**2)

        theta = np.arccos((F + W/2)/(F + R)) * 180/np.pi

        y_line = self.IO_L

        b.append(Line(0, 0, 0, y_line))
        b.append(Line(W, 0, W, y_line))

        if self.for_simulation:
            b.append(Line(0, y_line, W, y_line))
            return b

        b.append(Arc(x, y, R, -theta, 180+theta))

        b.append(Arc(-F, L, F, 0, theta))
        b.append(Arc(F+W, L, F, 180-theta, 180))

        return b

    def make_outlet(self, W):

        if self.for_simulation:

            angle = 10
            outlet = Block()

            outlet.append(Line(0, 0, W, 0))
            outlet.translate(O_BR, 0)
            outlet.rotate(angle)
            outlet.translate(-O_BR, -O_BR)

            outlet.append(Arc(-O_BR, -O_BR, O_BR, 0, angle))
            outlet.append(Arc(-O_BR, -O_BR, O_BR+W, 0, angle))

        else:
            outlet = self.make_io_port(W)

            outlet.rotate(90)

            outlet.translate(-self.O_S/2, 0)

            outlet.append(Line(-self.O_S/2, 0, -self.O_BR, 0))
            outlet.append(Line(-self.O_S/2, W, -self.O_BR, W))

            outlet.append(Arc(-self.O_BR, -self.O_BR, self.O_BR, 0, 90))
            outlet.append(Arc(-self.O_BR, -self.O_BR, self.O_BR+W, 0, 90))

        outlet.translate(0, self.O_BR)

        return outlet

    def make_outlets(self):

        outlet_small = self.make_outlet(self.C_L * self.O_LP)
        outlet_large = self.make_outlet(self.C_L * (1-self.O_LP))

        outlet_small.mirror_y(translating = False)
        outlet_small.rotate(180)
        outlet_small.translate(self.C_L, 0)

        b = Block()
        b.merge(outlet_small)
        b.merge(outlet_large)

        b.translate(0, self.O_D)
        b.append(Line(0, 0, 0, self.O_D))
        b.append(Line(self.C_L, 0, self.C_L, self.O_D))

        return b


    def generate(self, filename):

        doc, msp = DLD_Array.make_modelspace()

        # Generate patterns
        dld_array = self.make_dld_block()
        dld_array.to_modelspace(msp)


        if self.for_simulation:
            doc.saveas(filename + '.dxf')
            return


        if self.C_N != 0:
            corner_top = self.make_dld_corner()
            corner_top.translate(0, self.C_d)
            corner_top.to_modelspace(msp)

            corner_bottom = self.make_dld_corner()
            corner_bottom.mirror_y(translating = False)
            corner_bottom.translate(self.C_L+self.C_separation, 0)
            corner_bottom.to_modelspace(msp)

        # Repeat patterns

        for i in range(self.C_N*2):

            dld_array.rotate(180, True)
            dld_array.translate(self.C_L+self.C_separation, 0)
            dld_array.to_modelspace(msp)

            if i>1 and i%2 == 0:
                corner_top.translate(2*(self.C_L+self.C_separation), 0)
                corner_top.to_modelspace(msp)

            if i>1 and i%2 == 1:
                corner_bottom.translate(2*(self.C_L+self.C_separation), 0)
                corner_bottom.to_modelspace(msp)


        # Generate inlet

        inlet = self.make_io_port(self.C_L)
        inlet.mirror_y(translating = False)
        inlet.to_modelspace(msp)


        # Generate outlet

        outlets = self.make_outlets()
        outlets.translate(self.C_N*2*(self.C_L+self.C_separation), self.C_d)
        outlets.to_modelspace(msp)


        # Save DXF document.
        doc.saveas(filename + ".dxf")



### Geometry primitves

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
