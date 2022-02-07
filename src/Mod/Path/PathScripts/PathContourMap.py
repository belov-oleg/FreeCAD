# -*- coding: utf-8 -*-
# ***************************************************************************
# *   Copyright (c) 2022 Oleg Belov <obelov@audiology.ru>                   *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published yb the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

import Path
import math
import numpy

class ContourMap():
    def __init__( self, xmin, ymin, z, sample_interval, data ):
        if not isinstance(data, numpy.ndarray):
            raise TypeError

        self.sampleInterval = sample_interval
        self.xmin = xmin
        self.ymin = ymin
        self.z = z
        R, C = data.shape
        self.m = numpy.empty((R + 2, C + 2), dtype = numpy.int8)
        
        if data.dtype.name == "bool":
            self.m[1:-1, 1:-1] = data * 2
        elif data.dtype.name == "int8":
            self.m[1:-1, 1:-1] = data
        else:
            self.m[:, :] = 0
            self.m[1:-1, 1:-1][data > z] = 2
        
    def highlightBorders(self):
        neighbour_sum = self.m[1:-1, :-2] + self.m[1:-1, 2:]
        neighbour_sum += self.m[2:, 1:-1]
        neighbour_sum += self.m[:-2, 1:-1]
        numpy.maximum( self.m[1:-1, 1:-1], neighbour_sum > 0,
                      out = self.m[1:-1, 1:-1] )
        
    def shiftBorder(self, value):
        if value < 0:
            self.m[1:-1, 1:-1] = 2 - self.m[1:-1, 1:-1]
            self.shiftBorder( - value )
            self.m[1:-1, 1:-1] = 2 - self.m[1:-1, 1:-1]
        else:
            for i in range(0, int(value / self.sampleInterval)):
                self.highlightBorders()
            numpy.minimum( self.m[1:-1, 1:-1] * 2, 2,
                            out = self.m[1:-1, 1:-1] )
        
    def or_not(self, a):
        numpy.maximum( self.m[1:-1, 1:-1], (a.m[1:-1, 1:-1] == 0) * 2,
                      out = self.m[1:-1, 1:-1] )

    def subtract(self, a):
        numpy.minimum(self.m[1:-1, 1:-1], (a.m[1:-1, 1:-1] == 0) * 2,
                      out = self.m[1:-1, 1:-1] )
        
    def getBorder(self, climb = False ):  # return list of (i,j,z) or empty list
        index = numpy.where(self.m == 1)
        if len(index[0]) == 0:
            return []
        r, c = index[0][0], index[1][0]
        nodes = self._traceWaterLine(self.m, r, c, self.z, climb )
        if self.m[r, c] == 1:
            nodes = self._traceWaterLine(self.m, r, c, self.z, 
                                         not climb )[::-1] + nodes
        self.m[r, c] = 0
        if len(nodes) > 3:
            if nodes[0] == nodes[-2] and max(abs(nodes[0][0]-nodes[-1][0]),
                                             abs(nodes[0][1]-nodes[-1][1])) == 1:
                nodes.pop()
            elif nodes[1] == nodes[-1] and max(abs(nodes[0][0]-nodes[1][0]),
                                             abs(nodes[0][1]-nodes[1][1])) == 1:
                nodes.pop(0)
        return nodes
      
    #          ro  co  rd  cd
    DELTAS = (( 0,  1,  1,  1),
              ( 1,  0,  1,  1),
              ( 1,  0,  1, -1),
              ( 0, -1,  1, -1),
              ( 0, -1, -1, -1),
              (-1,  0, -1, -1),
              (-1,  0, -1,  1),
              ( 0,  1, -1,  1))

    def _traceWaterLine( self, m, r0, c0, z, climb ):
        r, c = r0, c0
        R, C = m.shape
        if climb:
            d = 4
        else:
            d = 0
        dro, dco, drd, dcd = self.DELTAS[d]
        ans = []
        while True:
            t = 0
            while True:
                # check diagonal
                if m[r + drd, c + dcd] == 1:
                    ans.append((c + (dcd + 1)//2 - 1, 
                                r + (drd + 1)//2 - 1,
                                z))
                    r += drd
                    c += dcd
                    break
                # check horizontal or vertical direction
                elif m[r + dro, c + dco] == 1:
                    r += dro
                    c += dco
                    break
                elif t == 7:
                    if r == r0 and c == c0 and ans != []:
                        ans.append(ans[0])
                    else:
                        ans.append((c + (dcd + 1)//2 - 1, 
                                    r + (drd + 1)//2 - 1,
                                    z))
                    m[r, c] = 0
                    return ans
                else:
                    t += 1
                    if climb:
                        d = (d - 1) % 8
                    else:
                        d = (d + 1) % 8
                    dro, dco, drd, dcd = self.DELTAS[d]
            m[r, c] = 0   # mark current cell
        m[r, c] = 0   # mark current cell
        return ans
      

    def projectToLine(self, p0, p1, p ):  
        # project p to p0..p1 in 2D, return (distance**2, strictly inside, point)
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        l  = math.sqrt(dx * dx + dy * dy)
        if l > 0:
            ux = dx / l
            uy = dy / l   # unity vector
            sp = ux * (p[0] - p0[0]) + uy * (p[1] - p0[1]) # scalar product
            cross = (p0[0] + ux * sp, p0[1] + uy * sp)
        else:
            sp = 0
            cross = p0
        return ((p[0] - cross[0]) ** 2 + (p[1] - cross[1]) ** 2, 
                sp > 0 and sp < l,
                cross)


    def nearestPoint(self, nodes, p):
        #  Find the nearest point to i, j. 
        #  Return (distance, (index) or (indices), point)
        best_i = (0)
        best_d = max(self.z.shape) ** 2 * 4
        for i in range(0, len(nodes)):
            # distance to node
            d2 = (nodes[i][0] - p[0]) ** 2 + (nodes[j][1] - p[1]) ** 2
            if d2 < best_d:
                best_d = d2
                best_i = (i,)
            # distance to edge
            d2, inside, pc = self.projectToLine( nodes[i-1], nodes[i], p ) 
            if d2 < best_d:
                best_d = d2
                best_i = (i - 1, i)
                best_p = pc
            if best_d == 0:
                break
        if len(best_i) == 1:
            return (best_d, best_i, nodes[best_i[0]])
        else:
            return (best_d, best_i, best_p)
    
    
    def prependRamp(self, nodes, z0, slope, start_point=None):
        # prepend the given list of nodes by the ramp from the level z0.
        # Slope is the maximum value of abs(dz/dl) along the ramp.
        # start_point is used only for cyclic path. If this parameter given
        # the ramp will start at the point nearest to the start point.
        # Return a new list of nodes.
        if len(nodes) < 2 or slope <= 0 or nodes[0][2] >= z0:
            return nodes
          
        ramp = []
        if nodes[0] == nodes[-1]:        # cyclic path
            if not start_point is None:
                distance, index, point = self.nearestPoint( nodes, start_point )
                
                nodes.pop(-1)
                if len(index) == 2:
                    # break line in this point
                    for i in range(0, index[1]):
                        nodes.append( nodes.pop(0))
                    nodes.insert(0, point)
                else:
                    for i in range(0, index[0]):
                        nodes.append( nodes.pop(0))
            # Now start from the first node
            while True:
                nn = nodes[0]
                nodes.append(nn)
                if nn[2] <= z0:
                    break
                dx = nodes[1][0]-nn[0]
                dy = nodes[1][1]-nn[1]
                distance = math.sqrt(dx*dx + dy*dy)
                if distance < 1:
                    nodes.pop(0)  # duplicate
                elif distance * slope > z0 - nn[2]:
                    k = (z0 - nn[2]) / slope / distance
                    ramp.append((nn[0] + k * dx, nn[1] + k * dy) + nn[2:])
                    break
                else:
                    nodes.pop(0)
                    z0 -= distance * slope
                    ramp.append(nn[0:2] + (z0,) + nn[3:])
                    
        else:  # linear path. Make ramp forward and backward along path
            di = 1
            zc = (nodes[0][2] + z0) / 2
            i = 0
            guard = int(2 * (z0 - zc) / slope)
            while z0 > nodes[0][2]:
                nn = nodes[i]
                dx = nodes[i+di][0]-nn[0]
                dy = nodes[i+di][1]-nn[1]
                distance = math.sqrt(dx*dx + dy*dy)
                if distance < 1:
                    guard -= 1
                    if guard < 0:
                        break
                    i += di
                elif z0 > zc and distance * slope > z0 - zc:  # make turn
                    k = (z0 - nn[2]) / slope / distance
                    ramp.append((nn[0] + k * dx, nn[1] + k * dy, zc) + nn[3:])
                    z0 = 2 * zc - z0
                    ramp.append(nn[0:2] + (z0,) + nn[3:])
                    di = -di
                else:
                    z0 = max(z0 - distance * slope, nodes[0][2])
                    i += di
                    ramp.append(nodes[i+di][:2] + (z0,) + nodes[i+di][3:])
                if i >= len(nodes) - 1:
                    di = -1
                elif i == 0:
                    di = 1
                
        return ramp + nodes    
                                                      
    
    def appendNodes(self, head, tail, max_distance ):
        # join two lists of nodes. On success head = head + junction + tail
        # If the tail is cyclic, it can be rotated so the shortest 
        # junction will be made.
        # Return true on success, false otherwise.
        ind = 0
        point = None
        if tail[0] == tail[-1]:    # cyclic list
            distance, index, point = self.nearest_point(tail, head[-1])
            if len(ind) == 1:
                point = None
            ind = index[-1]
        else:
            distance = math.sqrt((head[-1][0]-tail[0][0]) ** 2 +
                                 (head[-1][1]-tail[0][1]) ** 2)

        if distance > max_distance:
            return False
            
        #TODO Verify path
        
        if not Point is None:
            head.append( point + tail[ind][2:] )
        head.expand( tail[index[-1]:-1] )
        head.expand( tail[:index[-1]] )
        return True
                   

    def nodesToPath( self, nodes, clearanceHeight, speed ):
        # generate the path commands
        output = []
        if len(nodes) <= 1:
            return []

        xmin = self.xmin
        ymin = self.ymin
        si   = self.sampleInterval

        # Position cutter to begin loop
        output.append(
            Path.Command("G0", {"Z": clearanceHeight, "F": speed.vertRapid})
        )
        output.append(
            Path.Command("G0", {"X": xmin + nodes[0][0] * si,
                                "Y": ymin + nodes[0][1] * si,
                                "F": speed.horizRapid})
        )
        output.append(Path.Command("G1", {"Z": self.z, "F": speed.vertFeed}))

        last = len(nodes) - 1
        # Cycle through each point on loop
        prev = nodes[0]
        for i in range(1, last + 1):
            this = nodes[i]
            if i < last:
                nxt = nodes[i + 1]
                if ((nxt[0] - this[0]) * (this[1] - prev[1]) ==
                    (nxt[1] - this[1]) * (this[0] - prev[0]) and
                    (nxt[0] - this[0]) * (this[0] - prev[0]) > 0):
                    continue

            output.append(
                Path.Command("G1", {"X": xmin + this[0] * si,
                                    "Y": ymin + this[1] * si,
                                    "F": speed.horizFeed})
            )
            prev = this

        return output

