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
    # This is a friend object for PathLevelMap. This object keeps states
    # of corresponded cells on the particular layer. 5 states are defined now:
    # 0 - "air", free zone,
    # 1 - gray zone, a bound planned for tracing,
    # 2 - expanded, the material with the appropiate shift of the bound.
    # 3 - material. Don't touch it.
    def __init__( self, xmin, ymin, z, sample_interval, data ):
        self.m = None
        self.sampleInterval = sample_interval
        self.xmin = xmin
        self.ymin = ymin
        self.z = z
        self.tmp = None
        self.setContourMap( xmin, ymin, z, sample_interval, data )
            
    def setContourMap(self, xmin, ymin, z, sample_interval, data ):
        if not isinstance(data, numpy.ndarray):
            raise TypeError
        self.sampleInterval = sample_interval
        self.xmin = xmin
        self.ymin = ymin
        self.z = z
        R, C = data.shape
        if self.m is None or self.m.shape != (R + 2, C + 2):
            self.m = numpy.zeros((R + 2, C + 2), dtype = numpy.int8)
        
        if data.dtype.name == "bool":
            self.m[1:-1, 1:-1] = data * 3
        elif data.dtype.name == "int8":
            self.m[1:-1, 1:-1] = data
        else:
            self.m[:, :] = 0
            self.m[1:-1, 1:-1][data > z] = 3
        
    def highlightBorders(self, threshold = 0, new_state = 1):
        h_sum = self.m[1:-1, :-2] + self.m[1:-1, 2:]
        v_sum = self.m[2:, 1:-1] + self.m[:-2, 1:-1]
        if new_state > 1:
            numpy.maximum( self.m[1:-1, 1:-1], (h_sum + v_sum > threshold) * new_state,
                          out = self.m[1:-1, 1:-1] )
        else:
            numpy.maximum( self.m[1:-1, 1:-1], (h_sum + v_sum > threshold),
                          out = self.m[1:-1, 1:-1] )
            
    def bigShift(self):
        R, C = self.m.shape
        if self.tmp is None:
            self.tmp = numpy.empty((R + 7, C + 8), dtype = numpy.int8)
        self.tmp[:,:] = 0
        # bottom left corner of the sum is aligned to 3,3
        self.tmp[4:-5, 3:-7] = self.m[1:-1,1:-1]
        for j in range(0, R+1):
           self.tmp[j, 2:]  += self.tmp[j+1, 0:-2]
           self.tmp[j, :-2] += self.tmp[j+1, 2:]
           self.tmp[j, 3:]  += self.tmp[j+3, :-3]
           self.tmp[j, :-3] += self.tmp[j+3, 3:]
           self.tmp[j, 2:]  += self.tmp[j+5, :-2]
           self.tmp[j, :-2] += self.tmp[j+5, 2:]
           self.tmp[j, :]   += self.tmp[j+6, :]
        for j in range(0, R-1):   
           self.tmp[j, :] += self.tmp[j+1, :] + self.tmp[j+2, :]
           self.tmp[j, 0:-2] += self.tmp[j, 1:-1] + self.tmp[j, 2:]
        numpy.maximum( self.m[1:-1, 1:-1], (self.tmp[0:-9, 2:-8] > 0) * 2,
                      out = self.m[1:-1, 1:-1] )
        self.highlightBorders(new_state = 2)
        
    def shiftBorder(self, value):
        if value < 0:
            self.m[1:-1, 1:-1] = 3 - self.m[1:-1, 1:-1]
            self.shiftBorder( - value )
            self.m[1:-1, 1:-1] = numpy.minimum(3, 6 - self.m[1:-1, 1:-1] * 2)
        else:
            i = int(value / self.sampleInterval)
            while i > 5:
                self.bigShift()
                i -= 5
            while i > 0:
                self.highlightBorders(new_state = 2)
                i -= 1
                if (i % 2) != 0:
                    self.highlightBorders(threshold = 2,
                                          new_state = 2)
                    # Expand more in diagonal direction
        
    def or_not(self, a):
        numpy.maximum( self.m[1:-1, 1:-1], (a.m[1:-1, 1:-1] == 0) * 3,
                      out = self.m[1:-1, 1:-1] )

    def subtract(self, a):
        if not a is None:
            numpy.minimum(self.m[1:-1, 1:-1], (a.m[1:-1, 1:-1] == 0) * 3,
                          out = self.m[1:-1, 1:-1] )
        
    def getBorder(self, climb = False ):  # return list of (i,j,z) or empty list
        while True:
            index = numpy.where(self.m == 1)
            if len(index[0]) == 0:
                return []
            r, c = index[0][0], index[1][0]
            ccw = climb == (self.m[r+1, c+1] == 0)
            nodes = self._traceWaterLine(self.m, r, c, self.z, ccw )
            if self.m[r, c] == 1:
                nodes = self._traceWaterLine(self.m, r, c, self.z, 
                                            not ccw )[::-1] + nodes
            elif nodes[-1][0:2] == (r+1, c+1):
                nodes.append(nodes[0])
            
            self.m[r, c] = 0

            while len(nodes) > 2 and nodes[-1] == nodes[-2]:
                nodes.pop(-1)
                
            if len(nodes) > 2:
                break
        # adjust climb/conventional
        for i in range(1, len(nodes)):
            ca, ra = nodes[i-1][0:2]
            cb, rb = nodes[i][0:2]
            dmc = self.m[rb+1, cb+2] - self.m[rb+1, cb]
            # add 1 for each coordinate because in nodes r,c are shifted by -1
            dmr = self.m[rb+2, cb+1] - self.m[rb, cb+1]
            k = (cb - ca) * dmr - (rb - ra) * dmc
            if k != 0:
                if (k < 0) != climb:
                    return nodes[::-1]
                break
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

    def _traceWaterLine( self, m, r0, c0, z, ccw ):
        r, c = r0, c0
        R, C = m.shape
        ans = []
        if ccw:
            d = 0
            if m[r, c+1] == 1:         # horizontal direction expected
                ans.append((c, r, z))  # top right corner
        else:
            d = 2
            if m[r+1, c] == 1:         # vertical direction expected
                ans.append((c, r, z))
        dro, dco, drd, dcd = self.DELTAS[d]
        while True:
            t = 0
            while True:
                # check diagonal
                if m[r + drd, c + dcd] == 1:
                    ans.append((c + (dcd + 1)//2 - 1, 
                                r + (drd + 1)//2 - 1,
                                z))
                    if m[r + dro, c + dro] == 1:
                        m[r + dro, c + dro] = 0
                    r += drd
                    c += dcd
                    break
                # check horizontal or vertical direction
                elif m[r + dro, c + dco] == 1:
                    if t > 1:
                        ans.append((c + (2 * dco > dcd) - 1, 
                                    r + (2 * dro > drd) - 1,
                                    z))
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
                    # Look for continuation near the last turn
                    if len(ans) > 3:  # far enough from the start point
                        c = ans[-2][0] + 1
                        r = ans[-2][1] + 1
                        if dcd < 0:
                            c -= 1
                        if drd < 0:
                            r -= 1
                        if m[r,c] == 1:
                            break
                    return ans
                else:
                    if t == 0:
                        if dro == 0:
                            dd = drd * dco
                        else:
                            dd = -dcd * dro
                    t += 1
                    d = (d + dd) % 8
                    dro, dco, drd, dcd = self.DELTAS[d]
            m[r, c] = 0   # mark current cell
      

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
        #  Return (distance, index, None)
        #  or (distance, index_of_the_next_point, cross_point)
        #  Distance is in sampleInterval units
        best_i = (0)
        best_d = max(self.m.shape) ** 2 * 4
        best_p = None
        for i in range(0, len(nodes)):
            # distance to node
            d2 = (nodes[i][0] - p[0]) ** 2 + (nodes[i][1] - p[1]) ** 2
            if d2 < best_d:
                best_d = d2
                best_i = i
                best_p = None
            # distance to edge
            if i > 0: 
                d2, inside, pc = self.projectToLine( nodes[i-1], nodes[i], p ) 
                if inside and d2 < best_d:
                    best_d = d2
                    best_i = i - 1
                    best_p = pc
            if best_d == 0:
                break
        return (math.sqrt(best_d), best_i, best_p)
    
    
    def prependRamp(self, nodes, z0, slope, start_point=None):
        # prepend the given list of nodes by the ramp from the level z0.
        # Slope is the maximum value of abs(dz/dl) along the ramp.
        # start_point is used only for cyclic path. If this parameter given
        # the ramp will start at the point nearest to the start point.
        # Return a new list of nodes.
        if len(nodes) < 2 or slope <= 0 or nodes[0][2] >= z0:
            return nodes
        ramp = []
        # Try to locate cycle
        cycle = 1
        if nodes[0] == nodes[-1]:        
            # cyclic path, cycle points to the duplicate of the first node
            cycle = len(nodes) - 1
        while cycle < len(nodes):
            if nodes[0] == nodes[cycle]:
                break
            cycle += 1
        if cycle == len(nodes) - 1:   # all nodes form a one big cycle
            if not start_point is None:
                distance, index, point = self.nearestPoint( nodes, start_point )
                nodes = nodes[:-2]   # remove duplicated node from the tail
                                     # but not touch the original nodes
                # break line in this point
                if index > 0:
                    nodes = nodes[index:] + nodes[:index]
                if not point is None:
                    nodes.insert(0, point + nodes[0][2:])
            else:
                nodes = nodes[:-2]   # remove duplicated node from the tail
            
            # Now start from the first node.
            # Remember, that we should return to the first node.
            while True:
                nn = nodes[0]
                ramp.append(nn[0:2] + (z0,) + nn[3:]) # add point to the ramp
                nodes.append( nodes.pop(0))  # cyclic shift
                if nn[2] >= z0:
                    nodes.append(nodes[0])  # return to the first node
                    break
                dx = nodes[0][0]-nn[0]
                dy = nodes[0][1]-nn[1]
                distance = max(1, math.sqrt(dx*dx + dy*dy))
                if distance * slope > z0 - nn[2]:  # add a new node and terminate 
                    k = (z0 - nn[2]) / slope / distance
                    ramp.append((nn[0] + k * dx, nn[1] + k * dy) + nn[2:])
                    nodes.append(ramp[-1])
                    break
                else:
                    z0 -= distance * slope
            
            return ramp + nodes                                                        
                    
        else:  # linear path, probably with a cyclic fragment at the head. 
               # Go down along it aslant, or
               # make ramp moving forward and backward along path.
            di = 1
            i = 0   # A start point and a start direction
                    # The ramp will be constructed in reverse direction.
            if cycle < len(nodes):
                i = cycle
                di = -1
            else:
               cycle = None
            nn = nodes[i]
            for j in range(2, len(nodes)):  # Find cycle
                if nodes[j] == nn:
                    i = cycle = j
                    di = -1
            # now nodes[i] is a landing point
            z = nn[2]         
            while z < z0:
                # inspect the next point
                dx = nodes[i+di][0]-nn[0]
                dy = nodes[i+di][1]-nn[1]
                distance = max(1, math.sqrt(dx*dx + dy*dy))
                if distance * slope > z0 - z:  # make the final point
                    k = (z0 - z) / slope / distance
                    ramp.append((nn[0] + k * dx, nn[1] + k * dy, z0) + nn[3:])
                    break
                z += distance * slope
                i += di
                nn = nodes[i]                    
                ramp.append(nn[0:2] + (z,) + nn[3:])
                if i == 0:
                    if cycle:
                        i = cycle
                    else:
                        di = 1
                elif i >= len(nodes) - 1:
                    di = -1
            return ramp[::-1] + nodes                                                        
    
    def appendNodes(self, head, tail, max_distance ):
        # join two lists of nodes. On success head = head + junction + tail
        # If the tail is cyclic, it can be rotated so the shortest 
        # junction will be made.
        # Return true on success, false otherwise.
        max_distance /= self.sampleInterval
        index = 0
        point = None
        if tail[0] == tail[-1]:    # cyclic list
            distance, index, point = self.nearestPoint(tail, head[-1])
        else:
            distance = math.sqrt((head[-1][0]-tail[0][0]) ** 2 +
                                 (head[-1][1]-tail[0][1]) ** 2)

        if distance > max_distance:
            return False
            
        if point is None:
            if not self.freePass( head[-1], tail[0] ):
                return False
        else:
            if not self.freePass( head[-1], point ):
                return False
            head.append( point + tail[index][2:] )
        head.extend( tail[index:] )
        if index > 0:  #             Tail is cyclic
            head.extend( tail[1:index+1] )   # Tail.last == Tail.first, skip it
        return True
      
    def freePass(self, a, b):  
        # check if free stright pass by zones 0,1 or 2 from point a to point b exists 
        if a[0] > b[0]:
            a, b = b, a
        for i in range(int(math.ceil( a[0] )), int(math.ceil( b[0] ))):
            j = int(math.ceil(a[1] + (i - a[0]) * (b[1] - a[1]) / (b[0] - a[0])))
            if self.m[j+1, i+1] > 2:
                return False
        if a[1] > b[1]:
            a, b = b, a
        for j in range(int(math.ceil( a[1] )), int(math.ceil( b[1] ))):
            i = int(math.ceil(a[0] + (j - a[1]) * (b[0] - a[0]) / (b[1] - a[1])))
            if self.m[j+1, i+1] > 2:
                return False
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
        output.append(Path.Command("G1", {"Z": float(nodes[0][2]), "F": speed.vertFeed}))

        last = len(nodes) - 1
        # Cycle through each point of loop
        prev = nodes[0]
        for i in range(1, last + 1):
            this = nodes[i]
            if i < last:
                nxt = nodes[i + 1]
                if ((nxt[0] - this[0]) * (this[1] - prev[1]) ==
                    (nxt[1] - this[1]) * (this[0] - prev[0]) and
                    (nxt[0] - this[0]) * (this[0] - prev[0]) > 0 and 
                    (nxt[0] - this[0]) * (this[2] - prev[2]) ==
                    (nxt[2] - this[2]) * (this[0] - prev[0]) ):
                    #TODO check z here
                    continue
            output.append(
                Path.Command("G1", {"X": float(xmin + this[0] * si),
                                    "Y": float(ymin + this[1] * si),
                                    "Z": float(this[2]),
                                    "F": speed.horizFeed})
            )
            prev = this

        return output

    def copy(self, out = None):
        if not out is None:
            if out.m.shape != self.m.shape:
                out.m = self.m.copy()
            else:
                out.m[:,:] = self.m[:,:]
            out.xmin = self.xmin
            out.ymin = self.ymin
            out.z = self.z 
            out.sampleInterval = self.sampleInterval
            out.tmp = None
            return out
        return ContourMap(self.xmin, self.ymin, self.z, self.sampleInterval, 
                          self.m[1:-1, 1:-1] )
