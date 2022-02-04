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
        for i in range(0, int(value / self.sampleInterval)):
            self.highlightBorders()
        numpy.minimum( self.m[1:-1, 1:-1] * 2, 2,
                        out = self.m[1:-1, 1:-1] )
        
    def subtractMap(self, a):
        numpy.maximum( self.m[1:-1, 1:-1], (a.m[1:-1, 1:-1] == 0) * 2,
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

