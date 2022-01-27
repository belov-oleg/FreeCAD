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

import math
import numpy

# This is a square grid of elevations in given direction.
# Each cell [i,j] holds a mxaimum value of the model elevation in a square region
#   (xmin + i * sample_interval .. xmin + (i + 1) * sample_interval,
#    ymin + i * sample_interval .. ymin + (i + 1) * sample_interval)

class LevelMap():
    def __init__( self, xmin, xmax, ymin, ymax, zmin, sample_interval, border ):
        self.sampleInterval = sample_interval
        if xmax < xmin:
          xmax, xmin = xmin, xmax
        if ymax < ymin:
          ymax, ymin = ymin, ymax
        self.xmin = xmin
        cols = int(math.ceil((xmax - xmin) / sample_interval) + 1)
        self.xmax = self.xmin + cols * sample_interval
        self.ymin = ymin
        rows = int(math.ceil((ymax - ymin) / sample_interval) + 1)
        self.ymax = self.ymin + rows * sample_interval
        self.zmin = zmin
        self.matrix = None
        self.border = border
        self.z = numpy.full((rows + 2 * border, cols + 2 * border),  # elevation data
                             zmin - 1.0,
                             dtype = numpy.single
                            )
        self.kk = numpy.zeros(max(rows, cols) + 2 * border, dtype=int)
        
    def set_rotation(self):
        pass
        #TODO set matrix in accordance with sample_interval and given rotation.
        
    def reset(self):
        self.z[:] = self.zmin
      
    def rows(self):
        R, C = self.z.shape
        return R - 2 * self.border
      
    def columns(self):
        R, C = self.z.shape
        return C - 2 * self.border
      
    def level(self):
        return self.z[self.border:-self.border, self.border:-self.border]
      
    def add_facet( self, va, vb, vc ):
        #TODO if not self.matrix is Null: apply matrix
        xva = (va[0] - self.xmin) / self.sampleInterval + self.border
        xvb = (vb[0] - self.xmin) / self.sampleInterval + self.border
        xvc = (vc[0] - self.xmin) / self.sampleInterval + self.border
        yva = (va[1] - self.ymin) / self.sampleInterval + self.border
        yvb = (vb[1] - self.ymin) / self.sampleInterval + self.border
        yvc = (vc[1] - self.ymin) / self.sampleInterval + self.border
        
        b1 = self._add_edge( xva, yva, va[2], xvb, yvb, vb[2] )
        b2 = self._add_edge( xvb, yvb, vb[2], xvc, yvc, vc[2] )
        b3 = self._add_edge( xvc, yvc, vc[2], xva, yva, va[2] )
        if b1 or b2 or b3:
            self._add_triangle( xva, yva, va[2], xvb, yvb, vb[2], xvc, yvc, vc[2] )   
      
    def _add_edge( self, xa, ya, za, xb, yb, zb ):
        # This algorithm should be coded in C
        R, C = self.z.shape
        if ((xa < 0 and xb < 0) or (xa >= C and xb >= C) or 
            (ya < 0 and yb < 0) or (ya >= R and yb >= R)):
            return False
        
        if zb < za:
            xa, xb = xb, xa
            ya, yb = yb, ya
            za, zb = zb, za
        
        ia = int(math.floor( xa ))
        ja = int(math.floor( ya ))
        ib = int(math.floor( xb ))
        jb = int(math.floor( yb ))

        _x = xa - ia
        _y = ya - ja
        nx = abs(ia - ib)
        ny = abs(ja - jb)

        if nx > 0:
            dy_adx = (yb - ya) / abs(xb - xa)  # dy / abs(dx)
            dz_adx = (zb - za) / abs(xb - xa)
            if xb > xa:
                dx = 1
                y0 = ya + (1 - _x) * dy_adx  # y and z of the first crossection with 
                z0 = za + (1 - _x) * dz_adx  # a horizontal line of the grid
                if ia < 0:
                    y0 -= ia * dy_adx   # skip lines
                    z0 -= ia * dz_adx
                    nx += ia
                    ia = 0
                if ib >= C:
                    nx -= (ib - C + 1)
            else:
                dx = -1
                y0 = ya + _x * dy_adx  # y and z of the first crossection with 
                z0 = za + _x * dz_adx  # a horizontal line of the grid
                if ia >= C:
                    y0 += (ia - C + 1) * dy_adx   # skip lines
                    z0 += (ia - C + 1) * dz_adx
                    nx -= (ia - C + 1)
                    ia = C - 1
                if ib < 0:
                    nx += ib
              
            if y0 < 0 or ya < 0 or yb < 0 or y0 >= R or ya >= R or yb >= R: 
                # loop with range check
                for i in range(ia, ia + nx * dx, dx):
                    j = int(math.floor(y0))
                    if j >= 0 and j < R:
                        if self.z[j, i] < z0:
                            self.z[j, i] = z0
                    y0 += dy_adx
                    z0 += dz_adx
                    
            elif nx < 10:             # simple loop
                for i in range(ia, ia + nx * dx, dx):
                    j = int(math.floor(y0))
                    if self.z[j, i] < z0:
                        self.z[j, i] = z0
                    y0 += dy_adx
                    z0 += dz_adx
                    
            else:                     # vectorized, not necessary in C implementation
                if abs(dy_adx) > 0.0001:
                    self.kk[0:nx] = numpy.arange(y0, (y0 + dy_adx * (nx - 0.5)), 
                                                  dy_adx)
                    self.kk[0:nx] *= C
                else:
                    self.kk[0:nx] = int(y0) * C
                self.kk[0:nx] += numpy.arange(ia, ia + dx * nx, dx)
                if dz_adx < 0.001 / nx:
                    self.z.flat[self.kk[0:nx]] = numpy.maximum(
                                       self.z.flat[self.kk[0:nx]], z0)
                else:
                    self.z.flat[self.kk[0:nx]] = numpy.maximum(
                                  self.z.flat[self.kk[0:nx]], 
                                  numpy.arange(z0, z0 + dz_adx * (nx - 0.5), dz_adx))
                    
        if ny > 0:
            dx_ady = (xb - xa) / abs(yb - ya)  # dx / abs(dy)
            dz_ady = (zb - za) / abs(yb - ya)
            if yb > ya:
                dy = 1
                x0 = xa + (1 - _y) * dx_ady  # x and z of the first crossection with 
                z0 = za + (1 - _y) * dz_ady  # a horizontal line of the grid
                if ja < 0:
                    x0 -= ja * dx_ady   # skip lines
                    z0 -= ja * dz_ady
                    ny += ja
                    ja = 0
                if jb >= R:
                    ny -= (jb - R + 1)
            else:
                dy = -1
                x0 = xa + _y * dx_ady  # x and z of the first crossection with 
                z0 = za + _y * dz_ady  # a horizontal line of the grid
                if ja >= R:
                    x0 += (ja - R + 1) * dx_ady   # skip lines
                    z0 += (ja - R + 1) * dz_ady
                    ny -= (ja - R + 1)
                    ja = R - 1
                if jb < 0:
                    ny += jb
              
            if x0 < 0 or xa < 0 or xb < 0 or x0 >= C or xa >= C or xb >= C: 
                # loop with range check
                for j in range(ja, ja + ny * dy, dy):
                    i = int(math.floor(x0))
                    if i >= 0 and i < C:
                        if self.z[j, i] < z0:
                            self.z[j, i] = z0
                    x0 += dx_ady
                    z0 += dz_ady
                    
            elif ny < 10:             # simple loop        
                for j in range(ja, ja + ny * dy, dy):
                    i = int(math.floor(x0))
                    if self.z[j, i] < z0:
                        self.z[j, i] = z0
                    x0 += dx_ady
                    z0 += dz_ady
                    
            else:                     # vectorized, not necessary in C implementation
                base = ja * C + x0
                step = C * dy + dx_ady
                self.kk[0:ny] = numpy.arange(base, 
                                             base + step * (ny - 0.5), 
                                             step)
                if dz_ady < 0.001 / ny:
                    self.z.flat[self.kk[0:ny]] = numpy.maximum(
                                               self.z.flat[self.kk[0:ny]], z0)
                else:
                    self.z.flat[self.kk[0:ny]] = numpy.maximum(
                                  self.z.flat[self.kk[0:ny]], 
                                  numpy.arange(z0, z0 + dz_ady * (ny - 0.5), dz_ady))
                    

        if ib >=0 and jb >=0 and ib < C and jb < R and self.z[jb, ib] < zb:
            self.z[jb, ib] = zb
            
        return True    

    def _add_triangle( self, xa, ya, za, xb, yb, zb, xc, yc, zc ):
        # This algorithm should be coded in C
        # Only cells located inside the given triangle and not crossed with it's 
        # edges are marked.
        R, C = self.z.shape
        i0 = max(0, int(numpy.floor(min( xa, xb, xc ))) + 1)
        i1 = min(C - 1, int(numpy.floor(max( xa, xb, xc ))))
        j0 = max(0, int(numpy.floor(min( ya, yb, yc ))) + 1)
        j1 = min(R - 1, int(numpy.floor(max( ya, yb, yc ))))
        ni = i1 - i0
        nj = j1 - j0
        if ni <= 0 or nj <= 0:
            return
        
        if ni > nj:
            # sort ya <= yc <= yb
            if yb < ya:
                xb, xa = xa, xb
                yb, ya = ya, yb
                zb, za = za, zb
            if yc < ya:
                xc, xa = xa, xc
                yc, ya = ya, yc
                zc, za = za, zc
            if yc > yb:
                xc, xb = xb, xc
                yc, yb = yb, yc
                zc, zb = zb, zc
                
            jc = min(R - 1, int(numpy.floor(yc)))

        else:
            # sort xa <= xc <= xb
            if xb < xa:
                xb, xa = xa, xb
                yb, ya = ya, yb
                zb, za = za, zb
            if xc < xa:
                xc, xa = xa, xc
                yc, ya = ya, yc
                zc, za = za, zc
            if xc > xb:
                xc, xb = xb, xc
                yc, yb = yb, yc
                zc, zb = zb, zc
              
            ic = min(C - 1, int(numpy.floor(xc)))

        yab = ya - yb
        ybc = yb - yc
        yca = yc - ya
        xab = xa - xb
        xbc = xb - xc
        xca = xc - xa
        zac = za - zc
        zbc = zb - zc
            
        den =  xbc * yca - ybc * xca
        if abs(den) < 1.0:
            return

        dzdx = (zac * ybc + zbc * yca) / den
        dzdy = -(zac * xbc + zbc * xca) / den
        d    = 0

        if ni > nj:

            zfirst = za + (j0 - ya + int(dzdy > 0)) * dzdy
            
            a0 = a2 = xab / yab
            if j0 <= jc:
                a1 = xca / yca  # yca is not zero because j0-1 < jc
                if a0 > a1:
                    a0, a1 = a1, a0
            if jc < j1:
                a3 = xbc / ybc
                if a2 < a3:
                    a2, a3 = a3, a2

            for j in range(j0, j1):
                if j <= jc:
                    i0 = max(0, int(numpy.ceil(xa + (j - ya + int(a0 > 0)) * a0)))
                    i1 = min(C, int(numpy.floor(xa + (j - ya + int(a1 < 0)) * a1)))
                else:
                    i0 = 0
                    i1 = C
                
                if j >= jc:
                    i0 = max(i0, int(numpy.ceil(xb + (j - yb + int(a2 > 0)) * a2)))
                    i1 = min(i1, int(numpy.floor(xb + (j - yb + int(a3 < 0)) * a3)))

                z0 = zfirst + (i0 - xa + int(dzdx > 0)) * dzdx
                zfirst += dzdy
                
                if i1 - i0 < 10:
                    for i in range(i0, i1):
                        if self.z[j, i] < z0:
                            self.z[j, i] = z0
                            z0 += dzdx
                elif abs(dzdx) < 0.001 / (i1 - i0):
                    numpy.maximum(self.z[j, i0:i1], z0,
                                  out=self.z[j, i0:i1]
                                  )
                else:
                    numpy.maximum(self.z[j, i0:i1],
                                  numpy.arange(z0, z0 + dzdx * (i1 - i0 - 0.5), dzdx),
                                  out=self.z[j, i0:i1]
                                  )
          
        else:
              
            zfirst = za + (i0 - xa + int(dzdx > 0)) * dzdx
            
            a0 = a2 = yab / xab
            if i0 <= ic:
                a1 = yca / xca  # xca is not zero because i0-1 < ic
                if a0 > a1:
                    a0, a1 = a1, a0
            if ic < i1:
                a3 = ybc / xbc
                if a2 < a3:
                    a2, a3 = a3, a2

            for i in range(i0, i1):
                if i <= ic:
                    j0 = max(0, int(numpy.ceil(ya + (i - xa + int(a0 > 0)) * a0)))
                    j1 = min(R, int(numpy.floor(ya + (i - xa + int(a1 < 0)) * a1)))
                else:
                    j0 = 0
                    j1 = R
                
                if i >= ic:
                    j0 = max(j0, int(numpy.ceil(yb + (i - xb + int(a2 > 0)) * a2)))
                    j1 = min(j1, int(numpy.floor(yb + (i - xb + int(a3 < 0)) * a3)))

                z0 = zfirst + (j0 - ya + int(dzdy > 0)) * dzdy
                zfirst += dzdx

                if j1 - j0 < 10:
                    for j in range(j0, j1):
                        if self.z[j, i] < z0:
                            self.z[j, i] = z0
                            z0 += dzdy
                elif abs(dzdy) < 0.001 / (j1 - j0):
                    numpy.maximum(self.z[j0:j1, i], z0,
                                  out=self.z[j0:j1, i]
                                  )
                else:
                    numpy.maximum(self.z[j0:j1, i],
                                  numpy.arange(z0, z0 + dzdy * (j1 - j0 - 0.5), dzdy),
                                  out=self.z[j0:j1, i]
                                  )

    BLOCK_SIZE = (1, 3, 9, 27, 81, 243)

    def _fill_job_row( self, job, rt, irt, bs, j, i ):
        index = self.BLOCK_SIZE.index(bs)
        bs2 = bs // 2
        if i == 0 and (bs + j - 0.5) ** 2 + (bs2 - 0.5) ** 2 < rt ** 2:
            job.append((j + bs2 + 1, index, 0, 0))
            if j != -bs // 2 - 1:
                job.append((-(j + bs2 + 1), index, 0, 0))
            i = bs2
        while i + bs <= irt and (i + bs - 1) ** 2 + (j + bs - 1) ** 2 < (rt - 0.3) ** 2:
            job.append((j + bs2 + 1, index,   i + bs2 + 1, 0))
            job.append((j + bs2 + 1, index, -(i + bs2 + 1), 0))
            if j != -bs // 2 - 1:
                job.append((-(j + bs2 + 1), index,   i + bs2 + 1, 0))
                job.append((-(j + bs2 + 1), index, -(i + bs2 + 1), 0))
            i += bs
        bs = bs // 3
        if bs > 0:
            if j > 0:
                self._fill_job_row( job, rt, irt, bs, j, i )
            self._fill_job_row( job, rt, irt, bs, j + bs, i )
            self._fill_job_row( job, rt, irt, bs, j + 2 * bs, i )
        
                     
    def applyTool( self, radius, profile ):
        # profile is None for square end mill or sorted list of (radius, elevation)
        # Make a job as a sorted list of (y, source_index, x, elevation)
        rt = radius / self.sampleInterval
        border = self.border
        irt = min(border, int(numpy.ceil(rt)))
        job = []
        maxcol = max(border, (2000 - border) // 16 * 16)
        R, C = self.z.shape
        partial = [numpy.empty((R, min(maxcol + 2 * border, C)))]
        
        if C > maxcol:
            buff = numpy.empty((R, border))
            
        if profile is None:
            bs = 1
            while bs * 3 < math.sqrt(0.4) * rt and bs <= self.BLOCK_SIZE[-1]:
                partial.append(numpy.empty(partial[-1].shape))
                bs *= 3
            j = -(bs // 2) - 1
            while j < irt:
                while irt - j < bs:
                    bs //= 3
                self._fill_job_row( job, rt, irt, bs, j, 0 )
                j += bs
                
        else:
            pr = [p[0] for p in profile]
            pz = [p[1] for p in profile]
            for i in range(1, irt+1):
                z = - numpy.interp((i-1) * self.sample_interval, pr, pz)
                job.append(( i, 0,  0, z)) 
                job.append((-i, 0,  0, z)) 
                job.append(( 0, 0,  i, z)) 
                job.append(( 0, 0, -i, z)) 
                for j in range(1, i + 1):
                    r = math.sqrt((i - 1)**2 + (j - 1)**2)
                    if r < rt:
                        z = - numpy.interp(r * self.sample_interval, pr, pz)
                        job.append(( j, 0,  i, z)) 
                        job.append(( j, 0, -i, z)) 
                        job.append((-j, 0,  i, z)) 
                        job.append((-j, 0, -i, z)) 
                        if j < i:
                            job.append(( i, 0,  j, z)) 
                            job.append(( i, 0, -j, z)) 
                            job.append((-i, 0,  j, z)) 
                            job.append((-i, 0, -j, z)) 
        
        job.sort()
        
        offs = border
        while offs < C - border:
            cols = min( maxcol, C - border - offs )
            if offs > border:
                partial[0][:,:cols+2*border] = numpy.column_stack(
                                       (buff, self.z[:, offs:offs+cols+border])
                                     )
            else:
                partial[0][:,:cols+2*border] = self.z[:, offs-border:offs+cols+border]
            
            bs = 1
            for k in range(1, len(partial)):
                for m in range(bs, R-bs):
                    partial[k][m,:cols+2*border] = partial[k-1][m,:cols+2*border]
                    for i in range(-bs, bs+1, bs):
                         for j in range(-bs, bs+1, bs):
                             if i != 0 and j != 0:
                                 numpy.maximum(
                                     partial[k][m, bs:cols+2*border-bs],
                                     partial[k-1][m+j, bs+i:cols+2*border-bs+i],
                                     out = partial[k][m, bs:cols+2*border-bs]
                                     )
                bs *= 3

            if offs + cols < C - border:
                buff[:,:] = self.z[:, offs:offs+border]
                
            for m in range(border, R-border):
                for j, k, i, z in job:
                    numpy.maximum(
                        self.z[m, offs:offs+cols],
                        partial[k][m+j, border+i:cols+border+i] + z,
                        out = self.z[m, offs:offs+cols]
                        )
            offs += cols
                      
