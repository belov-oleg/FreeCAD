# This is a standalone test, use it from command line.
# This test will be rewrited as an automatic one and will be moved into the appropriate directory.
# For visualization of test results the appropriate scilab scripts should be used.
import PathLevelMap
import math
import numpy

def sphere( x0, y0, z0, r, n ):
    # return list of facets
    ans = []
    prev = [(0, 0, -r)] * n
    for par in range(1, n + 1):
        lat = math.pi * (par / n - 0.5)
        z = r * math.sin( lat )
        rp = r * math.cos( lat )
        
        pp = []
        for mer in range(0, n):
            lon = 2 * math.pi * (mer + par * 0.5) / n
            pp.append(( rp * math.cos(lon), rp * math.sin(lon), z))
        
        for i in range(0, n):
            if par != n:
                ans.append((pp[i-1], pp[i], prev[i]))
            if par != 1:
                ans.append((prev[i-1], pp[i-1], prev[i]))
            
        prev = pp
    return ans

def test_add_facet():
    lm = PathLevelMap.LevelMap( -15, 15, -15, 15, -9, 0.05, 10 )

    fp = open('s', 'w')
    for t in sphere( 0, 0, 0, 10, 10 ):
        fp.write( ("%.2f, " * 8 + "%.2f\n") % (t[0] + t[1] + t[2]))
        lm.add_facet( t[0], t[1], t[2] )   
    fp.close

    fp = open('lev', 'w')
    for i in range(0, lm.z.shape[0]):
        fp.write( ', '.join(['%.3f' % z for z in lm.z[:,i]]) + "\n")
    fp.close



def test_fill_job( border, rt ):
    BLOCK_SIZE = (1, 3, 9, 27, 81, 243)
    lm = PathLevelMap.LevelMap( -15, 15, -15, 15, -9, 0.1, border )
    irt = min(border, int(numpy.ceil(rt)))
    job = []
    bs = 3
    while bs * 3 < math.sqrt(0.4) * rt and bs <= BLOCK_SIZE[-1]:
        bs *= 3
    j = -(bs // 2) - 1
    while j < irt:
        while irt - j < bs:
            bs //= 3
        lm._fill_job_row( job, rt, irt, bs, j, 0 )
        j += bs
    fp = open('job', 'w')
    for ji in job:
        fp.write( "%i, %i, %i, %.3f\n" % ji )
    fp.close

test_fill_job( 20, 15)
