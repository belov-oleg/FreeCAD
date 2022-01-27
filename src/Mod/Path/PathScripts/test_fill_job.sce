cd /home/oleg/soft/research/freecad/FreeCAD/src/Mod/Path/PathScripts/

job = csvRead("job");

scf(1); clf;

r = 15;
rx = r * sin(0:0.01:6.3);
ry = r * cos(0:0.01:6.3);
plot( rx + (rx > 0) - 0.5, ry + (ry > 0) - 0.5);

for jj = job',  // (source_index, x, y, elevation)
    bs = 3 ^ jj(2) - 0.05 * (4 - jj(2));
    xx = [-bs, -bs, bs,  bs, -bs] * 0.5 + jj(1)
    yy = [-bs,  bs, bs, -bs, -bs] * 0.5 + jj(3)
    plot(xx, yy);
end
