cd /home/oleg/soft/research/freecad/FreeCAD/src/Mod/Path/PathScripts/

s = csvRead("s");

triangles = [s, s(:,1:3)];
x = triangles(:,1:3:$) * 0.95 + mean(triangles(:,1:3:9), 'c') * [0.05, 0.05, 0.05, 0.05]
y = triangles(:,2:3:$) * 0.95 + mean(triangles(:,2:3:9), 'c') * [0.05, 0.05, 0.05, 0.05]
z = triangles(:,3:3:$) * 0.95 + mean(triangles(:,3:3:9), 'c') * [0.05, 0.05, 0.05, 0.05]

scf(1);clf;param3d1(x',y',z');


lev = csvRead("lev");
[r,c]=size(lev);
mesh(-15+((0:r-1) + 0.5)*30/(r-1), -15+((0:c-1) + 0.5)*30/(c-1), lev');
h=gce(); //get handle on current entity (here the surface)
a=gca(); //get current axes
a.cube_scaling=0;
a.isoview=1;
h.color_flag=1; //color according to z
h.color_mode=-2;  //remove the facets boundary by setting color_mode to white color
h.color_flag=2; //color according to given colors
h.color_mode = -1; // put the facets boundary back by setting color_mode to black color
f=gcf();//get the handle of the parent figure
f.color_map=hsvcolormap(512);
