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

    
from PathScripts.PathLevelMap import LevelMap
import math
import numpy

class PathMesh():
    # It is a replacement for BoundBox, but also it contain already tessellated
    # set of faces.
    def __init__( self, faces, rotation=None ):
        MAX = 1e37
        self.XMin = MAX
        self.YMin = MAX
        self.ZMin = MAX
        self.XMax = -MAX
        self.YMax = -MAX
        self.ZMax = -MAX
        if rotation is None:
            self.matrix = None
        else:
            pass
            #TODO create matrix
        
        self.tessellated_faces = []
        for face in faces:
            vertices, facet_indices = face.tessellate( 0.01 )
            self.tessellated_faces.append((vertices, facet_indices))
            
            if not self.matrix is None:
                for v in vertices:
                    mr = self.matrix[0]
                    x = v[0] * mr[0] + v[1] * mr[1] + v[2] * mr[2]
                    mr = self.matrix[1]
                    y = v[0] * mr[0] + v[1] * mr[1] + v[2] * mr[2]
                    mr = self.matrix[2]
                    z = v[0] * mr[0] + v[1] * mr[1] + v[2] * mr[2]
                    self.XMin = min(self.XMin, x)
                    self.YMin = min(self.YMin, y)
                    self.ZMin = min(self.ZMin, z)
                    self.XMax = max(self.XMax, x)
                    self.YMax = max(self.YMax, y)
                    self.ZMax = max(self.ZMax, z)
                
            else:
                for v in vertices:
                    self.XMin = min(self.XMin, v[0])
                    self.YMin = min(self.YMin, v[1])
                    self.ZMin = min(self.ZMin, v[2])
                    self.XMax = max(self.XMax, v[0])
                    self.YMax = max(self.YMax, v[1])
                    self.ZMax = max(self.ZMax, v[2])

class LevelMapOp():
    def __init__( self, bound_box, sample_interval, final_depth, tool,
                 boundary_adjustment = 0):
        # tool is a tool object or tool radius
        # bound_box can be an instance of BoundBox or a tuple (faces, voids)
        # where faces and voids are lists of faces.
        
        self.area = None
        if isinstance(bound_box, tuple) and len(bound_box[0]) > 0:
            bound_box = PathMesh( bound_box[0] )
            self.area = bound_box
        
        max_size = max(bound_box.XMax - bound_box.XMin,
                       bound_box.YMax - bound_box.YMin) - 2 * boundary_adjustment
        if max_size <= 0:
            return
        sample_interval = max(sample_interval, 0.001 )
        while sample_interval < 0.0001 * max_size:
            sample_interval *= 2

        if type(tool) == int or type(tool) == float: 
            self.radius = tool
        else:
            self.radius = float(tool.Diameter) / 2

        border = max(2, int(numpy.ceil(self.radius / sample_interval)))

        extra = 0.5 * sample_interval - boundary_adjustment
        self.levelMap = LevelMap(bound_box.XMin - extra, 
                                 bound_box.XMax + extra,
                                 bound_box.YMin - extra, 
                                 bound_box.YMax + extra,
                                 final_depth, sample_interval, border
                                 )
  
    def raiseModel(self, model ):   
        # raise points of the map in accordance with the model
        if model.TypeId.startswith("Mesh"):
            if self.levelMap.includes(model.Mesh.BoundBox):
                for fa, fb, fc in model.Mesh.Facets.Points:
                    self.levelMap.add_facet(fa, fb, fc)  #TODO Not tested!!!
        else:
            if hasattr(model, 'Shape'):
                shape = model.Shape
            else:
                shape = model
            if self.levelMap.includes(shape.BoundBox):
                vertices, facet_indices = shape.tessellate(
                                    0.25 * self.levelMap.sampleInterval)

                for f in facet_indices:
                    self.levelMap.add_facet(vertices[f[0]], 
                                            vertices[f[1]], 
                                            vertices[f[2]])

    def raisePathMesh(self, path_mesh ):
        for vertices, facet_indices in path_mesh.tessellated_faces:
            for f in facet_indices:
                self.levelMap.add_facet(vertices[f[0]], 
                                        vertices[f[1]], 
                                        vertices[f[2]])
                

    def applyTool(self, tool ):
        radius = float(tool.Diameter) / 2
        if type(tool) == int or type(tool) == float:
            self.levelMap.applyTool( radius, None )
        
        elif hasattr(tool, 'ShapeName'):
            sample_interval = self.levelMap.sampleInterval
            tool_level_map = LevelMap(0, radius,
                                      0, sample_interval, -float(tool.Length),
                                      sample_interval, 0)
            tool_level_map.matrix = [[1/sample_interval, 0, 0],
                                     [0, 1/sample_interval, 0],
                                     [0, 0,               -1]]
            vertices, facet_indices = tool.Shape.tessellate(
                0.25 * sample_interval)
            for f in facet_indices:
                tool_level_map.add_facet(vertices[f[0]], 
                                         vertices[f[1]], 
                                         vertices[f[2]])
            profile = []
            for i in range (0, tool_level_map.columns()):
                profile.append((i * sample_interval, -tool_level_map.z[0, i]))
                
            self.levelMap.applyTool( radius, profile )
            
        else:
            # For end mill:
            self.levelMap.applyTool( radius, None ) 
            

    def getContourMap(self, z, dep_offset = 0 ):
        m = self.levelMap.getContourMap( z )
        m.z += dep_offset
        return m
      
