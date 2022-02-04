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

    
from PathScripts.PathLevelMap import PathLevelMap
import math

class LevelMapOp():
    def __init__( self, bound_box, sample_interval, final_depth, tool ):
        # tool is a tool object or tool radius
        max_size = max(bound_box.XMax - bound_box.XMin,
                       bound_box.YMax - bound_box.YMin)
        if max_size == 0:
            return
        sample_interval = max(sample_interval, 0.001 )
        while sample_interval < 0.0001 * max_size:
            sample_interval *= 2

        if type(tool) == int or type(tool) == float: 
            radius = tool
        else:
            radius = float(tool.Diameter) / 2

        border = max(2, int(numpy.ceil(radius / sample_interval))

        self.levelMap = LevelMap(bound_box.XMin - 0.5 * sample_interval, 
                                 bound_box.XMax + 0.5 * sample_interval,
                                 bound_box.YMin - 0.5 * sample_interval, 
                                 bound_box.YMin + 0.5 * sample_interval,
                                 obj.FinalDepth.Value, sample_interval, border
                                 )
  
    def raiseModel( model ):   
        # raise points of the map in accordance with the model
        if model.TypeId.startswith("Mesh"):
            for fa, fb, fc in model.Mesh.Facets.Points:
                self.levelMap.add_facet(fa, fb, fc)  #TODO Not tested!!!
        else:
            if hasattr(model, 'Shape'):
                shape = model.Shape
            else:
                shape = model
            vertices, facet_indices = shape.tessellate(0.25 * sample_interval)

            for f in facet_indices:
                self.levelMap.add_facet(vertices[f[0]], 
                                        vertices[f[1]], 
                                        vertices[f[2]])
                

    def applyTool( tool ):
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
            vertices, facet_indices = self.tool.Shape.tessellate(
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
            

    def getContourMap( z, dep_offset = 0 ):
        m = self.levelMap.getContourMap( layDep )
        m.z += dep_offset
        return m
      
