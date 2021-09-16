# -*- coding: utf-8 -*-
"""
/***************************************************************************
GeomorphicFloodArea
A QGIS plugin
GFA
 -------------------
version                : 2.0
author                 : Raffaele Albano
contact                : http://www2.unibas.it/raffaelealbano/?page_id=115
***************************************************************************/
    
/***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************/
This script initializes the plugin, making it known to QGIS.
"""

def classFactory(iface): 
  
  from .GeomorphicFloodIndex import GeomorphicFloodIndex 
  return GeomorphicFloodIndex(iface)
