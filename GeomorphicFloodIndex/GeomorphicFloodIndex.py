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

# Import the PyQt and QGIS libraries
from qgis.PyQt.QtCore import * 
from qgis.PyQt.QtGui import *
from qgis.PyQt.QtWidgets import QAction
from qgis.core import *
# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
#from FlowPathDown_BBDialog import FlowPathDown_BBDialog
from .Ui_GeomorphicFloodIndex import Ui_GeomorphicFloodIndex
from .doGeomorphicFloodIndex import *

class GeomorphicFloodIndex: 

  def __init__(self, iface):
    # Save reference to the QGIS interface
    self.iface = iface

  def initGui(self):  
    # Create action that will start plugin configuration
    self.action = QAction(QIcon(":/plugins/GeomorphicFloodIndex/icona.png"), \
        "GeomorphicFloodArea", self.iface.mainWindow())

    # connect the action to the run method
    # QObject.connect(self.action, SIGNAL("activated()"), self.run) 
    self.action.triggered.connect(self.run)

    # Add toolbar button and menu item
    self.iface.addToolBarIcon(self.action)

    self.iface.addPluginToMenu("GFA", self.action)

  def unload(self):
    # Remove the plugin menu item and icon
    self.iface.removePluginMenu("GFA",self.action)
    self.iface.removeToolBarIcon(self.action)
  

  # run method that performs all the real work
  def run(self): 
   
    dlg = GeomorphicFloodIndexDialog(self.iface)
    dlg.exec_()
        
