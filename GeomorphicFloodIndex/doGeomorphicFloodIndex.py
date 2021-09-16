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
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *
from qgis.core import *
import qgis.utils
from .Ui_GeomorphicFloodIndex import *
from qgis.PyQt.QtWidgets import QProgressBar, QDialog, QMessageBox, QFileDialog
import os, sys, time,  math
from osgeo import gdal, ogr
from osgeo.gdalconst import *
from scipy import ndimage
from osgeo import osr
import numpy
#from sklearn import metrics
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import webbrowser

import shlex
import subprocess
import platform
import tempfile
import shutil



class GeomorphicFloodIndexDialog(QDialog, Ui_GeomorphicFloodIndex):

    def __init__(self, iface):
        QDialog.__init__(self)
        self.iface = iface
        self.setupUi(self)
        self.NOVALUE=-340282346638528859811704183484516925440.000000
        self.checkBoxManualSet.setChecked(True)
        self.checkBoxCalibration.setChecked(False)
        # self.connect(self.toolButtonDemCon, SIGNAL("clicked()"), self.demconFile)
        self.toolButtonDemCon.clicked.connect(self.demconFile)
        # self.connect(self.toolButtonDemVoid, SIGNAL("clicked()"), self.demvoidFile)
        self.toolButtonDemVoid.clicked.connect(self.demvoidFile)
        # self.connect(self.toolButtonFlowDir, SIGNAL("clicked()"), self.flowdirFile)
        self.toolButtonFlowDir.clicked.connect(self.flowdirFile)
        # self.connect(self.toolButtonFlowAcc, SIGNAL("clicked()"), self.flowaccFile)
        self.toolButtonFlowAcc.clicked.connect(self.flowaccFile)
        # self.connect(self.toolButtonSFM, SIGNAL("clicked()"), self.sfmFile)
        self.toolButtonSFM.clicked.connect(self.sfmFile)

        # self.connect(self.checkBoxCalibration,SIGNAL("clicked()"),self.calibration_clicked)
        # self.checkBoxCalibration.stateChanged.connect(self.calibration_clicked)
        # self.connect(self.checkBoxManualSet,SIGNAL("clicked()"),self.manualset_clicked)
        # self.checkBoxManualSet.stateChanged.connect(self.manualset_clicked)

        # self.connect(self.btnOutput, SIGNAL("clicked()"), self.outFile)
        self.btnOutput.clicked.connect(self.outFile)
        # self.connect(self.btnOutputBin, SIGNAL("clicked()"), self.outFileBin)
        self.btnOutputBin.clicked.connect(self.outFileBin)
        # self.connect(self.btnOutputWD, SIGNAL("clicked()"), self.outFileWD)
        self.btnOutputWD.clicked.connect(self.outFileWD)

        # self.connect(self.buttonBox, SIGNAL("accepted()"),self.accept)
        self.buttonBox.accepted.connect(self.accept)
        
        # QObject.connect(self.buttonBox, SIGNAL("rejected()"),self, SLOT("reject()"))
        self.buttonBox.rejected.connect(self.reject)
        # QObject.connect(self.buttonBox, SIGNAL("helpRequested()"),self.call_help)
        self.buttonBox.helpRequested.connect(self.call_help)

        mapCanvas = self.iface.mapCanvas()
        # init dictionaries of items:
        self.rastItems = {}
        for i in range(mapCanvas.layerCount()):
            layer = mapCanvas.layer(i)
            if ( layer.type() == layer.RasterLayer ):
		# read  layers
                provider = layer.dataProvider()
                self.cmbDemCon.addItem(layer.source())
                self.cmbDemVoid.addItem(layer.source())
                self.cmbFlowDir.addItem(layer.source())
                self.cmbFlowAcc.addItem(layer.source())
                self.cmbSFM.addItem(layer.source())

        self.textEdit.clear()
        self.cmbFlowDirCoding.addItem("ESRI")
        self.cmbFlowDirCoding.addItem("HyGrid2k2")
        self.cmbFlowDirCoding.addItem("TauDEM")

        self.cmbSlopeType.addItem("slope_gdaldem")
        self.cmbSlopeType.addItem("slope")
        self.cmbSlopeType.setVisible ( False )
        self.lblSlope_3.setVisible ( False )
        self.cmbChannelType.addItem("channel_ASk")
        self.cmbChannelType.addItem("channel_FAt")

        self.checkBoxWD.setVisible ( False )
        self.btnOutputWD.setVisible ( False )
        self.lineOutputWD.setVisible ( False )
        self.doubleSpinBoxN.setValue(0.354429752)

    def my_roc_curve(self,floodhazard_1d, ln_hronH_norm_1d,passo):

        thresholds=numpy.linspace(1,-1 ,2/passo)
        nth=len(thresholds)
        fpr=numpy.zeros((1,nth)).flatten()
        tpr=numpy.zeros((1,nth)).flatten()
        self.progressBar.setRange(1,nth)
        for i in range(nth):
            self.progressBar.setValue(i)
            t=thresholds[i]
            tp_t=numpy.logical_and(ln_hronH_norm_1d>t,floodhazard_1d>0).sum()
            tn_t=numpy.logical_and(ln_hronH_norm_1d<=t,floodhazard_1d==0).sum()
            fp_t=numpy.logical_and(ln_hronH_norm_1d>t,floodhazard_1d==0).sum()
            fn_t=numpy.logical_and(ln_hronH_norm_1d<=t,floodhazard_1d>0).sum()

            fpr[i]=fp_t/float(fp_t+tn_t)
            tpr[i]=tp_t/float(tp_t+fn_t)

        return fpr,tpr,thresholds

    def my_auc(self,fpr,tpr):
        auc=0
        n=len(fpr)
        for i in range(n-1):
            auc+=(fpr[i+1]-fpr[i])*( tpr[i+1]+tpr[i] )

        auc*=0.5
        return auc

    def call_help(self):
       # qgis.utils.showPluginHelp()
       help_path = os.path.join(os.path.dirname(__file__), 'help', 'help.html')
       webbrowser.open(help_path)
    def outFile(self):
        "Display file dialog for output file"
        self.lineOutput.clear()
        outName, filter_string = QFileDialog.getSaveFileName(self, "GFI output file",".", "GeoTiff (*.tif)")

        if len(outName)>0:
                self.lineOutput.clear()
                self.lineOutput.insert(outName)

        return outName
    def outFileBin(self):
        "Display file dialog for output file"
        self.lineOutputBin.clear()
        outNameBin, filter_string = QFileDialog.getSaveFileName(self, "GFI output file",".", "GeoTiff (*.tif)")

        if len(outNameBin)>0:
                self.lineOutputBin.clear()
                self.lineOutputBin.insert(outNameBin)

        return outNameBin

    def outFileWD(self):
        "Display file dialog for output file"
        self.lineOutputWD.clear()
        outNameWD, filter_string = QFileDialog.getSaveFileName(self, "WD output file",".", "GeoTiff (*.tif)")

        if len(outNameWD)>0:
                self.lineOutputWD.clear()
                self.lineOutputWD.insert(outNameWD)

        return outNameWD

    def demconFile(self):
        "Display file dialog for input file"
        demconName, filter_string = QFileDialog.getOpenFileName(self, "FILL Dem input file",".", "ESRI ascii (*.txt);; GeoTiff (*.tif);;All files (*.*)")
        if len(demconName)>0:
                self.cmbDemCon.insertItem (0, demconName)
                self.cmbDemCon.setCurrentIndex(0)
        return demconName

    def demvoidFile(self):
        "Display file dialog for input file"
        demvoidName, filter_string = QFileDialog.getOpenFileName(self, "Dem input file",".", "ESRI ascii (*.txt);; GeoTiff (*.tif);;All files (*.*)")
        if len(demvoidName)>0:
                self.cmbDemVoid.insertItem (0, demvoidName)
                self.cmbDemVoid.setCurrentIndex(0)
        return demvoidName

    def flowdirFile(self):
        "Display file dialog for input file"
        flowdirName, filter_string = QFileDialog.getOpenFileName(self, "Flowdir input file",".", "ESRI ascii (*.txt);; GeoTiff (*.tif);;All files (*.*)")
        if len(flowdirName)>0:
                self.cmbFlowDir.insertItem (0, flowdirName)
                self.cmbFlowDir.setCurrentIndex(0)
        return flowdirName

    def flowaccFile(self):
        "Display file dialog for input file"
        flowaccName, filter_string = QFileDialog.getOpenFileName(self, "Flowacc input file",".", "ESRI ascii (*.txt);; GeoTiff (*.tif);;All files (*.*)")
        if len(flowaccName)>0:
                self.cmbFlowAcc.insertItem (0, flowaccName)
                self.cmbFlowAcc.setCurrentIndex(0)
        return flowaccName

    def sfmFile(self):
        "Display file dialog for input file"
        sfmName, filter_string = QFileDialog.getOpenFileName(self, "Standard Flood Map calibration file",".", "ESRI ascii (*.txt);; GeoTiff (*.tif);;All files (*.*)")
        if len(sfmName)>0:
                self.cmbSFM.insertItem (0, sfmName)
                self.cmbSFM.setCurrentIndex(0)
        return sfmName

    # def calibration_clicked(self, state):
        # self.checkBoxManualSet.nextCheckState()
    # def manualset_clicked(self):
        # self.checkBoxCalibration.nextCheckState()
###########################################################################################
    def from_flowacc_to_stream(self,tab_flowacc,threshold_cell):

        nrows=numpy.shape(tab_flowacc)[0]
        ncols=numpy.shape(tab_flowacc)[1]

        ar_flowacc=numpy.reshape(tab_flowacc,ncols*nrows)
        ar_stream=numpy.array(ar_flowacc)

        ar_stream[ar_flowacc>=threshold_cell]=1
        ar_stream[numpy.where((ar_flowacc<threshold_cell) & (ar_flowacc>-1))]=0

        total_cell=len(ar_flowacc[ar_flowacc>-1.])
        stream_cell=len(ar_stream[ar_stream==1])

       # self.textEdit.append( 'total_cell='+total_cell+'stream_cell='+stream_cell,'Drainage density=',float(stream_cell)/float(total_cell))

        tab_stream=numpy.reshape(ar_stream,(nrows,ncols))
        return tab_stream
        #f = file(file_stream_grid, 'w')
        #io.write_array(f, tab_stream)
        #f.close()



    def bingrid_to_label(self,tab, file_label='label.dat', write_file=False):
        """
        * Objective
          Replace values of the grid by the label of the cells.
          The label are assigned from 0 to nb_cell, from West to East, North to South.
        * Input
          file_bin_grid: A binary grid file (whatever it is).
        * Output
          tab: a 1D array with nb_cell components.
        """


        nrows=numpy.shape(tab)[0]
        ncols=numpy.shape(tab)[1]
        tab=numpy.reshape(tab,ncols*nrows)
        ind=numpy.where(tab>-99)
        tab_label=numpy.arange(len(tab[ind]))
        for i in tab_label:
            tab[ind[0][i]]=i
        tab=numpy.reshape(tab,(nrows,ncols))
        tab=tab.astype('int32')
        numpy.savetxt('tablabel.txt',tab)
        if write_file:
            f = file(file_label, 'w')
            io.write_array(f, tab)
            f.close()

        return tab


    def create_channel_slope_file(self,file_flowdir, file_DEM, file_slope_degree,di,dj):
        """Compute the channel slopes from 8D flow direction and a DEM.

        Calculate the slope from the centre of each cell in the catchment DEM
        to it's downstream neighbour. The calculated slopes are written to
        an ArcGIS Float grid file. This file is required by the TOPKAPI model
        for the generation of a parameter file.

        Parameters
        ----------
        file_flowdir : string
            Name of an ArcGIS binary Float file containing the flow
            direction raster.
        file_DEM : string
            Name of an ArcGIS binary Float file containing the DEM raster.
        file_slope_degree : string
            Name of an ArcGIS binary Float file for the output raster.

        """

        gdal.AllRegister() #driver registration
        self.textEdit.append( 'Opening ' + file_DEM)
        ds_dem = gdal.Open(file_DEM, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
        if ds_dem is None:
            self.textEdit.append( 'Could not open ' + file_DEM)
            sys.exit(1)

        self.textEdit.append( 'Opening ' + file_flowdir)
        ds_flowdir = gdal.Open(file_flowdir, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
        if ds_flowdir is None:
            self.textEdit.append( 'Could not open ' + file_flowdir)
            sys.exit(1)


        cols = ds_dem.RasterXSize# number of columns
        rows = ds_dem.RasterYSize# number of rows
        bands = ds_dem.RasterCount# number of bands

        geotransform = ds_dem.GetGeoTransform()#georeference functions
        originX = geotransform[0]#top left x
        originY = geotransform[3]#top left y
        cellsize = geotransform[5]#pixel resolution

        band_dem=ds_dem.GetRasterBand(1)#extract band 1
        dem=band_dem.ReadAsArray(0,0,cols,rows)
        band_flowdir=ds_flowdir.GetRasterBand(1)#extract band 1
        flowdir=band_flowdir.ReadAsArray(0,0,cols,rows)
        tab_label = self.bingrid_to_label(dem)

        tab_slope_degree = numpy.array(tab_label, numpy.float32)
        nrows=rows
        ncols=cols
        Xcell=cellsize
        self.progressBar.setRange(1,nrows)
        self.progressBar.show()
        for i in range(nrows-1):

            self.progressBar.setValue(i)


            for j in range(ncols-1):

                label = tab_label[i,j]
                direction = flowdir[i,j]

                if label > 0:

                    x=i+di[direction]
                    y=j+dj[direction]
                    if di[direction]==0 or dj[direction]==0:
                        dist=Xcell
                    else:
                        dist=Xcell*(2**0.5)


                    if tab_label[x,y] >= 0:
                        if dem[x,y] < dem[i,j]:
                            tab_slope_degree[i,j] = numpy.arctan((dem[i,j]-dem[x,y])
                                                                /dist) * 180./numpy.pi
                        if dem[x,y] == dem[i,j]:
                            tab_slope_degree[i,j] = 0.0
                        if dem[x,y] > dem[i,j]:

                            tab_slope_degree[i,j] = numpy.arctan((dem[i,j]-dem[x,y])
                                                                /dist) * 180./numpy.pi
                    else:
                        tab_slope_degree[i,j] = -9999

        self.progressBar.setValue(nrows)
        driver1=gdal.GetDriverByName('GTiff')
        driver1.Register()

        target_ds=driver1.Create(file_slope_degree,cols,rows,1,gdal.GDT_Float32)
        target_ds.SetGeoTransform(geotransform)
        band_target=target_ds.GetRasterBand(1)
        band_target.WriteArray(tab_slope_degree)
        band_target.FlushCache()
        band_target=None
        target_ds=None

        band_dem=None
        ds_dem=None

        band_flowdir=None
        ds_flowdir=None


    def loadOutputFile(self,outFile):
        "Load map in TOC"
        fileInfo = QFileInfo(outFile)
        baseName = fileInfo.baseName()
        rlayer = QgsRasterLayer(outFile, baseName)
        if not rlayer.isValid():
            self.textEdit.append("Layer failed to load!")
        if QGis.QGIS_VERSION_INT < 10900:
            rlayer.setDrawingStyle(QgsRasterLayer.SingleBandPseudoColor)
            rlayer.setColorShadingAlgorithm(QgsRasterLayer.FreakOutShader)

        QgsMapLayerRegistry.instance().addMapLayer(rlayer)


    def writeOutputGeoTiff(self,arrayData, transform, prj,rows, cols, outFile):
        "Write the given array data to the file 'outfile' with the given extent."
        try:
            format = "GTiff"
            driver = gdal.GetDriverByName( format )
            metadata = driver.GetMetadata()
            if gdal.DCAP_CREATE in metadata and metadata[gdal.DCAP_CREATE] == 'YES':
                pass
            else:
                QMessageBox.information(None,"info","Driver %s does not support Create() method." % format)
                return False
            outDataset = driver.Create(str(outFile), cols, rows, 1, gdal.GDT_Float32)
            outRaster=outDataset.GetRasterBand(1)
            outRaster.WriteArray(arrayData)
            outRaster.SetNoDataValue(self.NOVALUE)
            outDataset.SetGeoTransform(transform)
            outDataset.SetProjection(prj)
            return True
        except:
            QMessageBox.information(None,"Exiting gracefully","I can't write %s texture file." % outFile)
            return
###########################################################################################
    def ProcessRaster(self,paramDEMCON,paramDEMVOID,paramFLOWDIR, paramFLOWACC):
        "Register all of the GDAL drivers"
        gdal.AllRegister()
        # open the image
        paramDEMCON=str(paramDEMCON)
        paramDEMVOID=str(paramDEMVOID)
        paramFLOWDIR=str(paramFLOWDIR)
        paramFLOWACC=str(paramFLOWACC)

        dsDEM = gdal.Open(paramDEMCON,GA_ReadOnly)
        if dsDEM is None:
            QMessageBox.information(None,"Exiting gracefully","Could not open raster %s!" % paramDEM)


        # get image size
        rows  = dsDEM.RasterYSize
        cols  = dsDEM.RasterXSize
        bands = dsDEM.RasterCount
        self.textEdit.append('rows: '+str(rows)+' columns: '+str(cols))
        if bands >1:
            QMessageBox.information(None,"Exiting gracefully","Raster %s has %d bands!" % (paramDEM,bands))

        # get georeference info
        transform = dsDEM.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = transform[5]
####################################################

###################################################
    def accept(self):
        # Called when "OK" button pressed
        file_demcon=self.cmbDemCon.currentText()
        file_demvoid=self.cmbDemCon.currentText()
        file_flowdir=self.cmbFlowDir.currentText()
        file_flowacc=self.cmbFlowAcc.currentText()

        if os.path.isfile(file_demcon) and os.path.isfile(file_demvoid) and os.path.isfile(file_flowdir) and os.path.isfile(file_flowacc):



            t_ln_hronH= float(self.doubleSpinBoxTheshold.value())
            th_channel=float(self.doubleSpinBoxThesholdChannel.value())
            file_output=self.lineOutput.text()
            if len(file_output)==0:
                file_output="out.tif"
            flowDirCoding=self.cmbFlowDirCoding.currentText()
            channel_type=self.cmbChannelType.currentText()
            slope_type=self.cmbSlopeType.currentText()
            self.textEdit.append("Starting...")
            debug=0
            calibration=0

            if self.checkBoxDebug. isChecked():
                debug=1
            if self.checkBoxCalibration. isChecked():
                calibration=1
            #######################################################

            if (flowDirCoding=="ESRI"):
                        #row offset
                        di={1:0, 2:1, 4:1, 8:1, 16:0, 32:-1, 64:-1, 128:-1} #rows indexes start from upper to bottom!!!
                        #col offset
                        dj={1:1, 2:1, 4:0, 8:-1, 16:-1, 32:-1, 64:0, 128:1}
            elif (flowDirCoding=="HyGrid2k2"):
                #row offset
                di={0:0, 315:1, 270:1, 225:1, 180:0, 135:-1, 90:-1, 45:-1}
                #col offset
                dj={0:1, 315:1, 270:0, 225:-1, 180:-1, 135:-1, 90:0, 45:1}
            elif(flowDirCoding=="TauDEM"):
                #row offset
                di={1:0, 8:1, 7:1, 6:1, 5:0, 4:-1, 3:-1, 2:-1}
                #col offset
                dj={1:1, 8:1, 7:0, 6:-1, 5:-1, 4:-1, 3:0, 2:1}


            gdal.AllRegister()
            self.textEdit.append( 'Opening ' + file_demcon)
            ds_demcon = gdal.Open(file_demcon, GA_ReadOnly)
            if ds_demcon is None:
                self.textEdit.append( 'Could not open ' + file_demcon)
               # sys.exit(1)
                return 1

            self.textEdit.append( 'Opening ' + file_demvoid)
            ds_demvoid = gdal.Open(file_demvoid, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
            if ds_demvoid is None:
                self.textEdit.append( 'Could not open ' + file_demvoid)
                #sys.exit(1)
                return 1

            self.textEdit.append( 'Opening ' + file_flowdir)
            ds_flowdir = gdal.Open(file_flowdir, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
            if ds_flowdir is None:
                self.textEdit.append( 'Could not open ' + file_flowdir)
                #sys.exit(1)
                return 1

            self.textEdit.append( 'Opening ' + file_flowacc)
            ds_flowacc = gdal.Open(file_flowacc, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
            if ds_flowacc is None:
                self.textEdit.append( 'Could not open ' + file_flowacc)
                #sys.exit(1)
                return 1

            if calibration==1:
                file_floodhazard=self.cmbSFM.currentText()
                self.textEdit.append( 'Opening ' + file_floodhazard)
                ds_floodhazard = gdal.Open(file_floodhazard, GA_ReadOnly)# apri il file ed ottieni un oggetto dataset
                if ds_floodhazard is None:
                    self.textEdit.append( 'Could not open ' + file_floodhazard)
                    #sys.exit(1)
                    return 1
        ######  Check inpur file dimensions ######
            cols_demcon = ds_demcon.RasterXSize# number of columns
            rows_demcon = ds_demcon.RasterYSize# number of rows
            geotransform_demcon = ds_demcon.GetGeoTransform()#geotrasform functions list
            cellsize_demcon = geotransform_demcon[5]#pixel resolution

            cols_demvoid = ds_demvoid.RasterXSize# number of columns
            rows_demvoid = ds_demvoid.RasterYSize# number of rows
            geotransform_demvoid = ds_demvoid.GetGeoTransform()#geotrasform functions list
            cellsize_demvoid = geotransform_demvoid[5]#pixel resolution

            cols_flowdir = ds_flowdir.RasterXSize# number of columns
            rows_flowdir = ds_flowdir.RasterYSize# number of rows
            geotransform_flowdir = ds_flowdir.GetGeoTransform()#geotrasform functions list
            cellsize_flowdir = geotransform_flowdir[5]#pixel resolution

            cols_flowacc = ds_flowacc.RasterXSize# number of columns
            rows_flowacc = ds_flowacc.RasterYSize# number of rows
            geotransform_flowacc = ds_flowacc.GetGeoTransform()#geotrasform functions list
            cellsize_flowacc = geotransform_flowacc[5]#pixel resolution

            if cols_demcon!=cols_demvoid or rows_demcon!=rows_demvoid or cols_demcon!=cols_flowdir or rows_demcon!=rows_flowdir or cols_demcon!=cols_flowacc or rows_demcon!=rows_flowacc  :
                self.textEdit.append( 'Error : Input files have different number of rows or columns ' )
                return 1

            if cellsize_demcon!=cellsize_demvoid or cellsize_demcon!=cellsize_flowdir or cellsize_demcon!=cellsize_flowacc :
                self.textEdit.append( 'Error : Input files have different cell size ' )
                return 1

            if calibration==1:
                cols_calib = ds_floodhazard.RasterXSize# number of columns
                rows_calib = ds_floodhazard.RasterYSize# number of rows
                geotransform_calib = ds_floodhazard.GetGeoTransform()#geotrasform functions list
                cellsize_calib = geotransform_calib[5]#pixel resolution
                if cols_demcon!=cols_calib or rows_demcon!=rows_calib  :
                     self.textEdit.append( 'Error : Standard flood hazard map  has wrong number of rows or columns ' )
                     return 1

                if cellsize_demcon!=cellsize_calib :
                     self.textEdit.append( 'Error : Standard flood hazard map has wrong cell size ' )
                     return 1

            self.textEdit.append( 'Calculating slope... ' )

            temp_path = tempfile.mkdtemp()
            temp_slope = os.path.join(temp_path, "slope.tif")

            file_demcon1 = file_demcon

            if os.name == "nt":
                temp_slope = temp_slope.replace("\\", "\\\\")
                file_demcon1 = file_demcon.replace("\\", "\\\\")


            if os.path.isfile(temp_slope):
                os.remove(temp_slope)

            if slope_type=="slope_gdaldem":
                gdaldem_command="gdaldem slope "+ file_demcon1 + " " + temp_slope + " " # calcola lo slope con gdaldem

                cmd = shlex.split(gdaldem_command)

                if sys.platform == 'win32':
                    si = subprocess.STARTUPINFO()
                    si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                    subprocess.call(cmd, startupinfo=si)
                elif sys.platform == 'darwin':
                    gdaldem_command1 = "/Library/Frameworks/GDAL.framework/Programs/gdaldem slope "+ file_demcon + " " + temp_slope +  " "
                    cmd1 = shlex.split(gdaldem_command1)
                    subprocess.call(cmd1)
                else:
                    subprocess.call(cmd)


            else:
                self.create_channel_slope_file(file_flowdir, file_demcon, temp_slope,di,dj)

            self.textEdit.append( 'slope done')

            self.textEdit.append( 'Opening slope' )
            ds_G = gdal.Open(temp_slope, GA_ReadOnly)
            if ds_G is None:
                self.textEdit.append( 'Could not open slope' )
                sys.exit(1)

            cols = ds_demcon.RasterXSize# number of columns
            rows = ds_demcon.RasterYSize# number of rows
            bands = ds_demcon.RasterCount# number of bands

            geotransform = ds_demcon.GetGeoTransform()#geotrasform functions list
            prj_demcon = ds_demcon.GetProjectionRef()
            prj_demvoid = ds_demvoid.GetProjectionRef()
            prj_flowdir = ds_flowdir.GetProjectionRef()
            prj_flowacc = ds_flowacc.GetProjectionRef()
            
            
            if calibration == 1 :
                prj_floodhazard = ds_floodhazard.GetProjectionRef()
                if ( prj_demcon != prj_demvoid or prj_demcon != prj_flowdir or prj_demcon != prj_flowacc or prj_demcon != prj_floodhazard or  
                    prj_demvoid != prj_flowdir or prj_demvoid != prj_flowacc or prj_demvoid != prj_floodhazard or prj_flowdir != prj_flowacc or 
                    prj_flowdir != prj_floodhazard or prj_flowacc != prj_floodhazard ) :
                    self.textEdit.append('One or more input have different projection')
                    self.textEdit.append('Dem Projection :' + prj_demvoid)
                    self.textEdit.append('Filled DEM Projection :' + prj_demcon)
                    self.textEdit.append('Flow Direction Projection :' + prj_demcon)
                    self.textEdit.append('Flow Accumulation Projection :' + prj_flowacc)
                    self.textEdit.append('Flood Hazard Projection :' + prj_floodhazard)  
                    sys.exit(1)
            else :
                if ( prj_demcon != prj_demvoid or prj_demcon != prj_flowdir or prj_demcon != prj_flowacc or 
                    prj_demvoid != prj_flowdir or prj_demvoid != prj_flowacc or prj_flowdir != prj_flowacc ) :
                    self.textEdit.append('One or more input have different projection')
                    self.textEdit.append('Dem Projection :' + prj_demvoid)
                    self.textEdit.append('Filled DEM Projection :' + prj_demcon)
                    self.textEdit.append('Flow Direction Projection :' + prj_demcon)
                    self.textEdit.append('Flow Accumulation Projection :' + prj_flowacc)
                    sys.exit(1)
            
            
            prj_len = len(prj_demcon)
            self.textEdit.append(str(len(prj_demcon)))
            prj = prj_demcon
            if prj_len>0:
               self.textEdit.append('projection' + prj)               
            else:
               spatialRef = osr.SpatialReference()
               spatialRef.ImportFromEPSG(32633)
               projection = spatialRef.ExportToWkt()
               self.textEdit.append('import WGS84 utm 33N')
            
            originX = geotransform[0]#top left x
            originY = geotransform[3]#top left y
            #pixelWidth = geotransform[1]
            cellsize = geotransform[5]#pixel resolution

            self.textEdit.append( 'columns '+ str(cols))
            self.textEdit.append( 'rows '+ str(rows))
            self.textEdit.append( 'bands '+ str(bands))


            band_demcon=ds_demcon.GetRasterBand(1)
            inNoData=band_demcon.GetNoDataValue()
            demcon=band_demcon.ReadAsArray(0,0,cols,rows)
            

            band_demvoid=ds_demvoid.GetRasterBand(1)
            demvoid=band_demvoid.ReadAsArray(0,0,cols,rows)

            band_flowdir=ds_flowdir.GetRasterBand(1)
            flowdir=band_flowdir.ReadAsArray(0,0,cols,rows)

            band_flowacc=ds_flowacc.GetRasterBand(1)
            flowacc=band_flowacc.ReadAsArray(0,0,cols,rows)

            band_G=ds_G.GetRasterBand(1)
            G=band_G.ReadAsArray(0,0,cols,rows)

            G=numpy.deg2rad(G);# convert degree in radiant
            #nanvalue=-9999# nan  values
            nanvalue=inNoData
            # convert to float
            demcon=demcon.astype(numpy.float32)# convert to float
            demvoid=demvoid.astype(numpy.float32)# convert to float
            flowdir=flowdir.astype(numpy.float32)# convert to float
            flowacc=flowacc.astype(numpy.float32)# convert to float
            G=G.astype(numpy.float32)## convert to float

            G[demcon==nanvalue]=numpy.nan# in G change demcon=-9999 with nan      
            demvoid[demcon==nanvalue]=numpy.nan#change demcon=-9999 with nan
            flowdir[demcon==nanvalue]=numpy.nan#change demcon=-9999 with nan
            flowacc[demcon==nanvalue]=numpy.nan#change demcon=-9999 with nan
            demcon[demcon==nanvalue]=numpy.nan# change demcon=-9999 with nan

            dem=numpy.copy(demcon)# new array as demcon

            a=rows#a = rows number
            b=cols#b= columns number

            dem[numpy.isnan(dem)]=100000# nan is 100000

            self.textEdit.append( 'Calculating channel... ' )
            if channel_type=="channel_ASk":

                channel=numpy.zeros( (rows,cols))#matrix of rows and cols dimension

                tmp=flowacc*numpy.square(cellsize)*numpy.power( (G+0.0001), 1.7)# tmp= (FlowAcc.*cellsize^2.*(G+0.0001).^1.7<0.4*10^5)
                #c1=tmp<40000# true where matrix tmp è < 40000
                #c2=tmp>30000# true where matrix tmp è > 30000
                c1=tmp>100000
                #channel[c1*c2]=1# impose channel=1 where c1 and c2 are true
                channel[c1]=1# impose channel=1 where c1  true
                self.progressBar.setRange(1,rows-1)
                #self.progressBar.show()
                self.progressBar.activateWindow()
                for ctr in range(1,rows-2):# ctr counters rows
                    #self.textEdit.append( ctr)
                    self.progressBar.setValue(ctr)

                    for ctc in range(1,cols-2):# ctc counters colums
                        if channel[ctr,ctc]>0:# if channel > 0 in the row ctr e column ctc
                            x=ctr;
                            y=ctc;
                            while dem[x,y]<th_channel and x < a-2 and x>0 and y<b-2 and y>0:
                                fd=flowdir[x,y]
                                x=x+di[fd]
                                y=y+dj[fd]
                                if channel[x,y]==2:
                                    break
                                else:
                                    channel[x,y]=2# change channel to 2
            else:
                channel=self.from_flowacc_to_stream(flowacc,th_channel) # use threshold method

            channel[channel==nanvalue]=0 #chenge demcon=-9999 with 0
            dem=numpy.copy(demcon)# copy dem con in a new matrix
            self.textEdit.append( 'channel done' )
            H=numpy.zeros((rows,cols))# H is a matrix of rows, cols set to zero

            X=numpy.zeros((rows,cols))# X matrix with initial X coordinate
            Y=numpy.zeros((rows,cols))# Y matrix with initial Y coordinate
            MASK=numpy.zeros((rows,cols))# Mask matrix of visited cells
            ## if MASK(i,j) > 0 --> X(i,j),Y(i,j) initial x,y from i,j cells until CHANNEL(x,y)=1


            Ariver=numpy.zeros((rows,cols))# matrix of rows, cols set to zero

            self.textEdit.append( 'Calculating Ariver,H... ' )

            for ctr in range(1,rows-2):

                self.progressBar.setValue(ctr)
                for ctc in range(1,cols-2):
                    if dem[ctr,ctc]<10000:# if dem < 10000 in cell ctr,ctc
                        x=ctr
                        y=ctc
                        Ld=0
                        dm=0
                        sqrt2=numpy.sqrt(2)
                        while MASK[x,y]==0 and channel[x,y]==0 and x< a-2 and x>0 and y<b-2 and y>0 and   numpy.isnan(flowdir[x,y])==False :# il ciclo while finisce quando risalgo ad una cella con channel =1 oppure gia visitata (mask =1)
                            fd=flowdir[x,y]
                            x=x+di[fd]# set cell x,y on the basis of flowdir coding (user chosen)
                            y=y+dj[fd]

                        if  MASK[x,y]==1: # if cell was visited
                               Ariver[ctr,ctc]=flowacc[int(X[x,y]),int(Y[x,y])]# use X[x,y],Y[x,y]
                               H[ctr,ctc]=dem[ctr,ctc]-dem[int(X[x,y]),int(Y[x,y])]
                               MASK[ctr,ctc]=1# # if cell was visited, MASK = 1
                               X[ctr,ctc]=X[x,y]# value of x and y on the new cell
                               Y[ctr,ctc]=Y[x,y]
                        else:
                               Ariver[ctr,ctc]=flowacc[x,y]# ctr,ctc = floacc in cell x y
                               H[ctr,ctc]=dem[ctr,ctc]-dem[x,y]#ctr,ctc = dem in ctr,ctc - floacc in x y
                               MASK[ctr,ctc]=1
                               X[ctr,ctc]=x
                               Y[ctr,ctc]=y
            H[numpy.isnan(demvoid)]=numpy.nan#
            Ariver[ctr+1,:]=0;#last row = 0
            Ariver[:,ctc+1]=0;##last column = 0
            H[H==0]=0.00001# chenge H==0 to H=0.00001
            #%% Estimation of the water depth in the nearest element of the drainage network...
            #...connected to the location under exam using an hydraulic scaling relation (Leopold and maddock, 1953)
            #...h=aA^n
            self.textEdit.append( 'Ariver, H done' )
            #% Parameters estimated for the Ohio River basin:
            #a = 0.1035
           # n = 0.4057 (defoult)
            n=self.doubleSpinBoxN.value()
            self.textEdit.append('n ' +str(n))
            #% hydraulic scaling relation (A[km^2])
            self.textEdit.append( 'Calculating GFI... ' )
            #hr = a*numpy.power((((Ariver+1)*cellsize*cellsize)/1000000.0),n);# hr = a.*(((Ariver+1).*cellsize^2)./10^6).^n;
            hr = numpy.power((((Ariver+1)*cellsize*cellsize)/1000000.0),n);
            #%% Index GFI=ln[hr/H]
            hronH=hr/H;  #hr./H;
            ln_hronH=numpy.real(numpy.log(hronH));#ln_hronH=real(log(hronH));
            ## normalization
            idnan=numpy.isnan(ln_hronH)
            matrix_min=numpy.min(ln_hronH[idnan==False])
            matrix_max=numpy.max(ln_hronH[idnan==False])
            ln_hronH_norm= 2*( (ln_hronH-matrix_min)/(matrix_max-matrix_min)-0.5)
            id_hronH_norm_nan=numpy.where(numpy.isnan(ln_hronH_norm))
            nnan=len(id_hronH_norm_nan[0])
            self.textEdit.append( str(nnan) )
            if calibration==1:
                self.textEdit.append( 'Loading calibration files' )
                band_floodhazard=ds_floodhazard.GetRasterBand(1)
                floodhazard100=band_floodhazard.ReadAsArray(0,0,cols,rows)
                inNoData_floodhazard=band_floodhazard.GetNoDataValue()
                id_nan_floodhazard=numpy.where(floodhazard100==inNoData_floodhazard)
                floodhazard100[id_nan_floodhazard]=-9999

				
                self.textEdit.append( 'Starting calibration' )
                starttime_calib=time.time()

                ##Identification of the Drainage Network
                dem=numpy.copy(demcon)
                dem[numpy.isnan(dem)]=100000
                a=rows
                b=cols

                i1=(floodhazard100>0)*(floodhazard100<=1)
                i2=channel==2
                i3=numpy.multiply(i1,i2)
                channel[i3]=3

                for ctr in range(1,rows-2):
                    #print ctr
                    for ctc in range(1,cols-2):
                        if channel[ctr,ctc]==3:
                            x=ctr;
                            y=ctc;
                            while dem[x,y]<10000 and x < rows-2 and x>0 and y<cols-2 and y>0:
                                fd=flowdir[x,y]
                                x=x+di[fd]
                                y=y+dj[fd]
                                if channel[x,y]==3:
                                    break
                                else:
                                    channel[x,y]=3


                # MARGINAL HAZARD
                #RischioAdB=numpy.zeros((rows,cols))
                #RischioAdB[(floodhazard100>0)*(floodhazard100<=1)]=1
                #if debug==1:
                #    numpy.savetxt('RischioAdB_python.txt',RischioAdB)
                #Marginal_hazard=numpy.zeros((rows,cols))
                Marginal_hazard=numpy.copy(channel)
                for ctr in range(1,rows-2):
                    #print ctr
                    for ctc in range(1,cols-2):
                        if dem[ctr,ctc]<100000:
                            x=ctr;
                            y=ctc;
                            while channel[x,y] == 0 and x < rows-2 and x>0 and y<cols-2 and y>0 and flowdir[x,y]>0:
                                fd=flowdir[x,y]
                                x=x+di[fd]
                                y=y+dj[fd]
                                Marginal_hazard[ctr,ctc] =channel[x,y]
                Marginal_hazard[ctr+1,:]=0;
                Marginal_hazard[:,ctc+1]=0;
                if debug==1:
                    file_output_int=file_output[0:len(file_output)-4]+'_Marginal_hazard_python.txt'
                    numpy.savetxt(file_output_int, Marginal_hazard)
                Marginal_hazard[ Marginal_hazard!=3]=0
                Marginal_hazard[ Marginal_hazard==3]=1

                CalibrationArea=Marginal_hazard
                CalibrationArea[numpy.isnan(ln_hronH_norm)]=0
                if debug==1:
                    file_output_int=file_output[0:len(file_output)-4]+'_CalibrationArea_python.txt'
                    numpy.savetxt(file_output_int, CalibrationArea)

                floodhazard=numpy.zeros((rows,cols));
                floodhazard[ (floodhazard100>0) ]=1;

                if debug==1:
                    file_output_int=file_output[0:len(file_output)-4]+'_floodhazard_python.txt'
                    numpy.savetxt(file_output_int, floodhazard)
                ### make 1d vectors
                CalibrationArea_1d=CalibrationArea.flatten()
                floodhazard_1d=floodhazard.flatten()
                floodhazard100_1d=floodhazard100.flatten()
                ln_hronH_norm_1d=ln_hronH_norm.flatten()
                ### select the calibration area
                floodhazard_1d=floodhazard_1d[CalibrationArea_1d==1]
                ln_hronH_norm_1d=ln_hronH_norm_1d[CalibrationArea_1d==1]

                ### calculate the roc curve
               # fpr, tpr, thresholds = metrics.roc_curve(floodhazard_1d, ln_hronH_norm_1d)
                fpr, tpr, thresholds =self.my_roc_curve(floodhazard_1d, ln_hronH_norm_1d,0.001)
             
                #roc_auc = metrics.auc(fpr, tpr)
                roc_auc=self.my_auc(fpr,tpr)
                n_th=len(thresholds)
                
                ## optimal threshold
                F=fpr+(1-tpr)
                t_ln_hronH=thresholds[F.argmin()]
                step_th=t_ln_hronH-thresholds[F.argmin()-1]
                endtime_calib=time.time()
                self.textEdit.append('calibration completed in ' + str(endtime_calib-starttime_calib))
                self.textEdit.append('area under roc: ' +str(roc_auc))
                
                self.textEdit.append('optimal threshold: ' +str(t_ln_hronH))
                self.textEdit.append('n threshold: ' +str(n_th))
                self.textEdit.append('step threshold: ' +str(step_th))
                file_output_txt=file_output[0:len(file_output)-4]+'_performance.txt'
                tmp=numpy.zeros((6,1));
                tmp[0,0]=t_ln_hronH
                tmp[1,0]=fpr[F.argmin()]
                tmp[2,0]=tpr[F.argmin()]
                tmp[3,0]=F[F.argmin()]
                tmp[4,0]=roc_auc
                id_calibration=numpy.where(CalibrationArea_1d>0)
                n_calibration=len(id_calibration[0])
                n_tot=len(CalibrationArea_1d)
                p_calibration=numpy.double(n_calibration)/numpy.double(n_tot)
                tmp[5,0]=p_calibration*100

                string_out='Threshold= '+ str(tmp[0,0])+'\nRfp= '+str(tmp[1,0])+'\nRtp= '+str(tmp[2,0])+'\nobject function='+str(tmp[3,0])+'\nauc='+str(tmp[4,0])+'\ncalibration area='+str(tmp[5,0])
                text_file = open(file_output_txt, "w")
                text_file.write(string_out)
                text_file.close()
                 #numpy.savetxt(file_output_txt,tmp)
            # ln_hronH threshold

            ln_hronHbin=numpy.zeros((rows,cols))
            ln_hronHbin[ln_hronH_norm>t_ln_hronH]=1;
            ln_hronHbin[numpy.isnan(demcon)]=0;
            # ripulitura pixel isolati
            
            s = [[1,1,1],[1,1,1],[1,1,1]]
            label_im, nb_labels = ndimage.label(ln_hronHbin, structure=s)
            
            sizes = ndimage.sum(ln_hronHbin, label_im, range(nb_labels + 1))
            
            id_ok=numpy.where(sizes>=8)[0]
            n_ok=len(id_ok)
            ln_hronHbin=numpy.zeros((rows,cols))
            for ct_ok in range(0,n_ok):
                    print(ct_ok)
                    ct_label=id_ok[ct_ok]
                    print(ct_label)
                    ln_hronHbin[label_im==ct_label] = 1

           
            
            self.textEdit.append( 'GFI done' )

            if self.checkBoxWD. isChecked():
                 dem=numpy.copy(demcon)
                 sizew=2;#window size in pixel
                ### dem filtering
                 for ctr in range(sizew,rows-sizew):# ctr contatore delle righe
                    #print ctr
                    for ctc in range(sizew,cols-sizew):# ctc contatore delle colonne
                        if channel[ctr,ctc]>0:
                            window=dem[ctr-sizew:ctr+sizew,ctc-sizew:ctc+sizew];
                            dem[ctr,ctc]=numpy.mean(window)
                # wdepth sul channel
                 GFIwd=numpy.zeros((rows,cols))
                 for ctr in range(2,rows-1):
                    #print ctr
                    for ctc in range(2,cols-2):
                        if ln_hronHbin[ctr,ctc]>0:
                            x=ctr;
                            y=ctc;
                            while channel[x,y] == 0 and x < rows-2 and x>0 and y<cols-2 and y>0 and numpy.isnan(flowdir[x,y])==False:
                                fd=flowdir[x,y]
                                x=x+di[fd]
                                y=y+dj[fd]
                                if channel[x,y]>0:
                                    if GFIwd[x,y]< H[ctr,ctc]:
                                       GFIwd[x,y]= H[ctr,ctc]
                ## wdepth sulle altre celle
                 for ctr in range(2,rows-1):
                    #print ctr
                    for ctc in range(2,cols-2):
                        if ln_hronHbin[ctr,ctc]>0:
                            x=ctr;
                            y=ctc;
                            while channel[x,y] == 0 and x < rows-2 and x>0 and y<cols-2 and y>0 and numpy.isnan(flowdir[x,y])==False:
                                fd=flowdir[x,y]
                                x=x+di[fd]
                                y=y+dj[fd]
                                if channel[x,y]>0:
                                       GFIwd[ctr,ctc]= GFIwd[x,y]-H[ctr,ctc]
                 if debug ==1:
                     file_output_int=file_output[0:len(file_output)-4]+'_GFIwd_python.txt'
                     numpy.savetxt(file_output_int,GFIwd)
            if debug ==1:
                file_output_int_base=file_output[0:len(file_output)-4]
                numpy.savetxt(file_output_int_base+'_G_python.txt',G)
                self.writeOutputGeoTiff(G, geotransform,prj, rows, cols, file_output_int_base+'_G_python.tif')
                numpy.savetxt(file_output_int_base+'_Ariver_python.txt',Ariver)
                numpy.savetxt(file_output_int_base+'_H_python.txt',H)
                numpy.savetxt(file_output_int_base+'_channel_python.txt',channel)
                numpy.savetxt(file_output_int_base+'_hr_python.txt',hr)
                numpy.savetxt(file_output_int_base+'_hronH_python.txt',hronH)
                numpy.savetxt(file_output_int_base+'_ln_hronH_python.txt',ln_hronH)
                numpy.savetxt(file_output_int_base+'_ln_hronHbin_python.txt',ln_hronHbin)
                numpy.savetxt(file_output_int_base+'_ln_hronH_norm_python.txt',ln_hronH_norm)
                numpy.savetxt(file_output_int_base+'_ln_hronH_norm_1d_python.txt',ln_hronH_norm_1d)
            self.writeOutputGeoTiff(ln_hronH, geotransform,prj, rows, cols, file_output)
            file_outputbin=self.lineOutputBin.text()
            if len(file_outputbin)==0:
                file_outputbin=file_output[0:len(file_output)-4]+'bin.tif'

            if self.checkBoxFM.isChecked():
                self.writeOutputGeoTiff(ln_hronHbin, geotransform,prj, rows, cols, file_outputbin)


            file_outputWD=self.lineOutputWD.text()
            if len(file_outputWD)==0:
                file_outputWD=file_output[0:len(file_output)-4]+'wd.tif'
            if self.checkBoxWD.isChecked():
                self.writeOutputGeoTiff(GFIwd, geotransform,prj, rows, cols, file_outputWD)

            file_output_norm=file_output[0:len(file_output)-4]+'norm.tif'
            self.writeOutputGeoTiff(ln_hronH_norm, geotransform,prj, rows, cols, file_output_norm)

            self.progressBar.setValue(self.progressBar.maximum())
            self.textEdit.append( 'completed ')

            if self.checkBoxAdd. isChecked():
              self.loadOutputFile(file_output)
              if self.checkBoxFM.isChecked():
                    self.loadOutputFile(file_outputbin)
              if self.checkBoxWD.isChecked():
                    self.loadOutputFile(file_outputWD)
            band_target=None
            target_ds=None

            band_targetbin=None
            target_dsbin=None

            band_demcon=None
            ds_demcon=None

            band_demvoid=None
            ds_demvoid=None

            band_flowdir=None
            ds_flowdir=None

            band_flowacc=None
            ds_flowacc=None

            band_G=None
            ds_G=None
            #######################################################

            shutil.rmtree(temp_path)







