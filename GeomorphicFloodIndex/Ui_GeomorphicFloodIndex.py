# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Ui_GeomorphicFloodIndex.ui'
#
# Created: Tue Mar 21 13:17:23 2017
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from qgis.PyQt import QtCore, QtGui, QtWidgets

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)

class Ui_GeomorphicFloodIndex(object):
    def setupUi(self, GeomorphicFloodIndex):
        GeomorphicFloodIndex.setObjectName(_fromUtf8("GeomorphicFloodIndex"))
        GeomorphicFloodIndex.resize(900, 700)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GeomorphicFloodIndex.sizePolicy().hasHeightForWidth())
        GeomorphicFloodIndex.setSizePolicy(sizePolicy)
        GeomorphicFloodIndex.setMinimumSize(QtCore.QSize(900, 700))
        GeomorphicFloodIndex.setMaximumSize(QtCore.QSize(900, 700))
        self.groupBox_2 = QtWidgets.QGroupBox(GeomorphicFloodIndex)
        self.groupBox_2.setGeometry(QtCore.QRect(9, 360, 881, 341))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.groupBox_4 = QtWidgets.QGroupBox(self.groupBox_2)
        self.groupBox_4.setGeometry(QtCore.QRect(480, 10, 201, 61))
        self.groupBox_4.setStyleSheet(_fromUtf8("border-color: rgb(4, 4, 4);"))
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.progressBar = QtWidgets.QProgressBar(self.groupBox_4)
        self.progressBar.setGeometry(QtCore.QRect(13, 20, 159, 25))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setAlignment(QtCore.Qt.AlignCenter)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.groupBox_3 = QtWidgets.QGroupBox(self.groupBox_2)
        self.groupBox_3.setGeometry(QtCore.QRect(0, 210, 361, 91))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.textEdit = QtWidgets.QTextEdit(self.groupBox_3)
        self.textEdit.setGeometry(QtCore.QRect(10, 20, 331, 61))
        self.textEdit.setObjectName(_fromUtf8("textEdit"))
        self.checkBoxWD = QtWidgets.QCheckBox(self.groupBox_2)
        self.checkBoxWD.setGeometry(QtCore.QRect(70, 120, 127, 17))
        self.checkBoxWD.setObjectName(_fromUtf8("checkBoxWD"))
        self.checkBoxFM = QtWidgets.QCheckBox(self.groupBox_2)
        self.checkBoxFM.setGeometry(QtCore.QRect(70, 70, 351, 17))
        self.checkBoxFM.setObjectName(_fromUtf8("checkBoxFM"))
        self.checkBoxDebug = QtWidgets.QCheckBox(self.groupBox_2)
        self.checkBoxDebug.setGeometry(QtCore.QRect(20, 190, 311, 17))
        self.checkBoxDebug.setObjectName(_fromUtf8("checkBoxDebug"))
        self.checkBoxAdd = QtWidgets.QCheckBox(self.groupBox_2)
        self.checkBoxAdd.setGeometry(QtCore.QRect(20, 170, 311, 17))
        self.checkBoxAdd.setObjectName(_fromUtf8("checkBoxAdd"))
        self.btnOutput = QtWidgets.QToolButton(self.groupBox_2)
        self.btnOutput.setGeometry(QtCore.QRect(20, 30, 42, 25))
        self.btnOutput.setMinimumSize(QtCore.QSize(0, 25))
        self.btnOutput.setMaximumSize(QtCore.QSize(16777215, 25))
        self.btnOutput.setObjectName(_fromUtf8("btnOutput"))
        self.lineOutput = QtWidgets.QLineEdit(self.groupBox_2)
        self.lineOutput.setGeometry(QtCore.QRect(70, 30, 401, 25))
        self.lineOutput.setMinimumSize(QtCore.QSize(280, 25))
        self.lineOutput.setObjectName(_fromUtf8("lineOutput"))
        self.lblSlope_7 = QtWidgets.QLabel(self.groupBox_2)
        self.lblSlope_7.setGeometry(QtCore.QRect(70, 10, 197, 19))
        self.lblSlope_7.setObjectName(_fromUtf8("lblSlope_7"))
        self.btnOutputBin = QtWidgets.QToolButton(self.groupBox_2)
        self.btnOutputBin.setGeometry(QtCore.QRect(20, 90, 42, 25))
        self.btnOutputBin.setMinimumSize(QtCore.QSize(0, 25))
        self.btnOutputBin.setMaximumSize(QtCore.QSize(16777215, 25))
        self.btnOutputBin.setObjectName(_fromUtf8("btnOutputBin"))
        self.lineOutputBin = QtWidgets.QLineEdit(self.groupBox_2)
        self.lineOutputBin.setGeometry(QtCore.QRect(70, 90, 401, 25))
        self.lineOutputBin.setMinimumSize(QtCore.QSize(280, 25))
        self.lineOutputBin.setObjectName(_fromUtf8("lineOutputBin"))
        self.buttonBox = QtWidgets.QDialogButtonBox(self.groupBox_2)
        self.buttonBox.setGeometry(QtCore.QRect(160, 300, 251, 27))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Help|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.lblSlope_3 = QtWidgets.QLabel(self.groupBox_2)
        self.lblSlope_3.setGeometry(QtCore.QRect(440, 210, 281, 16))
        self.lblSlope_3.setObjectName(_fromUtf8("lblSlope_3"))
        self.cmbSlopeType = QtWidgets.QComboBox(self.groupBox_2)
        self.cmbSlopeType.setGeometry(QtCore.QRect(540, 210, 85, 20))
        self.cmbSlopeType.setEditable(False)
        self.cmbSlopeType.setModelColumn(0)
        self.cmbSlopeType.setObjectName(_fromUtf8("cmbSlopeType"))
        self.btnOutputWD = QtWidgets.QToolButton(self.groupBox_2)
        self.btnOutputWD.setGeometry(QtCore.QRect(20, 140, 42, 25))
        self.btnOutputWD.setMinimumSize(QtCore.QSize(0, 25))
        self.btnOutputWD.setMaximumSize(QtCore.QSize(16777215, 25))
        self.btnOutputWD.setObjectName(_fromUtf8("btnOutputWD"))
        self.lineOutputWD = QtWidgets.QLineEdit(self.groupBox_2)
        self.lineOutputWD.setGeometry(QtCore.QRect(70, 140, 401, 25))
        self.lineOutputWD.setMinimumSize(QtCore.QSize(280, 25))
        self.lineOutputWD.setObjectName(_fromUtf8("lineOutputWD"))
        self.groupBox = QtWidgets.QGroupBox(GeomorphicFloodIndex)
        self.groupBox.setGeometry(QtCore.QRect(9, 30, 881, 161))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.cmbDemVoid = QtWidgets.QComboBox(self.groupBox)
        self.cmbDemVoid.setGeometry(QtCore.QRect(203, 55, 600, 25))
        self.cmbDemVoid.setMinimumSize(QtCore.QSize(280, 25))
        self.cmbDemVoid.setMaximumSize(QtCore.QSize(600, 25))
        self.cmbDemVoid.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cmbDemVoid.setObjectName(_fromUtf8("cmbDemVoid"))
        self.cmbDemCon = QtWidgets.QComboBox(self.groupBox)
        self.cmbDemCon.setGeometry(QtCore.QRect(203, 24, 600, 25))
        self.cmbDemCon.setMinimumSize(QtCore.QSize(280, 25))
        self.cmbDemCon.setMaximumSize(QtCore.QSize(600, 25))
        self.cmbDemCon.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cmbDemCon.setObjectName(_fromUtf8("cmbDemCon"))
        self.lblFlowDir = QtWidgets.QLabel(self.groupBox)
        self.lblFlowDir.setGeometry(QtCore.QRect(11, 88, 181, 16))
        self.lblFlowDir.setObjectName(_fromUtf8("lblFlowDir"))
        self.lblTRIGGER = QtWidgets.QLabel(self.groupBox)
        self.lblTRIGGER.setGeometry(QtCore.QRect(11, 57, 171, 16))
        self.lblTRIGGER.setObjectName(_fromUtf8("lblTRIGGER"))
        self.cmbFlowDir = QtWidgets.QComboBox(self.groupBox)
        self.cmbFlowDir.setGeometry(QtCore.QRect(203, 86, 600, 25))
        self.cmbFlowDir.setMinimumSize(QtCore.QSize(280, 25))
        self.cmbFlowDir.setMaximumSize(QtCore.QSize(600, 25))
        self.cmbFlowDir.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cmbFlowDir.setObjectName(_fromUtf8("cmbFlowDir"))
        self.lblFLOWACC = QtWidgets.QLabel(self.groupBox)
        self.lblFLOWACC.setGeometry(QtCore.QRect(11, 119, 181, 16))
        self.lblFLOWACC.setObjectName(_fromUtf8("lblFLOWACC"))
        self.lblDEMCON = QtWidgets.QLabel(self.groupBox)
        self.lblDEMCON.setGeometry(QtCore.QRect(11, 26, 100, 16))
        self.lblDEMCON.setMaximumSize(QtCore.QSize(100, 20))
        self.lblDEMCON.setAutoFillBackground(False)
        self.lblDEMCON.setMargin(1)
        self.lblDEMCON.setObjectName(_fromUtf8("lblDEMCON"))
        self.cmbFlowAcc = QtWidgets.QComboBox(self.groupBox)
        self.cmbFlowAcc.setGeometry(QtCore.QRect(203, 117, 600, 25))
        self.cmbFlowAcc.setMinimumSize(QtCore.QSize(280, 25))
        self.cmbFlowAcc.setMaximumSize(QtCore.QSize(600, 25))
        self.cmbFlowAcc.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cmbFlowAcc.setObjectName(_fromUtf8("cmbFlowAcc"))
        self.toolButtonDemCon = QtWidgets.QToolButton(self.groupBox)
        self.toolButtonDemCon.setGeometry(QtCore.QRect(840, 30, 25, 19))
        self.toolButtonDemCon.setObjectName(_fromUtf8("toolButtonDemCon"))
        self.toolButtonDemVoid = QtWidgets.QToolButton(self.groupBox)
        self.toolButtonDemVoid.setGeometry(QtCore.QRect(840, 60, 25, 19))
        self.toolButtonDemVoid.setObjectName(_fromUtf8("toolButtonDemVoid"))
        self.toolButtonFlowDir = QtWidgets.QToolButton(self.groupBox)
        self.toolButtonFlowDir.setGeometry(QtCore.QRect(840, 90, 25, 19))
        self.toolButtonFlowDir.setObjectName(_fromUtf8("toolButtonFlowDir"))
        self.toolButtonFlowAcc = QtWidgets.QToolButton(self.groupBox)
        self.toolButtonFlowAcc.setGeometry(QtCore.QRect(840, 120, 25, 19))
        self.toolButtonFlowAcc.setObjectName(_fromUtf8("toolButtonFlowAcc"))
        self.groupBox_5 = QtWidgets.QGroupBox(GeomorphicFloodIndex)
        self.groupBox_5.setGeometry(QtCore.QRect(10, 190, 411, 151))
        self.groupBox_5.setObjectName(_fromUtf8("groupBox_5"))
        self.layoutWidget = QtWidgets.QWidget(self.groupBox_5)
        self.layoutWidget.setGeometry(QtCore.QRect(290, 20, 99, 126))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.cmbFlowDirCoding = QtWidgets.QComboBox(self.layoutWidget)
        self.cmbFlowDirCoding.setEditable(False)
        self.cmbFlowDirCoding.setModelColumn(0)
        self.cmbFlowDirCoding.setObjectName(_fromUtf8("cmbFlowDirCoding"))
        self.verticalLayout.addWidget(self.cmbFlowDirCoding)
        self.cmbChannelType = QtWidgets.QComboBox(self.layoutWidget)
        self.cmbChannelType.setEditable(False)
        self.cmbChannelType.setModelColumn(0)
        self.cmbChannelType.setObjectName(_fromUtf8("cmbChannelType"))
        self.verticalLayout.addWidget(self.cmbChannelType)
        self.doubleSpinBoxThesholdChannel = QtWidgets.QDoubleSpinBox(self.layoutWidget)
        self.doubleSpinBoxThesholdChannel.setDecimals(2)
        self.doubleSpinBoxThesholdChannel.setMinimum(0.0)
        self.doubleSpinBoxThesholdChannel.setMaximum(1000000000.0)
        self.doubleSpinBoxThesholdChannel.setSingleStep(1e-05)
        self.doubleSpinBoxThesholdChannel.setProperty("value", 10000.0)
        self.doubleSpinBoxThesholdChannel.setObjectName(_fromUtf8("doubleSpinBoxThesholdChannel"))
        self.verticalLayout.addWidget(self.doubleSpinBoxThesholdChannel)
        self.doubleSpinBoxN = QtWidgets.QDoubleSpinBox(self.layoutWidget)
        self.doubleSpinBoxN.setDecimals(4)
        self.doubleSpinBoxN.setMinimum(0.0)
        self.doubleSpinBoxN.setMaximum(100000.0)
        self.doubleSpinBoxN.setSingleStep(1e-05)
        self.doubleSpinBoxN.setProperty("value", 0.4057)
        self.doubleSpinBoxN.setObjectName(_fromUtf8("doubleSpinBoxN"))
        self.verticalLayout.addWidget(self.doubleSpinBoxN)
        self.lblSlope_4 = QtWidgets.QLabel(self.groupBox_5)
        self.lblSlope_4.setGeometry(QtCore.QRect(10, 60, 281, 16))
        self.lblSlope_4.setObjectName(_fromUtf8("lblSlope_4"))
        self.lblSlope_5 = QtWidgets.QLabel(self.groupBox_5)
        self.lblSlope_5.setGeometry(QtCore.QRect(10, 90, 291, 16))
        self.lblSlope_5.setObjectName(_fromUtf8("lblSlope_5"))
        self.lblSlope_6 = QtWidgets.QLabel(self.groupBox_5)
        self.lblSlope_6.setGeometry(QtCore.QRect(11, 121, 291, 16))
        self.lblSlope_6.setObjectName(_fromUtf8("lblSlope_6"))
        self.lblSlope_2 = QtWidgets.QLabel(self.groupBox_5)
        self.lblSlope_2.setGeometry(QtCore.QRect(10, 30, 271, 16))
        self.lblSlope_2.setObjectName(_fromUtf8("lblSlope_2"))
        self.groupBox_6 = QtWidgets.QGroupBox(GeomorphicFloodIndex)
        self.groupBox_6.setGeometry(QtCore.QRect(430, 190, 461, 151))
        self.groupBox_6.setObjectName(_fromUtf8("groupBox_6"))
        self.doubleSpinBoxTheshold = QtWidgets.QDoubleSpinBox(self.groupBox_6)
        self.doubleSpinBoxTheshold.setGeometry(QtCore.QRect(400, 20, 51, 20))
        self.doubleSpinBoxTheshold.setDecimals(2)
        self.doubleSpinBoxTheshold.setMinimum(-5.0)
        self.doubleSpinBoxTheshold.setMaximum(5.0)
        self.doubleSpinBoxTheshold.setSingleStep(1e-05)
        self.doubleSpinBoxTheshold.setProperty("value", -0.53)
        self.doubleSpinBoxTheshold.setObjectName(_fromUtf8("doubleSpinBoxTheshold"))
        self.cmbSFM = QtWidgets.QComboBox(self.groupBox_6)
        self.cmbSFM.setGeometry(QtCore.QRect(10, 110, 381, 25))
        self.cmbSFM.setMinimumSize(QtCore.QSize(280, 25))
        self.cmbSFM.setMaximumSize(QtCore.QSize(600, 25))
        self.cmbSFM.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cmbSFM.setObjectName(_fromUtf8("cmbSFM"))
        self.toolButtonSFM = QtWidgets.QToolButton(self.groupBox_6)
        self.toolButtonSFM.setGeometry(QtCore.QRect(410, 110, 25, 19))
        self.toolButtonSFM.setObjectName(_fromUtf8("toolButtonSFM"))
        self.checkBoxCalibration = QtWidgets.QCheckBox(self.groupBox_6)
        self.checkBoxCalibration.setGeometry(QtCore.QRect(11, 73, 351, 17))
        self.checkBoxCalibration.setObjectName(_fromUtf8("checkBoxCalibration"))
        self.lblDEMCON_2 = QtWidgets.QLabel(self.groupBox_6)
        self.lblDEMCON_2.setGeometry(QtCore.QRect(100, 90, 261, 16))
        self.lblDEMCON_2.setAutoFillBackground(False)
        self.lblDEMCON_2.setMargin(1)
        self.lblDEMCON_2.setObjectName(_fromUtf8("lblDEMCON_2"))
        self.lblSlope = QtWidgets.QLabel(self.groupBox_6)
        self.lblSlope.setGeometry(QtCore.QRect(220, 20, 181, 16))
        self.lblSlope.setObjectName(_fromUtf8("lblSlope"))
        self.checkBoxManualSet = QtWidgets.QCheckBox(self.groupBox_6)
        self.checkBoxManualSet.setGeometry(QtCore.QRect(10, 20, 181, 17))
        self.checkBoxManualSet.setObjectName(_fromUtf8("checkBoxManualSet"))

        self.retranslateUi(GeomorphicFloodIndex)
        QtCore.QMetaObject.connectSlotsByName(GeomorphicFloodIndex)
        GeomorphicFloodIndex.setTabOrder(self.cmbDemCon, self.cmbFlowDir)
        GeomorphicFloodIndex.setTabOrder(self.cmbFlowDir, self.btnOutput)
        GeomorphicFloodIndex.setTabOrder(self.btnOutput, self.buttonBox)
        GeomorphicFloodIndex.setTabOrder(self.buttonBox, self.textEdit)
        GeomorphicFloodIndex.setTabOrder(self.textEdit, self.lineOutput)

    def retranslateUi(self, GeomorphicFloodIndex):
        GeomorphicFloodIndex.setWindowTitle(_translate("GeomorphicFloodIndex", "GeomorphicFloodArea", None))
        self.groupBox_2.setTitle(_translate("GeomorphicFloodIndex", "Output ", None))
        self.groupBox_4.setTitle(_translate("GeomorphicFloodIndex", "Progress...", None))
        self.groupBox_3.setTitle(_translate("GeomorphicFloodIndex", "Log", None))
        self.checkBoxWD.setText(_translate("GeomorphicFloodIndex", "GFI waterdepth", None))
        self.checkBoxFM.setText(_translate("GeomorphicFloodIndex", "GFI derived flood-prone areas map", None))
        self.checkBoxDebug.setText(_translate("GeomorphicFloodIndex", "Create intermediate files", None))
        self.checkBoxAdd.setText(_translate("GeomorphicFloodIndex", "Add results to canvas", None))
        self.btnOutput.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.lblSlope_7.setText(_translate("GeomorphicFloodIndex", "GFI raster", None))
        self.btnOutputBin.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.lblSlope_3.setText(_translate("GeomorphicFloodIndex", "SLOPE TYPE", None))
        self.btnOutputWD.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.groupBox.setTitle(_translate("GeomorphicFloodIndex", "Input ", None))
        self.cmbDemVoid.setToolTip(_translate("GeomorphicFloodIndex", "parametro Cv", None))
        self.cmbDemCon.setToolTip(_translate("GeomorphicFloodIndex", "Parametro a", None))
        self.lblFlowDir.setText(_translate("GeomorphicFloodIndex", "Flow direction D8", None))
        self.lblTRIGGER.setText(_translate("GeomorphicFloodIndex", "Filled DEM", None))
        self.cmbFlowDir.setToolTip(_translate("GeomorphicFloodIndex", "Parametro n", None))
        self.lblFLOWACC.setText(_translate("GeomorphicFloodIndex", "Flow accumulation", None))
        self.lblDEMCON.setText(_translate("GeomorphicFloodIndex", "DEM", None))
        self.cmbFlowAcc.setToolTip(_translate("GeomorphicFloodIndex", "parametro Cv", None))
        self.toolButtonDemCon.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.toolButtonDemVoid.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.toolButtonFlowDir.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.toolButtonFlowAcc.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.groupBox_5.setTitle(_translate("GeomorphicFloodIndex", "Set Methodology Options", None))
        self.lblSlope_4.setText(_translate("GeomorphicFloodIndex", "Drainage network identification method", None))
        self.lblSlope_5.setText(_translate("GeomorphicFloodIndex", "Drainage network identification threshold", None))
        self.lblSlope_6.setText(_translate("GeomorphicFloodIndex", "Hydraulic scaling relation exponent", None))
        self.lblSlope_2.setText(_translate("GeomorphicFloodIndex", "Flow direction coding", None))
        self.groupBox_6.setTitle(_translate("GeomorphicFloodIndex", "Set Calibration Options", None))
        self.cmbSFM.setToolTip(_translate("GeomorphicFloodIndex", "Parametro a", None))
        self.toolButtonSFM.setText(_translate("GeomorphicFloodIndex", "...", None))
        self.checkBoxCalibration.setText(_translate("GeomorphicFloodIndex", "Calibrate Threshold", None))
        self.lblDEMCON_2.setText(_translate("GeomorphicFloodIndex", "Standard flood hazard map", None))
        self.lblSlope.setText(_translate("GeomorphicFloodIndex", "GFI classifier threshold", None))
        self.checkBoxManualSet.setText(_translate("GeomorphicFloodIndex", "Manually set threshold", None))

