# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'nullDimWindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_NullDimWindow(object):
    def setupUi(self, NullDimWindow):
        NullDimWindow.setObjectName("NullDimWindow")
        NullDimWindow.resize(947, 473)
        self.centralwidget = QtWidgets.QWidget(NullDimWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setMaximumSize(QtCore.QSize(341, 16))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_3.addWidget(self.label_3)
        spacerItem = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.verticalLayout_3.addItem(spacerItem)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setTitle("Essential Parameters")
        self.groupBox.setAlignment(QtCore.Qt.AlignCenter)
        self.groupBox.setFlat(False)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.browseFile = QtWidgets.QPushButton(self.groupBox)
        self.browseFile.setObjectName("browseFile")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.browseFile)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.cellManager = QtWidgets.QPushButton(self.groupBox)
        self.cellManager.setObjectName("cellManager")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.cellManager)
        self.labelFileName = QtWidgets.QLabel(self.groupBox)
        self.labelFileName.setObjectName("labelFileName")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.SpanningRole, self.labelFileName)
        self.showCurrentCell = QtWidgets.QPushButton(self.groupBox)
        self.showCurrentCell.setObjectName("showCurrentCell")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.SpanningRole, self.showCurrentCell)
        self.gridLayout.addLayout(self.formLayout_2, 0, 0, 1, 4)
        self.label_8 = QtWidgets.QLabel(self.groupBox)
        self.label_8.setTextFormat(QtCore.Qt.RichText)
        self.label_8.setObjectName("label_8")
        self.gridLayout.addWidget(self.label_8, 8, 2, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.groupBox)
        self.label_12.setTextFormat(QtCore.Qt.RichText)
        self.label_12.setObjectName("label_12")
        self.gridLayout.addWidget(self.label_12, 12, 2, 1, 1)
        self.input_x_H2O = QtWidgets.QLineEdit(self.groupBox)
        self.input_x_H2O.setAutoFillBackground(False)
        self.input_x_H2O.setClearButtonEnabled(False)
        self.input_x_H2O.setObjectName("input_x_H2O")
        self.gridLayout.addWidget(self.input_x_H2O, 9, 3, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 11, 0, 1, 1)
        self.input_pressure = QtWidgets.QLineEdit(self.groupBox)
        self.input_pressure.setAutoFillBackground(False)
        self.input_pressure.setText("")
        self.input_pressure.setClearButtonEnabled(False)
        self.input_pressure.setObjectName("input_pressure")
        self.gridLayout.addWidget(self.input_pressure, 3, 1, 1, 3)
        self.label_9 = QtWidgets.QLabel(self.groupBox)
        self.label_9.setTextFormat(QtCore.Qt.RichText)
        self.label_9.setObjectName("label_9")
        self.gridLayout.addWidget(self.label_9, 9, 2, 1, 1)
        self.input_x_O2 = QtWidgets.QLineEdit(self.groupBox)
        self.input_x_O2.setAutoFillBackground(False)
        self.input_x_O2.setClearButtonEnabled(False)
        self.input_x_O2.setObjectName("input_x_O2")
        self.gridLayout.addWidget(self.input_x_O2, 11, 3, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 5, 0, 1, 1)
        self.input_x_H2 = QtWidgets.QLineEdit(self.groupBox)
        self.input_x_H2.setAutoFillBackground(False)
        self.input_x_H2.setClearButtonEnabled(False)
        self.input_x_H2.setObjectName("input_x_H2")
        self.gridLayout.addWidget(self.input_x_H2, 8, 3, 1, 1)
        self.line_2 = QtWidgets.QFrame(self.groupBox)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout.addWidget(self.line_2, 10, 0, 1, 4)
        self.label_11 = QtWidgets.QLabel(self.groupBox)
        self.label_11.setTextFormat(QtCore.Qt.RichText)
        self.label_11.setObjectName("label_11")
        self.gridLayout.addWidget(self.label_11, 11, 2, 1, 1)
        self.input_x_N2 = QtWidgets.QLineEdit(self.groupBox)
        self.input_x_N2.setAutoFillBackground(False)
        self.input_x_N2.setClearButtonEnabled(False)
        self.input_x_N2.setObjectName("input_x_N2")
        self.gridLayout.addWidget(self.input_x_N2, 12, 3, 1, 1)
        self.input_molar_flowrate_fuel = QtWidgets.QLineEdit(self.groupBox)
        self.input_molar_flowrate_fuel.setAutoFillBackground(False)
        self.input_molar_flowrate_fuel.setText("")
        self.input_molar_flowrate_fuel.setClearButtonEnabled(False)
        self.input_molar_flowrate_fuel.setObjectName("input_molar_flowrate_fuel")
        self.gridLayout.addWidget(self.input_molar_flowrate_fuel, 4, 1, 1, 3)
        self.line = QtWidgets.QFrame(self.groupBox)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 6, 0, 1, 4)
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 3, 0, 1, 1)
        self.input_molar_flowrate_oxygen = QtWidgets.QLineEdit(self.groupBox)
        self.input_molar_flowrate_oxygen.setAutoFillBackground(False)
        self.input_molar_flowrate_oxygen.setText("")
        self.input_molar_flowrate_oxygen.setClearButtonEnabled(False)
        self.input_molar_flowrate_oxygen.setObjectName("input_molar_flowrate_oxygen")
        self.gridLayout.addWidget(self.input_molar_flowrate_oxygen, 5, 1, 1, 3)
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 7, 0, 2, 1)
        self.gridLayout_3.addWidget(self.groupBox, 0, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setText("")
        self.label_15.setObjectName("label_15")
        self.gridLayout_3.addWidget(self.label_15, 1, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_3, 0, 0, 3, 1)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_4.addItem(spacerItem1, 0, 1, 1, 1)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_13 = QtWidgets.QLabel(self.groupBox_2)
        self.label_13.setObjectName("label_13")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_13)
        self.input_meanTemperature = QtWidgets.QLineEdit(self.groupBox_2)
        self.input_meanTemperature.setAutoFillBackground(False)
        self.input_meanTemperature.setClearButtonEnabled(False)
        self.input_meanTemperature.setObjectName("input_meanTemperature")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.input_meanTemperature)
        self.verticalLayout.addLayout(self.formLayout)
        self.label_14 = QtWidgets.QLabel(self.groupBox_2)
        self.label_14.setTextFormat(QtCore.Qt.RichText)
        self.label_14.setObjectName("label_14")
        self.verticalLayout.addWidget(self.label_14)
        self.input_jList = QtWidgets.QLineEdit(self.groupBox_2)
        self.input_jList.setAutoFillBackground(False)
        self.input_jList.setClearButtonEnabled(False)
        self.input_jList.setObjectName("input_jList")
        self.verticalLayout.addWidget(self.input_jList)
        self.gridLayout_2.addWidget(self.groupBox_2, 0, 0, 1, 1)
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_16 = QtWidgets.QLabel(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setAlignment(QtCore.Qt.AlignCenter)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout.addWidget(self.label_16)
        self.input_experimentName = QtWidgets.QLineEdit(self.groupBox_3)
        self.input_experimentName.setAutoFillBackground(False)
        self.input_experimentName.setClearButtonEnabled(False)
        self.input_experimentName.setObjectName("input_experimentName")
        self.horizontalLayout.addWidget(self.input_experimentName)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.label_17 = QtWidgets.QLabel(self.groupBox_3)
        self.label_17.setTextFormat(QtCore.Qt.RichText)
        self.label_17.setObjectName("label_17")
        self.verticalLayout_2.addWidget(self.label_17)
        self.input_jExperiment = QtWidgets.QLineEdit(self.groupBox_3)
        self.input_jExperiment.setAutoFillBackground(False)
        self.input_jExperiment.setClearButtonEnabled(False)
        self.input_jExperiment.setObjectName("input_jExperiment")
        self.verticalLayout_2.addWidget(self.input_jExperiment)
        self.label_18 = QtWidgets.QLabel(self.groupBox_3)
        self.label_18.setTextFormat(QtCore.Qt.RichText)
        self.label_18.setObjectName("label_18")
        self.verticalLayout_2.addWidget(self.label_18)
        self.input_vExperiment = QtWidgets.QLineEdit(self.groupBox_3)
        self.input_vExperiment.setAutoFillBackground(False)
        self.input_vExperiment.setClearButtonEnabled(False)
        self.input_vExperiment.setObjectName("input_vExperiment")
        self.verticalLayout_2.addWidget(self.input_vExperiment)
        self.gridLayout_2.addWidget(self.groupBox_3, 0, 1, 1, 1)
        self.generateCurve = QtWidgets.QPushButton(self.centralwidget)
        self.generateCurve.setObjectName("generateCurve")
        self.gridLayout_2.addWidget(self.generateCurve, 1, 0, 1, 2)
        self.gridLayout_4.addLayout(self.gridLayout_2, 0, 2, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem2, 1, 2, 1, 1)
        self.logText = QtWidgets.QTextBrowser(self.centralwidget)
        self.logText.setObjectName("logText")
        self.gridLayout_4.addWidget(self.logText, 2, 2, 1, 1)
        self.gridLayout_4.setColumnStretch(0, 5)
        self.gridLayout_4.setColumnStretch(1, 1)
        self.gridLayout_4.setColumnStretch(2, 10)
        self.verticalLayout_4.addLayout(self.gridLayout_4)
        NullDimWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(NullDimWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 947, 21))
        self.menubar.setObjectName("menubar")
        NullDimWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(NullDimWindow)
        self.statusbar.setObjectName("statusbar")
        NullDimWindow.setStatusBar(self.statusbar)

        self.retranslateUi(NullDimWindow)
        QtCore.QMetaObject.connectSlotsByName(NullDimWindow)

    def retranslateUi(self, NullDimWindow):
        _translate = QtCore.QCoreApplication.translate
        NullDimWindow.setWindowTitle(_translate("NullDimWindow", "MainWindow"))
        self.label_3.setText(_translate("NullDimWindow", "Isothermal rSOC 0-D Model (Hydrogen Cell)"))
        self.label.setText(_translate("NullDimWindow", "Select existing Cell from file"))
        self.browseFile.setText(_translate("NullDimWindow", "Browse file"))
        self.label_2.setText(_translate("NullDimWindow", "Create/Edit Cells"))
        self.cellManager.setText(_translate("NullDimWindow", "Cell Manager"))
        self.labelFileName.setText(_translate("NullDimWindow", "Cell File: Not selected"))
        self.showCurrentCell.setText(_translate("NullDimWindow", "Show currently used parameters"))
        self.label_8.setText(_translate("NullDimWindow", "<html><head/><body><p>H<span style=\" vertical-align:sub;\">2</span></p></body></html>"))
        self.label_12.setText(_translate("NullDimWindow", "<html><head/><body><p>N<span style=\" vertical-align:sub;\">2</span></p></body></html>"))
        self.label_4.setText(_translate("NullDimWindow", "Molar flow rate Fuel (mol/s)"))
        self.label_10.setText(_translate("NullDimWindow", "Molar composition of Oxidant"))
        self.label_9.setText(_translate("NullDimWindow", "<html><head/><body><p>H<span style=\" vertical-align:sub;\">2</span>O</p></body></html>"))
        self.label_6.setText(_translate("NullDimWindow", "Molar flow rate Oxidant (mol/s)"))
        self.label_11.setText(_translate("NullDimWindow", "<html><head/><body><p>O<span style=\" vertical-align:sub;\">2</span></p></body></html>"))
        self.label_5.setText(_translate("NullDimWindow", "Operation Pressure (bar)"))
        self.label_7.setText(_translate("NullDimWindow", "Molar composition of Fuel"))
        self.groupBox_2.setTitle(_translate("NullDimWindow", "Voltage-Current Curve and Performance"))
        self.label_13.setText(_translate("NullDimWindow", "Mean operation temperature (K)"))
        self.label_14.setText(_translate("NullDimWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Current density list (A/cm<span style=\" vertical-align:super;\">2</span>). Example: 0, 0.1, 0.2, 0.5, 1.3</p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">OR  </span> 0 : 2.0 : 10 (start : end : total points)</p></body></html>"))
        self.groupBox_3.setTitle(_translate("NullDimWindow", "Compare with Experimental Data"))
        self.label_16.setText(_translate("NullDimWindow", "Label"))
        self.label_17.setText(_translate("NullDimWindow", "<html><head/><body><p>Current density list (A/cm<span style=\" vertical-align:super;\">2</span>): 1, 2, 3, 4, ...</p></body></html>"))
        self.label_18.setText(_translate("NullDimWindow", "<html><head/><body><p>Voltage list (V): 1, 2, 3, 4, ...</p></body></html>"))
        self.generateCurve.setText(_translate("NullDimWindow", "Generate Voltage-Current Curve"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    NullDimWindow = QtWidgets.QMainWindow()
    ui = Ui_NullDimWindow()
    ui.setupUi(NullDimWindow)
    NullDimWindow.show()
    sys.exit(app.exec_())
