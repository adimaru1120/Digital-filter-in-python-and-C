import os
import numpy as np
from pylab import *
import scipy.signal as signal
# graphs
import pyqtgraph as pg
# C code handling
import ctypes as ct
# GUI
from pyqtgraph.Qt import QtGui, QtCore
import levinson as levi

import sys
from functools import partial

os.chdir('D:/studia/magisterka/praca magisterska/program')
lib = ct.cdll.LoadLibrary('D:/studia/magisterka/praca magisterska/program/dsp/filtr.dll')

# processing params
fs = 44100

class FilterDesign(QtGui.QMainWindow):

    b = 1
    a = 1
    coefB = 1
    coefA = 1
    order = 0
    q = 15
    sectionsNo = 0
    sos = 1
    

    label1 = "Rectangular"
    label2 = "Floating point"
    label3 = "direct"
    
    ftype = 2 #FIR = 0, IIR =1
    fstructure = 0 #FIR: direct = 0, lattice =1; IIR: direct = 0, lattice = 1


    def __init__(self):
        super(FilterDesign, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("Filter Design")
        self.UIcomponents()

    def UIcomponents(self):

        self.computeButton = QtGui.QPushButton('Compute coefficients', self)
        self.shwButton = QtGui.QPushButton('Show Specification', self)
        self.importButton = QtGui.QPushButton('Import to C', self)
        
        self.computeButton.resize(110, 20)
        self.computeButton.move(130,20)

        self.shwButton.resize(100, 20)
        self.shwButton.move(245,20)
        
        self.importButton.resize(100, 20)
        self.importButton.move(350,20)


        self.computeButton.clicked.connect(self.Compute)
        self.shwButton.clicked.connect(self.ShwSpecif)
        self.importButton.clicked.connect(self.ImportFilter)

        self.responseTypeBox = pg.QtGui.QComboBox(self)
        self.responseTypeBox.addItem("Lowpass")
        self.responseTypeBox.addItem("Highpass")
        self.responseTypeBox.addItem("Bandpass")
        self.responseTypeBox.addItem("Bandstop")
        self.responseTypeBox.resize(70,20)
        self.responseTypeBox.move(50, 20)

        self.responseTypeBox.currentIndexChanged.connect(self.Characteristic_Change)

        self.responseBox = pg.QtGui.QComboBox(self)
        self.responseBox.addItem("FIR")
        self.responseBox.addItem("IIR")
        self.responseBox.resize(70,20)
        self.responseBox.move(50, 50)

        self.responseBox.currentIndexChanged.connect(self.Response_Change)

        self.typeFIRBox = pg.QtGui.QComboBox(self)
        #self.typeFIRBox.addItem("Equiripple")
        #self.typeFIRBox.addItem("Least-squared")
        self.typeFIRBox.addItem("Rectangular")
        self.typeFIRBox.addItem("Blackman")
        self.typeFIRBox.addItem("Hamming")
        self.typeFIRBox.addItem("Hann")
        self.typeFIRBox.addItem("Bartlett")
        self.typeFIRBox.resize(110, 20)
        self.typeFIRBox.move(130, 50)

        self.typeIIRBox = pg.QtGui.QComboBox(self)
        self.typeIIRBox.addItem("Butterworth")
        self.typeIIRBox.addItem("Chebyshev type 1")
        self.typeIIRBox.addItem("Chebyshev type 2")
        self.typeIIRBox.addItem("Elliptic")
        self.typeIIRBox.resize(120,20)
        self.typeIIRBox.move(245, 50)
        
        self.structureBox = pg.QtGui.QComboBox(self)
        self.structureBox.addItem("direct")
        self.structureBox.addItem("lattice")
        self.structureBox.addItem("biquad")
        self.structureBox.resize(70,20)
        self.structureBox.move(50, 80)
        
        self.coding = pg.QtGui.QComboBox(self)
        self.coding.addItem("Floating point")
        self.coding.addItem("Q15")
        self.coding.addItem("Q12")
        self.coding.addItem("Q10")
        self.coding.addItem("Q 8")
        self.coding.resize(120,20)
        self.coding.move(130, 80)
        

        #parametr 1 - FIR: order, IIR: Fpass1

        self.param1 = pg.QtGui.QLabel('Order:', self)
        self.param1.move(50, 115)

        self.param1Val = QtGui.QSpinBox(self)
        self.param1Val.setRange(1, 499)
        self.param1Val.resize(80, 20)
        self.param1Val.move(100, 120)
        
        #parametr 2 - FIR: -, IIR: Fpass2

        self.param2 = pg.QtGui.QLabel('Fpass2', self)
        self.param2.move(290, 115)

        self.param2Val = QtGui.QSpinBox(self)
        self.param2Val.setRange(1, 22050)
        self.param2Val.resize(80, 20)
        self.param2Val.move(350, 120)
        
        #parametr 3 - FIR: Fc1, IIR: Fstop1

        self.param3 = pg.QtGui.QLabel('Fc1:', self)
        self.param3.move(50, 155)

        self.param3Val = QtGui.QSpinBox(self)
        self.param3Val.setRange(1, 22050)
        self.param3Val.resize(80, 20)
        self.param3Val.move(100, 160)

        #parametr 4 - FIR: -, IIR: Fstop2
        
        self.param4 = pg.QtGui.QLabel('Fstop2', self)
        self.param4.move(290, 155)

        self.param4Val = QtGui.QSpinBox(self)
        self.param4Val.setRange(1, 22050)
        self.param4Val.resize(80, 20)
        self.param4Val.move(350, 160)

        #parametr 5, FIR: Fc2, IIR: Apass & parametr 6: FIR: -, IIR:Astop

        self.param5 = pg.QtGui.QLabel('Fc2:', self)
        self.param5.move(50, 195)

        self.param5Val = QtGui.QSpinBox(self)
        self.param5Val.setRange(1, 22050)
        self.param5Val.resize(80, 20)
        self.param5Val.move(100, 200)

        self.param6 = pg.QtGui.QLabel('Astop', self)
        self.param6.move(290, 195)

        self.param6Val = QtGui.QSpinBox(self)
        self.param6Val.setRange(1, 150)
        self.param6Val.resize(80, 20)
        self.param6Val.move(350, 200)

        #Init UIcomponents disable
        
        #bandpass/bandstop spinbox disable
        self.param2Val.setEnabled(False)
        self.param4Val.setEnabled(False)
        self.param5Val.setEnabled(False)
        self.param6Val.setEnabled(False)

        self.structureBox.setEnabled(False)
        
        self.param2.setHidden(True)
        self.param4.setHidden(True)
        self.param6.setHidden(True)
        
        
        #IIR windows comboboox disable 
        self.typeIIRBox.setEnabled(False)


        self.show()

    def Characteristic_Change(self):        #LP, HP, BP, BS
        if(self.responseTypeBox.currentText() == "Bandpass" or
           self.responseTypeBox.currentText() == "Bandstop"):
            
            if(self.responseBox.currentText() == "FIR"):
                self.param5Val.setEnabled(True)
            elif(self.responseBox.currentText() == "IIR"):
                self.param2Val.setEnabled(True)
                self.param4Val.setEnabled(True)

        else:
            if(self.responseBox.currentText() == "FIR"):
                self.param5Val.setEnabled(False)
                self.param2Val.setEnabled(False)
                self.param4Val.setEnabled(False)
            elif(self.responseBox.currentText() == "IIR"):
                self.param2Val.setEnabled(False)
                self.param4Val.setEnabled(False)
                


    def Response_Change(self):              # FIR/IIR
        if (self.responseBox.currentText() == "FIR"):
        
            self.typeFIRBox.setEnabled(True)
            self.typeIIRBox.setEnabled(False)
            self.param5Val.setEnabled(False)
            self.param6Val.setEnabled(False)
            self.structureBox.setEnabled(False)
            
            self.param1.setText("Order")
            self.param3.setText("Fc1")
            self.param5.setText("Fc2")
            
            self.param2.setHidden(True)
            self.param4.setHidden(True)
            self.param6.setHidden(True)
            
            self.param1Val.setMaximum(499)
            self.param5Val.setMaximum(22050)
        else:
            self.typeFIRBox.setEnabled(False)
            self.typeIIRBox.setEnabled(True)
            self.param6Val.setEnabled(True)
            self.param5Val.setEnabled(True)
            self.structureBox.setEnabled(True)
            
            self.param1.setText("Fpass1")
            self.param3.setText("Fstop1")
            self.param5.setText("Apass")
            
            self.param2.setHidden(False)
            self.param4.setHidden(False)
            self.param6.setHidden(False)
            
            self.param1Val.setMaximum(22050)
            self.param5Val.setMaximum(150)
            
        

    def Compute(self):
        if(self.responseBox.currentText() == "FIR"):
            if (self.responseTypeBox.currentText() == "Lowpass"):
                #if (self.typeFIRBox.currentText() == "Equiripple"):
                #    pass
                #elif (self.typeFIRBox.currentText() == "Least-squared"):
                #    pass
                if (self.typeFIRBox.currentText() == "Rectangular"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'boxcar')
                    self.a = 1                                                                                               
                elif (self.typeFIRBox.currentText() == "Blackman"):                                                                        
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'blackman')
                    self.a = 1                                                                                               
                elif (self.typeFIRBox.currentText() == "Hamming"):                                                                         
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'hamming')
                    self.a = 1                                                                                               
                elif (self.typeFIRBox.currentText() == "Hann"):                                                                            
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'hann')
                    self.a = 1                                                                                               
                elif (self.typeFIRBox.currentText() == "Bartlett"):                                                                        
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'bartlett')
                    self.a = 1
                
            elif (self.responseTypeBox.currentText() == "Highpass"):
                #if (self.typeFIRBox.currentText() == "Equiripple"):
                #    pass
                #elif (self.typeFIRBox.currentText() == "Least-squared"):
                #    pass
                if (self.typeFIRBox.currentText() == "Rectangular"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'boxcar', pass_zero=False)
                    self.a = 1
                elif (self.typeFIRBox.currentText() == "Blackman"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'blackman', pass_zero=False)
                    self.a = 1
                elif (self.typeFIRBox.currentText() == "Hamming"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'hamming', pass_zero=False)
                    self.a = 1
                elif (self.typeFIRBox.currentText() == "Hann"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'hann', pass_zero=False)
                    self.a = 1
                elif (self.typeFIRBox.currentText() == "Bartlett"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = self.param3Val.value(), fs = fs, window = 'bartlett', pass_zero=False)
                    self.a = 1
                
            elif (self.responseTypeBox.currentText() == "Bandpass"):
                #if (self.typeFIRBox.currentText() == "Equiripple"):
                #    pass
                #elif (self.typeFIRBox.currentText() == "Least-squared"):
                #    pass
                if (self.typeFIRBox.currentText() == "Rectangular"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'boxcar', pass_zero=False)
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Blackman"):                                                                       
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'blackman', pass_zero=False)
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Hamming"):                                                                        
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'hamming', pass_zero=False)
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Hann"):                                                                           
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'hann', pass_zero=False)
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Bartlett"):                                                                       
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'bartlett', pass_zero=False)
                    self.a = 1
                
            elif (self.responseTypeBox.currentText() == "Bandstop"):
                #if (self.typeFIRBox.currentText() == "Equiripple"):
                #    pass
                #elif (self.typeFIRBox.currentText() == "Least-squared"):
                #    pass
                if (self.typeFIRBox.currentText() == "Rectangular"):
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'boxcar')
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Blackman"):                                                                       
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'blackman')
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Hamming"):                                                                        
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'hamming')
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Hann"):                                                                           
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'hann')
                    self.a = 1                                                                                              
                elif (self.typeFIRBox.currentText() == "Bartlett"):                                                                       
                    self.b = signal.firwin(numtaps = self.param1Val.value() + 1, cutoff = [self.param3Val.value(), self.param5Val.value()], fs = fs, window = 'bartlett')
                    self.a = 1
            
            self.ftype = 0

            self.label3 = "direct"
            if(self.coding.currentText() == "Floating point"):
                self.coefB = self.b
                self.coefA = self.a
                self.fstructure = 0
            else:
                self.q = int(self.coding.currentText()[1:3])
                
                for i in range(len(self.b)):
                    self.coefB[i] = round(self.b[i] * (2**self.q))
                            
                self.fstructure = 1
                
            self.label1 = "FIR " + self.typeFIRBox.currentText() + " window"
            self.label2 = self.coding.currentText()
            
        elif(self.responseBox.currentText() == "IIR"):
            if (self.responseTypeBox.currentText() == "Lowpass"):
                if (self.typeIIRBox.currentText() == "Butterworth"):
                    self.order, Wn = signal.buttord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                    gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.butter(N = self.order, Wn = Wn, btype = 'low', fs = fs)
                    
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 1"):
                    self.order, Wn = signal.cheb1ord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby1(N = self.order,Wn = Wn, rp = self.param5Val.value(), btype = 'low', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 2"):
                    self.order, Wn = signal.cheb2ord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby2(N = self.order, rs = self.param6Val.value(), Wn = Wn, btype = 'low', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Elliptic"):
                    self.order, Wn = signal.ellipord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.ellip(N = self.order, rp = self.param5Val.value(), rs = self.param6Val.value(), Wn = Wn , btype = 'low', fs = fs)
                    
            elif (self.responseTypeBox.currentText() == "Highpass"):
                if (self.typeIIRBox.currentText() == "Butterworth"):
                    self.order, Wn = signal.buttord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                    gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.butter(N = self.order,Wn = Wn, btype = 'high', fs = fs)
                    
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 1"):
                    self.order, Wn = signal.cheb1ord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby1(N = self.order,Wn = Wn, rp = self.param5Val.value(), btype = 'high', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 2"):
                    self.order, Wn = signal.cheb2ord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby2(N = self.order, rs = self.param6Val.value(), Wn = Wn, btype = 'high', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Elliptic"):
                    self.order, Wn = signal.ellipord(wp = self.param1Val.value(), ws = self.param3Val.value(),
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.ellip(N = self.order, rp = self.param5Val.value(), rs = self.param6Val.value(),
                                                  Wn = Wn ,btype = 'high', fs = fs)

            elif (self.responseTypeBox.currentText() == "Bandpass"):
                if (self.typeIIRBox.currentText() == "Butterworth"):
                    self.order, Wn = signal.buttord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                    ws = [self.param3Val.value(), self.param4Val.value()],
                                                    gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.butter(N = self.order,Wn = Wn, btype = 'bandpass', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 1"):
                    self.order, Wn = signal.cheb1ord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby1(N = self.order,Wn = Wn, rp = self.param5Val.value(), btype = 'bandpass', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Chebyshev type 2"):
                    self.order, Wn = signal.cheb2ord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby2(N = self.order, rs = self.param6Val.value(), Wn = Wn, btype = 'bandpass', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Elliptic"):
                    self.order, Wn = signal.ellipord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.ellip(N = self.order, rp = self.param5Val.value(), rs = self.param6Val.value(),
                                                  Wn = Wn ,btype = 'bandpass', fs = fs)

            elif (self.responseTypeBox.currentText() == "Bandstop"):
                if(self.typeIIRBox.currentText() == "Butterworth"):
                    self.order, Wn = signal.buttord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                    ws = [self.param3Val.value(), self.param4Val.value()],
                                                    gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                                           
                    self.b, self.a = signal.butter(N = self.order, Wn = Wn, btype = 'bandstop', fs = fs)
                
                elif(self.typeIIRBox.currentText() == "Chebyshev type 1"):
                    self.order, Wn = signal.cheb1ord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                                            
                    self.b, self.a = signal.cheby1(N = self.order, Wn = Wn, rp = self.param5Val.value(), btype = 'bandstop', fs = fs)
                
                elif (self.typeIIRBox.currentText() == "Chebyshev type 2"):
                    self.order, Wn = signal.cheb2ord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.cheby2(N = self.order, rs = self.param6Val.value(), Wn = Wn, btype = 'bandstop', fs = fs)
                    
                elif (self.typeIIRBox.currentText() == "Elliptic"):
                    self.order, Wn = signal.ellipord(wp = [self.param1Val.value(), self.param2Val.value()],
                                                     ws = [self.param3Val.value(), self.param4Val.value()],
                                                     gpass = self.param5Val.value(), gstop = self.param6Val.value(), fs = fs)
                    
                    self.b, self.a = signal.ellip(N = self.order, rp = self.param5Val.value(), rs = self.param6Val.value(),
                                                  Wn = Wn ,btype = 'bandstop', fs = fs)
                                                       
            
            self.ftype = 1
            
            if(self.structureBox.currentText() == "direct"):
                    self.coefA = self.a
                    self.coefB = self.b
                    self.fstructure = 0
                
            elif(self.structureBox.currentText() == "lattice"):
                R, U, self.coefA, e = levi.rlevinson(self.a,1)
                self.coefA = -self.coefA
                V = np.zeros(len(self.a))
    
                for i in range(len(self.a)-1,-1, -1):
                    subterm = U@V
                    V[i] = self.b[i] - subterm[i]
                    
                self.coefB = V
                self.fstructure = 1
            
            elif(self.structureBox.currentText() == "biquad"):
                self.sos = signal.tf2sos(self.b, self.a, 'nearest')
                
                self.sectionsNo = self.sos.shape[0]
                self.fstructure = 4
 
            if(self.coding.currentText() != "Floating point"):                
                self.q = int(self.coding.currentText()[1:3])
                
                if(self.structureBox.currentText() == "biquad"):
                    sosI = np.zeros((self.sos.shape[0],6))
                    for i in range(self.sos.shape[0]):
                        for j in range(len(self.sos[0])):
                            sosI[i][j] = round(self.sos[i][j] * (2**self.q))/(2**self.q)
                            
                    self.sos = sosI
                    self.fstructure = 4
                    
                elif(self.structureBox.currentText() == "direct"):
                    for i in range(len(self.coefB)):
                        self.coefB[i] = round(self.coefB[i] * (2**self.q))
                    
                    for i in range(len(self.coefA)):
                        self.coefA[i] = round(self.coefA[i] * (2**self.q))
                    self.fstructure = 2
                    
                elif(self.structureBox.currentText() == "lattice"):
                    for i in range(len(self.coefB)):
                        self.coefB[i] = round(self.coefB[i] * (2**self.q))/(2**self.q)
                    
                    for i in range(len(self.coefA)):
                        self.coefA[i] = round(self.coefA[i] * (2**self.q))/(2**self.q)
                        
                    self.fstructure = 1
            
            self.label1 = self.typeIIRBox.currentText()
            self.label2 = self.coding.currentText()
            self.label3 = self.structureBox.currentText()

        if(self.structureBox.currentText() == "biquad"):
            print("sos: ", self.sos)
            print("sectionsNo", self.sectionsNo)
        else:
            print("Numerator: " ,self.coefB)
            print("Denominator: " ,self.coefA)
        print("strukture", self.fstructure)
        

    def ShwSpecif(self):
        #figure(1).clf()
        figure(2).clf()
        figure(3).clf()
        #plot frequency
        w,h = signal.freqz(self.b, self.a)
        h_dB = 20 * log10(abs(h))
        f = w * fs / (2 * np.pi)
        figure(1)
        plot(f,h_dB)
        ylim(-150, 5)
        ylabel('Magnitude (db)')
        xlabel(r'Frequency (Hz)')
        title(r'Magnitude response')
        
        
        #plot phase
        figure(2)
        h_phase = np.unwrap(np.angle(h))
        plot(f,h_phase)
        ylabel('Phase (radians)')
        xlabel(r'Frequency (Hz)')
        title(r'Phase response')
        
        #plot zero-poles
        z, p, k = signal.tf2zpk(self.b, self.a)
        figure(3)
        theta = linspace(-pi,pi,201)
        plot(cos(theta),sin(theta))
        scatter(real(z),imag(z), color = 'red', label = 'zero')
        scatter(real(p),imag(p), color = 'green', label = 'poles')
        legend()
        title(r'Zero-Poles')

        show()
        
    def ImportFilter(self):
        wsk_ftype = ct.c_int.in_dll(lib, 'ftype')
        wsk_ftype.value = self.ftype
        
        wsk_fstructure = ct.c_int.in_dll(lib, 'fstructure')
        wsk_fstructure.value = self.fstructure        
        
        if(self.coding.currentText() == "Floating point"):
            if(self.structureBox.currentText() == "biquad"):
                
                wsk_sectionNo = ct.c_int.in_dll(lib, 'sectionsNo')
                wsk_sectionNo.value = self.sectionsNo
                
                type_for_Bn = ct.c_float*6
                Bn = type_for_Bn()
                lib.importBiquad.argtype = ct.POINTER(ct.c_float), ct.c_int
                type_for_ci = ct.c_int
                ci= type_for_ci()
                
                for i in range(self.sectionsNo):
                    ci = i  
                    for j in range(len(self.sos[0])):
                        Bn[j] = self.sos[i][j]
                    
                    lib.importBiquad(Bn, ci)
                        
            else:
                wsk_lenB = ct.c_int.in_dll(lib, 'lenB')
                wsk_lenB.value = len(self.coefB)
                
                type_for_Bn = (ct.c_float * len(self.coefB))
                Bn = type_for_Bn()
            
                for i in range(len(self.coefB)):
                    Bn[i] = self.coefB[i]
                
                if(self.ftype == 0):
                    lib.importNumerator.argtypes = ct.POINTER(ct.c_float),ct.c_size_t
                    lib.importNumerator(Bn, len(Bn))    #wywolanie funcji z srodowiska C
                    
                else:
                    wsk_lenA = ct.c_int.in_dll(lib, 'lenA')
                    wsk_lenA.value = len(self.coefA)
                    
                    type_for_Bdn = (ct.c_float * len(self.coefA))
                    Bdn = type_for_Bdn()
                    
                    for i in range(len(self.coefA)):
                        Bdn[i] = self.coefA[i]
                    
                    lib.importNumerator.argtypes = ct.POINTER(ct.c_float),ct.c_size_t
                    lib.importNumerator(Bn, len(Bn))
                    
                    lib.importDenominator.argtypes = ct.POINTER(ct.c_float), ct.c_size_t
                    lib.importDenominator(Bdn, len(Bdn))
        else:
            wsk_q = ct.c_int.in_dll(lib, 'q')
            wsk_q.value = self.q
            
            if(self.structureBox.currentText() == "biquad"):
                wsk_sectionNo = ct.c_int.in_dll(lib, 'sectionsNo')
                wsk_sectionNo.value = self.sectionsNo
                
                type_for_Bn = ct.c_float*6
                Bn = type_for_Bn()
                lib.importBiquad.argtype = ct.POINTER(ct.c_float), ct.c_int
                type_for_ci = ct.c_int
                ci= type_for_ci()
                
                for i in range(self.sectionsNo):
                    ci = i  
                    for j in range(len(self.sos[0])):
                        Bn[j] = self.sos[i][j]
                    lib.importBiquad(Bn, ci)
                        
            elif(self.structureBox.currentText() == "lattice"):
                wsk_lenB = ct.c_int.in_dll(lib, 'lenB')
                wsk_lenB.value = len(self.coefB)
                
                type_for_Bn = (ct.c_float * len(self.coefB))
                Bn = type_for_Bn()
            
                for i in range(len(self.coefB)):
                    Bn[i] = self.coefB[i]
                
                wsk_lenA = ct.c_int.in_dll(lib, 'lenA')
                wsk_lenA.value = len(self.coefA)
                    
                type_for_Bdn = (ct.c_float * len(self.coefA))
                Bdn = type_for_Bdn()
                    
                for i in range(len(self.coefA)):
                    Bdn[i] = self.coefA[i]
                    
                lib.importNumerator.argtypes = ct.POINTER(ct.c_float),ct.c_size_t
                lib.importNumerator(Bn, len(Bn))
                    
                lib.importDenominator.argtypes = ct.POINTER(ct.c_float), ct.c_size_t
                lib.importDenominator(Bdn, len(Bdn))
            else:
                type_for_BnI = (ct.c_int * len(self.coefB))
                BnI = type_for_BnI()
                
                for i in range(len(self.coefB)):
                    BnI[i] = int(self.coefB[i])
                
                if(self.ftype == 0):
                    lib.importNumeratorI.argtypes = ct.POINTER(ct.c_int),ct.c_size_t
                    lib.importNumeratorI(BnI, len(BnI))    #wywolanie funcji z srodowiska C
                    
                else:
                    wsk_lenA = ct.c_int.in_dll(lib, 'lenA')
                    wsk_lenA.value = len(self.coefA)
                    
                    type_for_BdnI = (ct.c_int * len(self.coefA))
                    BdnI = type_for_BdnI()
                    
                    for i in range(len(self.coefA)):
                        BdnI[i] = int(self.coefA[i])
                    
                    lib.importNumeratorI.argtypes = ct.POINTER(ct.c_int),ct.c_size_t
                    lib.importNumeratorI(BnI, len(BnI))
                    
                    lib.importDenominatorI.argtypes = ct.POINTER(ct.c_int), ct.c_size_t
                    lib.importDenominatorI(BdnI, len(BdnI))
