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
        self.typeFIRBox.addItem("Boxcar")
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
        self.coding.addItem("floating point")
        self.coding.addItem("Q15")
        self.coding.addItem("Q12")
        self.coding.addItem("Q10")
        self.coding.addItem("Q 8")
        self.coding.resize(120,20)
        self.coding.move(130, 80)
        

        #parametr 1 - FIR: order, IIR: Fpass1

        self.Order_Fpass1 = pg.QtGui.QLabel('Order:', self)
        self.Order_Fpass1.move(50, 115)

        self.OrderVal = QtGui.QSpinBox(self)
        self.OrderVal.setRange(1, 499)
        self.OrderVal.resize(80, 20)
        self.OrderVal.move(100, 120)
        
        self.Fpass1Val = QtGui.QSpinBox(self)
        self.Fpass1Val.setRange(1, 499)
        self.Fpass1Val.resize(80, 20)
        self.Fpass1Val.move(100, 120)
        
        self.Fpass1Val.setHidden(True)
        
        #parametr 2 - FIR: -, IIR: Fpass2

        self.Fpass2 = pg.QtGui.QLabel('Fpass2', self)
        self.Fpass2.move(290, 115)
        
        self.Fpass2.setHidden(True)

        self.Fpass2Val = QtGui.QSpinBox(self)
        self.Fpass2Val.setRange(1, 22050)
        self.Fpass2Val.resize(80, 20)
        self.Fpass2Val.move(350, 120)
        
        self.Fpass2Val.setHidden(True)
        
        #parametr 3 - FIR: Fc1, IIR: Fstop1

        self.Fc1_Fstop1 = pg.QtGui.QLabel('Fc1:', self)
        self.Fc1_Fstop1.move(50, 155)

        self.Fc1Val = QtGui.QSpinBox(self)
        self.Fc1Val.setRange(1, 22050)
        self.Fc1Val.resize(80, 20)
        self.Fc1Val.move(100, 160)
        
        self.Fstop1Val = QtGui.QSpinBox(self)
        self.Fstop1Val.setRange(1, 22050)
        self.Fstop1Val.resize(80, 20)
        self.Fstop1Val.move(100, 160)
        
        self.Fstop1Val.setHidden(True)

        #parametr 4 - FIR: -, IIR: Fstop2
        
        self.Fstop2 = pg.QtGui.QLabel('Fstop2', self)
        self.Fstop2.move(290, 155)
        
        self.Fstop2.setHidden(True)

        self.Fstop2Val = QtGui.QSpinBox(self)
        self.Fstop2Val.setRange(1, 22050)
        self.Fstop2Val.resize(80, 20)
        self.Fstop2Val.move(350, 160)
        
        self.Fstop2Val.setHidden(True)

        #parametr 5, FIR: Fc2, IIR: Apass & parametr 6: FIR: -, IIR:Astop

        self.Fc2_Apass = pg.QtGui.QLabel('Fc2:', self)
        self.Fc2_Apass.move(50, 195)

        self.Fc2Val = QtGui.QSpinBox(self)
        self.Fc2Val.setRange(1, 22050)
        self.Fc2Val.resize(80, 20)
        self.Fc2Val.move(100, 200)
        
        self.ApassVal = QtGui.QSpinBox(self)
        self.ApassVal.setRange(1, 22050)
        self.ApassVal.resize(80, 20)
        self.ApassVal.move(100, 200)
        
        self.ApassVal.setHidden(True)

        self.Astop = pg.QtGui.QLabel('Astop', self)
        self.Astop.move(290, 195)
        
        self.Astop.setHidden(True)

        self.AstopVal = QtGui.QSpinBox(self)
        self.AstopVal.setRange(1, 150)
        self.AstopVal.resize(80, 20)
        self.AstopVal.move(350, 200)
        
        self.AstopVal.setHidden(True)

        #Init UIcomponents disable
        
        #bandpass/bandstop spinbox disable
        self.Fc2Val.setEnabled(False)

        self.structureBox.setEnabled(False)        
        
        #IIR windows comboboox disable 
        self.typeIIRBox.setEnabled(False)


        self.show()

    def Characteristic_Change(self):        #LP, HP, BP, BS
        if(self.responseTypeBox.currentText() == "Bandpass" or
           self.responseTypeBox.currentText() == "Bandstop"):
            
            if(self.responseBox.currentText() == "FIR"):
                self.Fc2Val.setEnabled(True)
            elif(self.responseBox.currentText() == "IIR"):
                self.Fpass2Val.setEnabled(True)
                self.Fstop2Val.setEnabled(True)

        else:
            if(self.responseBox.currentText() == "FIR"):
                self.Fc2Val.setEnabled(False)

            elif(self.responseBox.currentText() == "IIR"):
                self.Fpass2Val.setEnabled(False)
                self.Fstop2Val.setEnabled(False)

    def Response_Change(self):              # FIR/IIR
        if (self.responseBox.currentText() == "FIR"):
        
            self.typeFIRBox.setEnabled(True)
            self.typeIIRBox.setEnabled(False)
            self.structureBox.setEnabled(False)
            
            self.OrderVal.setHidden(False)
            self.Fc1Val.setHidden(False)
            self.Fc2Val.setHidden(False)
            
            self.Fstop2.setHidden(True)
            self.Fpass2.setHidden(True)
            
            self.Fpass1Val.setHidden(True)
            self.Fstop1Val.setHidden(True)
            
            self.Fpass2Val.setHidden(True)
            self.Fstop2Val.setHidden(True)
            
            self.ApassVal.setHidden(True)
            
            self.Astop.setHidden(True)
            self.AstopVal.setHidden(True)
            
            self.Order_Fpass1.setText("Order")
            self.Fc1_Fstop1.setText("Fc1")
            self.Fc2_Apass.setText("Fc2")
            
            self.OrderVal.setMaximum(499)
            self.Fc2Val.setMaximum(fs/2)
        else:
            self.typeFIRBox.setEnabled(False)
            self.typeIIRBox.setEnabled(True)
            self.structureBox.setEnabled(True)
            
            self.OrderVal.setHidden(True)
            self.Fc1Val.setHidden(True)
            self.Fc2Val.setHidden(True)
            
            self.Fstop2.setHidden(False)
            self.Fpass2.setHidden(False)
            
            self.Fpass1Val.setHidden(False)
            self.Fstop1Val.setHidden(False)
            
            self.Fpass2Val.setHidden(False)
            self.Fstop2Val.setHidden(False)
            
            self.ApassVal.setHidden(False)
            
            self.Astop.setHidden(False)
            self.AstopVal.setHidden(False)
            
            self.Order_Fpass1.setText("Fpass1")
            self.Fc1_Fstop1.setText("Fstop1")
            self.Fc2_Apass.setText("Apass")
            
            self.Fpass1Val.setMaximum(fs/2)
            self.AstopVal.setMaximum(150)
            
        self.Characteristic_Change()   

    def Compute(self):
        if(self.responseBox.currentText() == "FIR"):
            tempName = self.typeFIRBox.currentText()
            windowName = tempName[0].lower() + tempName[1:]
            if (self.responseTypeBox.currentText() == "Lowpass"):
                self.b = signal.firwin(numtaps = self.OrderVal.value() + 1, cutoff = self.Fc1Val.value(), fs = fs, window = windowName)
   
            elif (self.responseTypeBox.currentText() == "Highpass"):
                self.b = signal.firwin(numtaps = self.OrderVal.value() + 1, cutoff = self.Fc1Val.value(), fs = fs, window = windowName, pass_zero=False)

            elif (self.responseTypeBox.currentText() == "Bandpass"):
                self.b = signal.firwin(numtaps = self.OrderVal.value() + 1, cutoff = [self.Fc1Val.value(), self.Fc2Val.value()], fs = fs, window = windowName, pass_zero=False)

            elif (self.responseTypeBox.currentText() == "Bandstop"):
                self.b = signal.firwin(numtaps = self.OrderVal.value() + 1, cutoff = [self.Fc1Val.value(), self.Fc2Val.value()], fs = fs, window = windowName)                                                                                              
                
            self.a = 1
            self.ftype = 0

            self.label3 = "direct"
            if(self.coding.currentText() == "floating point"):
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
            
            if(self.typeIIRBox.currentText() == "Butterworth"):
                iirTypeName = "butter"
            elif(self.typeIIRBox.currentText() == "Chebyshev type 1"):
                iirTypeName = "cheby1"
            elif(self.typeIIRBox.currentText() == "Chebyshev type 2"):
                iirTypeName = "cheby2"
            else:
                iirTypeName = "ellip"
            
            if (self.responseTypeBox.currentText() == "Lowpass" or self.responseTypeBox.currentText() == "Highpass"):
                
                if(self.responseTypeBox.currentText() == "Lowpass" and self.Fpass1Val.value() > self.Fstop1Val.value() or
                   self.responseTypeBox.currentText() == "Highpass" and self.Fpass1Val.value() < self.Fstop1Val.value()):
                    print("złe wartości") #tymczasowe rozwiazanie
                    self.b = 0
                    self.a = 0
                else:
                    self.b, self.a = signal.iirdesign(wp = self.Fpass1Val.value(), ws = self.Fstop1Val.value(), 
                                                      gpass = self.ApassVal.value(), gstop = self.AstopVal.value(), fs = fs, ftype = iirTypeName)
                    
            elif (self.responseTypeBox.currentText() == "Bandpass" or self.responseTypeBox.currentText() == "Bandstop"):
            
                x = {'fp1' : self.Fpass1Val.value(), 'fs1' : self.Fstop1Val.value(), 'fp2' : self.Fpass2Val.value(), 'fs2' : self.Fstop2Val.value()}
                
                sort_x = {k: v for k, v in sorted(x.items(), key = lambda item: item[1])}
                
                if((self.responseTypeBox.currentText() == "Bandpass" and (list(sort_x.keys()).index("fs1") != 0 or list(sort_x.keys()).index("fp1") != 1 or
                    list(sort_x.keys()).index("fp2") != 2 or list(sort_x.keys()).index("fs2") != 3)) or (self.responseTypeBox.currentText() == "Bandstop" and
                    (list(sort_x.keys()).index("fp1") != 0 or list(sort_x.keys()).index("fs1") != 1 or list(sort_x.keys()).index("fs2") != 2 or 
                    list(sort_x.keys()).index("fp2") != 3))):
                   
                    print("złe wartości") #tymczasowe rozwiazanie
                    self.b = 0
                    self.a = 0
                else:
                    self.b, self.a = signal.iirdesign(wp = [self.Fpass1Val.value(), self.Fpass2Val.value()], ws = [self.Fstop1Val.value(), self.Fstop2Val.value()],
                                                      gpass = self.ApassVal.value(), gstop = self.AstopVal.value(), fs = fs, ftype = iirTypeName)
            
            self.ftype = 1
            
            if(self.coding.currentText() != "floating point"):                
                self.q = int(self.coding.currentText()[1:3])
            
            if(self.structureBox.currentText() == "direct"):
                if(self.coding.currentText() == "floating point"):    
                    self.coefA = self.a
                    self.coefB = self.b
                    self.fstructure = 0
                else:
                    self.coefB = self.QuantizationLattDirect(self.b, self.q)
                    self.coefA = self.QuantizationLattDirect(self.a, self.q)
                    self.fstructure = 2
                
            elif(self.structureBox.currentText() == "lattice"):
                Kr, V = self.LatticeCoefsCompute(self.a, self.b)
                if(self.coding.currentText() == "floating point"): 
                    self.coefA = Kr
                    self.coefB = V
                    self.fstructure = 1
                else:
                    self.coefB = np.divide(self.QuantizationLattDirect(V, self.q), 2**self.q)
                    self.coefA = np.divide(self.QuantizationLattDirect(Kr, self.q), 2**self.q)
                    self.fstructure = 1

            elif(self.structureBox.currentText() == "biquad"):
                tempSOS = signal.tf2sos(self.b, self.a, 'nearest')
                self.sectionsNo = tempSOS.shape[0]
                if(self.coding.currentText() == "floating point"):
                    self.sos = tempSOS
                    self.fstructure = 4
                else:
                    self.sos = np.divide(self.QuantizationBiquad(tempSOS, self.q), 2**self.q)
                    self.fstructure = 4

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
        figure(1).clf()
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
        
        if(self.coding.currentText() == "floating point"):
            if(self.structureBox.currentText() == "biquad"):
                self.ImportBiquadToCFloatingP(self.sos, self.sectionsNo) 
            else:
                self.ImportToCFloatingP(self.ftype, self.coefB, self.coefA)
        else:
            wsk_q = ct.c_int.in_dll(lib, 'q')
            wsk_q.value = self.q
            
            if(self.structureBox.currentText() == "biquad"):
                self.ImportBiquadToCFloatingP(self.sos, self.sectionsNo)
                        
            elif(self.structureBox.currentText() == "lattice"):
                self.ImportToCFloatingP(self.ftype, self.coefB, self.coefA)
            else:
                self.ImportToCFixedP(self.ftype, self.coefB, self.coefA)
    
    def ImportToCFixedP(self, ftype, coefB, coefA):
        type_for_BnI = (ct.c_int * len(coefB))
        BnI = type_for_BnI()
                
        for i in range(len(coefB)):
            BnI[i] = int(self.coefB[i])
                
        if(ftype == 0):
            lib.importNumeratorI.argtypes = ct.POINTER(ct.c_int),ct.c_size_t
            lib.importNumeratorI(BnI, len(BnI))    #wywolanie funcji z srodowiska C
                    
        else:
            wsk_lenA = ct.c_int.in_dll(lib, 'lenA')
            wsk_lenA.value = len(coefA)
            
            type_for_BdnI = (ct.c_int * len(coefA))
            BdnI = type_for_BdnI()
            
            for i in range(len(self.coefA)):
                BdnI[i] = int(coefA[i])
            
            lib.importNumeratorI.argtypes = ct.POINTER(ct.c_int),ct.c_size_t
            lib.importNumeratorI(BnI, len(BnI))
            
            lib.importDenominatorI.argtypes = ct.POINTER(ct.c_int), ct.c_size_t
            lib.importDenominatorI(BdnI, len(BdnI))

    def ImportToCFloatingP(self, ftype, coefB, coefA):
        wsk_lenB = ct.c_int.in_dll(lib, 'lenB')
        wsk_lenB.value = len(coefB)
        
        type_for_Bn = (ct.c_float * len(coefB))
        Bn = type_for_Bn()
        
        for i in range(len(coefB)):
            Bn[i] = coefB[i]
        
        if(ftype == 0):
            lib.importNumerator.argtypes = ct.POINTER(ct.c_float),ct.c_size_t
            lib.importNumerator(Bn, len(Bn))    #wywolanie funcji z srodowiska C
            
        else:
            wsk_lenA = ct.c_int.in_dll(lib, 'lenA')
            wsk_lenA.value = len(coefA)
            
            type_for_Bdn = (ct.c_float * len(coefA))
            Bdn = type_for_Bdn()
            
            for i in range(len(coefA)):
                Bdn[i] = coefA[i]
            
            lib.importNumerator.argtypes = ct.POINTER(ct.c_float),ct.c_size_t
            lib.importNumerator(Bn, len(Bn))
            
            lib.importDenominator.argtypes = ct.POINTER(ct.c_float), ct.c_size_t
            lib.importDenominator(Bdn, len(Bdn))
    
    def ImportBiquadToCFloatingP(self, coefSOS, sectionsNo):
        wsk_sectionNo = ct.c_int.in_dll(lib, 'sectionsNo')
        wsk_sectionNo.value = sectionsNo
        
        type_for_Bn = ct.c_float*6
        Bn = type_for_Bn()
        lib.importBiquad.argtype = ct.POINTER(ct.c_float), ct.c_int
        type_for_ci = ct.c_int
        ci = type_for_ci()
        
        for i in range(sectionsNo):
            ci = i  
            for j in range(len(coefSOS[0])):
                Bn[j] = coefSOS[i][j]
            lib.importBiquad(Bn, ci)           
    
    def ImportBiquadToCFixedP(self, coefSOS, sectionsNo):
        wsk_sectionNo = ct.c_int.in_dll(lib, 'sectionsNo')
        wsk_sectionNo.value = sectionsNo
        
        type_for_Bn = ct.c_int*6
        Bn = type_for_Bn()
        lib.importBiquadI.argtype = ct.POINTER(ct.c_int), ct.c_int
        type_for_ci = ct.c_int
        ci = type_for_ci()
        
        for i in range(sectionsNo):
            ci = i  
            for j in range(len(coefSOS[0])):
                Bn[j] = int(coefSOS[i][j])
            lib.importBiquadI(Bn, ci)
            
    def QuantizationLattDirect(self, coefs, q):
        temp = [0] * len(coefs)
        for i in range(len(coefs)):
            temp[i] = round(coefs[i] * (2**q))
        return temp
        
    def QuantizationBiquad(self, coefs, q):
        sosI = np.zeros((coefs.shape[0],6))
        for i in range(coefs.shape[0]):
            for j in range(len(coefs[0])):
                sosI[i][j] = round(coefs[i][j] * (2**q))
                            
        return sosI
          
    def LatticeCoefsCompute(self, a, b):
        R, U, tempKr, e = levi.rlevinson(a,1)
        Kr = -tempKr
        V = np.zeros(len(a))
    
        for i in range(len(a)-1,-1, -1):
            subterm = U@V
            V[i] = b[i] - subterm[i]
            
        return Kr, V