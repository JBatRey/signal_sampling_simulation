from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QPixmap
import PyQt5.QtCore

import qdarktheme
import sys
import mainwin
import matplotlib.pyplot as plt
import numpy as np



from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import random

from backend import *

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

FSPAN = 100000


class myApp(QtWidgets.QMainWindow, mainwin.Ui_MainWindow):
    def __init__(self, parent=None):
        super(myApp, self).__init__(parent)
        self.setupUi(self)

        self.resize(1000,1000)

        self.checkBox_nodo4.setChecked(True)
        self.checkBox_nodo4.stateChanged.connect(self.newvalues)
        self.checkBox_nodo3.setChecked(True)
        self.checkBox_nodo3.stateChanged.connect(self.newvalues)
        self.checkBox_nodo2.setChecked(True)
        self.checkBox_nodo2.stateChanged.connect(self.newvalues)
        self.checkBox_nodo1.setChecked(True)
        self.checkBox_nodo1.stateChanged.connect(self.newvalues)

        
        self.label = QtWidgets.QLabel(self)
        self.pixmap = QPixmap('./imgs/sistema.jpg')
        self.pixmap = self.pixmap.scaled(1000, 1000, PyQt5.QtCore.Qt.KeepAspectRatio, PyQt5.QtCore.Qt.SmoothTransformation)
        self.label.setPixmap(self.pixmap)
        
        self.label.resize(self.pixmap.width(),
                          self.pixmap.height())

        self.layout_a.addWidget(self.label)

        self.slider_fsampling.setRange(1, 40000)
        self.slider_fsampling.setSingleStep(1)
        self.slider_fsampling.setValue(20000)
        self.label_fsampleo.setText('f sampleo: 20000Hz')
        self.slider_fsampling.valueChanged.connect(self.print_fsampling_value)
        self.slider_fsampling.sliderReleased.connect(self.newvalues)

        self.slider_fsignal.setRange(1,20000)
        self.slider_fsignal.setSingleStep(1)
        self.slider_fsignal.setValue(7000)
        self.label_fsignal.setText('f señal: 7000Hz')
        self.slider_fsignal.valueChanged.connect(self.print_fsignal_value)
        self.slider_fsignal.sliderReleased.connect(self.newvalues)
        
        self.filter_slider.setRange(1, 40000)
        self.filter_slider.setSingleStep(1)
        self.filter_slider.setValue(10000)
        self.label_filter.setText('f filtro: 10000Hz')
        self.filter_slider.valueChanged.connect(self.print_filter_value)
        self.filter_slider.sliderReleased.connect(self.newvalues)

        self.dc_slider.setRange(1,50)
        self.dc_slider.setSingleStep(1)
        self.dc_slider.setValue(25)
        self.label_dc.setText('DC llave analógica: 25%')
        self.dc_slider.valueChanged.connect(self.print_dc_value)
        self.dc_slider.sliderReleased.connect(self.newvalues)

        self.combo_sigtype.currentIndexChanged.connect(self.newvalues)

        self.radio_freq.setChecked(True)


        self.radio_time.clicked.connect(self.newvalues)
        self.radio_freq.clicked.connect(self.newvalues)


        # a figure instance to plot on
        self.figure = plt.figure()
  
        # this is the Canvas Widget that
        # displays the 'figure'it takes the
        # 'figure' instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
  
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)  
                  
        # adding tool bar to the layout
        self.layout_graph.addWidget(self.toolbar)
          
        # adding canvas to the layout
        self.layout_graph.addWidget(self.canvas)

        self.newvalues()
          
          

    def print_fsampling_value(self):

        self.label_fsampleo.setText('f sampleo: ' + str(self.slider_fsampling.value()))

    def print_fsignal_value(self):

        self.label_fsignal.setText('f señal: ' + str(self.slider_fsignal.value()))
    
    def print_filter_value(self):

        self.label_filter.setText('f filtro: ' + str(self.filter_slider.value()))

    def print_dc_value(self):

        self.label_dc.setText('DC llave analógica: ' + str(self.dc_slider.value()) + '%')

    def newvalues(self):

        fSampling = self.slider_fsampling.value()
        fSignal = self.slider_fsignal.value()
        fFilter = self.filter_slider.value()

        showTimeDomain = self.radio_time.isChecked()
        signalType = self.combo_sigtype.currentText()

        if signalType == 'Sine':
            xti, yti = create_sine(fSignal)
            xfi, yfi = create_sine_fourier(fSignal, FSPAN) 
        elif signalType == 'Cosine':
            xti, yti = create_cosine(fSignal)
            xfi, yfi = create_cosine_fourier(fSignal, FSPAN) 
        if signalType == 'Sawtooth':
            xti, yti = create_sawtooth(fSignal)
            xfi, yfi = create_sawtooth_fourier(fSignal, FSPAN) 
        if signalType == 'Square':
            xti, yti = create_square(fSignal)
            xfi, yfi = create_square_fourier(fSignal, FSPAN) 

                
        xto, yto = xti, yti
        xfo, yfo = xfi, yfi 


        if self.checkBox_nodo1.isChecked():
            xfo, yfo = lp_filter_cheby1(xfo, yfo, fFilter, 2*fFilter, 1, 40)
            xto, yto = signal_from_fourier(yfo, 1/fSignal, 2)




        if self.checkBox_nodo2.isChecked():
            cSH = True
        else:
            cSH = False

        if self.checkBox_nodo3.isChecked():
            cA = True
        else:
            cA = False
        
        tau = (self.dc_slider.value()/100)/fSampling

        if (cSH and cA):
            xfo, yfo = create_instantaneous_fourier(xfo, yfo, tau, fSampling)
            xto, yto = graph_instant(xto, yto, fSampling, self.dc_slider.value())
        elif (cSH and not cA):
            xfo, yfo = create_instantaneous_fourier(xfo, yfo, 1/fSampling, fSampling)
            xto, yto = graph_instant(xto, yto, fSampling, 100)
        elif (not cSH and cA):
            xfo, yfo = create_natural_fourier(xfo, yfo, tau, fSampling)
            xto, yto = graph_natural(xto, yto, fSampling, self.dc_slider.value())
        else:
            pass

        if self.checkBox_nodo4.isChecked():
            xfo, yfo = lp_filter_cheby1(xfo, yfo, fFilter, 2*fFilter, 1, 40)
            xto, yto = signal_from_fourier(yfo, 1/fSignal, 2)

        # clearing old figure
        self.figure.clear()
            
        # create an axis
        ax = self.figure.add_subplot(111)
  
        # plot data
        #ax.plot(xto, np.zeros(xto.size), 'k')
        if showTimeDomain != True:
            ax.plot(xfo,abs(yfo))
            ax.set_xlim(-40000, 40000)
        else: 
            ax.plot(xto, yto)
            ax.set_ylim(-1,1)
            
        ax.grid()
  
        # refresh canvas
        self.canvas.draw()
    

def main():
    app = QApplication(sys.argv)
    qdarktheme.setup_theme()

    form = myApp()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()