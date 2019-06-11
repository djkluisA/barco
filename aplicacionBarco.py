#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys    
import usv    #fichero de Juan Jimenez
import boom    #fichero de Juan Jimenez
import time    
import corrientes   #fichero de Juan Jimenez
import bezier_cvr as bz   #fichero de Juan Jimenez
import pickle
import copy_reg
import types
import os
from numpy.linalg.linalg import norm
from copy import copy
from PyQt4.Qt import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from matplotlib.path import Path
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as patches
from os import mkdir, stat, path



#tabla de explicaciones de todos los atributos que aparecen para cada clase.
tabla_explicacion={
"identificador":u"Identificador único(número entero)",
"tipo":u"Define el  sistema de propulsión",
"tiempo":u"Tiempo inicial de simulación",
"alfab":u"Aceleración angular (rad/s^2)",
"wb":u"Velocidad angular (rad/s)",
"theta":u"rumbo actual (rad)",
"Fm":u"fuerza de propulsión (N)",
"setl":u"Valor de consigna del motor",
"setw":u"Valor de consigna para el rumbo",
"interrv":u"integral del error de velocidad",
"interrumbo":u"Integral del error de rumbo",
"thewj":u"Ángulo del timón o del waterjet",
"thewjmax":u"Ángulo máximo de giro del timón o waterjet ",
"Ac":u"Distancia del eje del timón al centro de masa",
"M":u"Torque aplicado al barco por la fuerza de propulsión",
"mb":u"Masa del barco (Kg)",
"mul2":u"Coeficiente de resistencia al avance (Término cuadrático)",
"mul":u"Coeficiente de ressitencia al avance (Término lineal)",
"mut":u"Coeficiente de resistencia al desplazamiento lateral (Término lineal)",
"mut2":u"Coeficiente de resitencia al desplazamiento lateral (Término cuadrático)",
"mua":u"Coeficiente de resistencia al giro (Término lineal)",
"mua2":u"Coeficiente de resistencia al giro (Término cuadrático)",
"ls":u"Eslora ",
"pmax":u"Potencia máxima del motor",
"Tsf":u"Constante de tiempo del motor del barco",
"Tsw":u"Constante de tiempo del motor del sistema de guiado",
"cfr":u"factor de forma para los barcos guiados por timón",
"link":u"Enlace radio entre barcos",
"krpid[0]":u"Constante proporcional control PID rumbo",
"krpid[1]":u"Constante integral control PID rumbo",
"krpid[2]":u"Constante derivativa control PID rumbo",
"Ib":u"Momento de inercia",
"kvpid[0]":u"Constante proporcional control PID velocidad",
"kvpid[1]":u"Constante integral control PID velocidad",
"kvpid[2]":u"Constante derivativa control PID velocidad",
"krd[0]":u"Constante proporcional control PID distancia entre USVs",
"krd[1]":u"Constante derivativa control PID distancia entre USVs",
"krd[2]":u"Constante integral control PID distancia entre USVs",
"krl[0]":u"Constante proporcional control PID alineación entre USVs",
"krl[1]":u"Constante derivativa control PID alineación entre USVs",
"krl[2]":u"Constante integral control PID alineación entre USVs",
"ab[0]":u"aceleración x",
"ab[1]":u"aceleración y",
"pb[0]":u"posición x",
"pb[1]":u"posición y",
"vb[0]":u"velocidad x",
"vb[1]":u"velocidad y",
"Fb[0]":u"fuerza externa ejercida x",
"Fb[1]":u"fuerza de externa ejercida y",
"mA[0]":u"Masas anadidas x",
"mA[1]":u"Masas anadidas y",
"interdist[0]":u"Integral del error de distancia entre USVs x",
"interdist[1]":u"Integral del error de distnacia entre USVs y",
"cable[0]":u"longitud Cable remolque izquierdo",
"cable[1]":u"longitud Cable remolque derecho",
"cable[2]":u"Constante de elasticidad",
"Anch[0]":u"Extremo izquierdo de la barrera anclado a tierra 0/1",
"Anch[1]":u"Extremo deredcho de la barrera anclado a tierra 0/1",
"Fd[0]":u"Fuerza aplicada al extremo derecho de la barrera x",
"Fd[1]":u"Fuerza aplicada al extremo derecho de la barrera y",
"Fi[0]":u"Fuerza aplicada al extremo izquierdo de la barrera x",
"Fi[1]":u"Fuerza aplicada al extremo izquierdo de la barrera y",
"esl":u"Número de elementos (Eslabones)",
"objeto":u"Identificador único de la barrera (nombre)",
"s":u"superficie de un elemento",
"q":u"Espesor de un elemento de la barrera",
"s2":u"Coeficiente de resistencia al desplazamiento transversal de un elemento (Término cuadrático)",
"q2":u"Coeficiente de resistencia al desplazamiento longitudonal de un elemento (Término cuadrático)",
"A":u"Coeficiente de resistencia al giro de un elemento",
"L":u"Semilongitud de un elemento",
"m":u"Masa de un elemento",
"mA":"Masa anadida a un elemento",
"step":u"Número de iteraciones",
"delta":u"Paso de integración (tiempo)",
"tam":u"Tamano del experimento (número total de iteraciones)",
"d_rumbo":u"consigna de rumbo",
"d_sway":u"Consigna de distancia lateral entre USVs",
"d_surge":u"Consigna de alineamiento entre USVs"
}












#tabla usada para mostrar los atributos para la clase barco cadena y en la interfaz simulacion.

class MyTable(QTableWidget):
    
    def __init__(self, data,c,d,objeto_clase):
        QTableWidget.__init__(self,c,d)
        
        self.objeto_clase=objeto_clase
        self.data = data
        self.setRowCount(c)
        self.setColumnCount(d)
        self.setmydata()
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.cellChanged.connect(self.tableItemChanged)
        self.cellClicked.connect(self.tableItemExplic)
        self.verticalHeader().setStretchLastSection(True)
        self.horizontalHeader().setStretchLastSection(True)
		
	#función que abre una ventana de explicación para los atributos	
    def tableItemExplic(self,x,y):
        t=self.item(x,y)
        titulo=self.item(x,0)
        if y==2:
			self.contener = MyContenido(t.text())
			self.contener.setGeometry(QRect(100, 100, 400, 50))
			self.contener.setWindowIcon(QIcon('images.ico'))
			self.contener.setWindowTitle(titulo.text())
			self.contener.show()
        
    #función que cambia el valor de los atributos al cambiar cada una de las celdas de la tabla creada 
    #para cada interfaz: interfaz barco, interfaz cadena, interfaz simular	
    def tableItemChanged(self,x,y):
        
        t=self.item(x,y-1)
        
        
        for key in self.data.keys():
              if key==t.text():
                if str(t.text()).find("rumbo")!=-1:
                   self.data[key]=float(self.item(x,y).text())*np.pi/180
                else:
                   self.data[key]=self.item(x,y).text()
                if self.objeto_clase!=None:
                   
                   if isinstance(self.objeto_clase,boom.cadena):
                    if str(t.text()).find("[0]")!=-1  or str(t.text()).find("[1]")!=-1 or str(t.text()).find("[2]")!=-1:
                      self.objeto_clase.setAtributo(str(t.text()),self.item(x,y).text())
                    else:
                      if str(t.text()).find("objeto")==-1:
                        if str(t.text()).find("esl")!=-1:
                          setattr(self.objeto_clase,str(t.text()),int(self.item(x,y).text()))
                          self.objeto_clase.setEslabones()
                          
                        
                        else: 
                          setattr(self.objeto_clase,str(t.text()),float(self.item(x,y).text()))
                
                   if isinstance(self.objeto_clase,usv.barco ):
                    if str(t.text()).find("[0]")!=-1  or str(t.text()).find("[1]")!=-1 or str(t.text()).find("[2]")!=-1:
                        self.objeto_clase.setAtributo(str(t.text()),self.item(x,y).text())
                    else:
                        if str(t.text()).find("identificador")==-1:
                           if str(t.text()).find("theta")!=-1 or str(t.text()).find("thewjmax")!=-1:
                              setattr(self.objeto_clase,str(t.text()),float(self.item(x,y).text())*np.pi/180)
                           else:
                              if str(t.text()).find("tipo")!=-1:
                                setattr(self.objeto_clase,str(t.text()),str(self.item(x,y).text()))
                              else:
                                setattr(self.objeto_clase,str(t.text()),float(self.item(x,y).text()))
                  
                  
                
        
	#función que construye la tabla(las columnas, las filas en donde hace que la segunda columna sea #inmutable)
    def setmydata(self):
        
        
        for n, key in enumerate(self.data.keys()):
              newitem1 = QTableWidgetItem(str(key))
              newitem1.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
              newitem1.setFlags(Qt.ItemIsEnabled)
              newitem1.setBackground(Qt.lightGray)
              self.setItem(n,0,newitem1)
              newitem2 = QTableWidgetItem(str(self.data[key]))
              newitem2.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
              
              self.setItem(n, 1, newitem2)
              try:
                newitem3=QTableWidgetItem(tabla_explicacion[key])
                newitem3.setFlags(Qt.ItemIsEnabled)
                newitem3.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
                newitem3.setBackground(Qt.lightGray)
                self.setItem(n, 2, newitem3)
              except:
                print "error al cargar la tabla"
			  
              
        
        afont = QFont()
        afont.setFamily("Arial Black")
        afont.setPointSize(8)
        item1 = QTableWidgetItem(u"Parámetro")
        item1.setFont(afont)
        item1.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        self.setHorizontalHeaderItem(0,item1)
        

        item2 = QTableWidgetItem(u'Valor')
        item2.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        item2.setFont(afont)
        self.setHorizontalHeaderItem(1,item2)

        item3 = QTableWidgetItem(u'Explicación')
        item3.setTextAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        item3.setFont(afont)
        self.setHorizontalHeaderItem(2,item3)
        self.verticalHeader().setVisible(False)
 

#interfaz del contenido que muestra la ventana de explicación  
class MyContenido(QWidget):
    def __init__(self,cadena):
     QWidget.__init__(self)
     with open("style.css") as f:
          self.setStyleSheet(f.read())
     self.mostrar = QLabel()
     self.mostrar.setText(cadena)
     layout=QHBoxLayout()
     layout.addWidget(self.mostrar)
     self.setLayout(layout)

#esta clase crea el check para escoger los distintos atributos a comunicar
class CheckableComboBox(QComboBox):
    def __init__(self,dato):
        QComboBox.__init__(self)
        self.view().pressed.connect(self.handleItemPressed)
        self.setModel(QStandardItemModel(self))
        self.dato=dato
           
	#función que cambia los checks de estado
    def handleItemPressed(self, index):
        item = self.model().itemFromIndex(index)			
        if item.checkState() == Qt.Checked:
            item.setCheckState(Qt.Unchecked)
            t=item.row()
            self.dato.pop(int(t))
            
            
        else:
            item.setCheckState(Qt.Checked)
            t=item.row()
            r=self.itemText(int(t))
            self.dato[int(t)]=str(r)
            

#clase que crea una instancia de la clase check arriba indicada.
#Esta clase es la que se configura para comunicar los distintos atributos entre las instancias de los barcos.        
class MyRadioLink(QWidget):
    def __init__(self,barcoizquierdo,barcoderecho):
     QWidget.__init__(self)
     with open("stylecombo.css") as f:
          self.setStyleSheet(f.read())
     self.b = QWidget(self)
     self.barcoi=barcoizquierdo
     self.barcod=barcoderecho
     self.dato=dict()
     self.Combobox = CheckableComboBox(self.dato)
     self.btn = QPushButton("Guardar RadioLink")
     self.connect(self.btn, SIGNAL("clicked()"), self.guardar)
     self.Combobox.addItem("Darle Check al Atributo a Comunicar")
     tablerolink = {1:'ide',2:'tipo',3:'tiempo', 4:'alfab', 5:'wb',6:'theta',7:'Fm',8:'setl',9:'setw',
                   10:'interrv',11:'interrumbo',12:'thewj',13:'thewjmax',14:'Ac',15:'M',16:'mb',17:'mul2',18:'mul',
                   19:'mut',20:'mut2',21:'Ib',22:'mua',23:'mua2',24:'ls',25:'pmax',26:'Tsf',27:'Tsw',28:'cfr',29:'link',30:'krpid',31:'krd',32:'krl',33:'kvpid',34:'ab',35:'pb',36:'vb',37:'Fb',38:'mA',39:'interdist'}
        
     for i in range(len(tablerolink)):
        self.Combobox.addItem("" + str(tablerolink[i+1]))
        item = self.Combobox.model().item(i+1,0)
        item.setCheckState(Qt.Unchecked)

     self.mostrarError = QLabel()
     self.mostrarError.setOpenExternalLinks(False)
     self.layout = QVBoxLayout(self)
     self.layout.addWidget(self.Combobox)
     self.layout.addWidget(self.btn)
     self.layout.addWidget(self.mostrarError)
     self.setLayout(self.layout)
    #este es el botón que se le dá para guardar la configuración elegida
    def guardar(self):
     try:
      
       li=[]
       ld=[]
       for atributos in self.dato.values():
	    
         li.append(self.barcoi.getAtributo(atributos))
         ld.append(self.barcod.getAtributo(atributos))
       self.barcoi.rlink = [li,[]]
       self.barcod.rlink = [ld,[]]
       app.main.combina.btn.setStyleSheet("color:white;")
       app.main.combina.btn.setEnabled(True)
   
       self.close()
     except:
        self.mostrarError.setText("Falla, no hay barcos ni cadenas, cerrar ventana")
        self.mostrarError.setOpenExternalLinks(True)
        


	 

#clase de la interfaz barco.
class MyBarco(QWidget):
    
    def __init__(self,listabarcos,barco,tab,krpid_tab,krd_tab,krl_tab,kvpid_tab,
		                          ab_tab,pb_tab,vb_tab,Fb_tab,mA_tab,interdist_tab):
        QWidget.__init__(self)
        self.b = QWidget(self)
        with open("barcostyle.css") as f:
          self.setStyleSheet(f.read())
        tablero1=copy(tab)
        self.barco=barco
        self.existe_barco=False
        self.listabarcos=listabarcos
        tablero_krpid=copy(krpid_tab)
        tablero_krd=copy(krd_tab)
        tablero_krl=copy(krl_tab)
        tablero_kvpid=copy(kvpid_tab)
        tablero_ab=copy(ab_tab)
        tablero_pb=copy(pb_tab)
        tablero_vb=copy(vb_tab)
        tablero_Fb=copy(Fb_tab)
        tablero_mA=copy(mA_tab)
        tablero_interdist=copy(interdist_tab)
		
        self.Etiqueta = QLabel(self.b)
        
        self.Etiqueta.setOpenExternalLinks(True)
        self.Etiqueta.setText("Estructuras Variables Registros")
        self.EtiquetaE = QLabel(self.b)
        
        self.EtiquetaE.setOpenExternalLinks(True)
        self.EtiquetaE.setText("Variables")
        self.mostrar = QLabel(self.b)
        self.mostrar.setStyleSheet("font-size: 20px; color: red;")
        self.mostrar.setOpenExternalLinks(False)
        self.btn = QPushButton("Mostrar Dibujo Barco", self.b)
        self.connect(self.btn, SIGNAL("clicked()"), self.dibujo)
        self.btn.setVisible(False)
        
        self.btn2 = QPushButton("Guardar Barco", self.b)
        
        self.connect(self.btn2, SIGNAL("clicked()"), self.guardar)
        if self.barco!=None: 
           self.btn.setVisible(True)
           
        self.tabla1=MyTable(tablero1,29,3,barco)
        self.tabla2=MyTable(tablero_krpid,3,3,barco)
        self.tabla3=MyTable(tablero_krd,3,3,barco)
        self.tabla4=MyTable(tablero_krl,3,3,barco)
        self.tabla5=MyTable(tablero_kvpid,3,3,barco)
        self.tabla6=MyTable(tablero_ab,2,3,barco)
        self.tabla7=MyTable(tablero_pb,2,3,barco)
        self.tabla8=MyTable(tablero_vb,2,3,barco)
        self.tabla9=MyTable(tablero_Fb,2,3,barco)
        self.tabla10=MyTable(tablero_mA,2,3,barco)
        self.tabla11=MyTable(tablero_interdist,2,3,barco)
        layout11 = QVBoxLayout()
        layout16 = QVBoxLayout()
        layout20 = QVBoxLayout()
        layout13=QHBoxLayout()
        layout15=QHBoxLayout()
        layout14=QHBoxLayout()
        layout17=QHBoxLayout()
        layout18=QHBoxLayout()
        layout19=QHBoxLayout()
        layout13.addWidget(self.tabla2)
        layout13.addWidget(self.tabla3)
        layout11.addWidget(self.Etiqueta)
        layout11.addLayout(layout13)
        layout14.addWidget(self.tabla4)
        layout14.addWidget(self.tabla5)
        layout11.addLayout(layout14)
        layout17.addWidget(self.tabla6)
        layout17.addWidget(self.tabla7)
        layout11.addLayout(layout17)
        layout18.addWidget(self.tabla8)
        layout18.addWidget(self.tabla9)
        layout11.addLayout(layout18)
        layout19.addWidget(self.tabla10)
        layout19.addWidget(self.tabla11)
        layout11.addLayout(layout19)
        layout20.addWidget(self.EtiquetaE)
        layout20.addWidget(self.tabla1)
        layout15.addLayout(layout20)
        layout15.addLayout(layout11)
        layout16.addLayout(layout15)
        layout16.addWidget(self.btn2)
        layout16.addWidget(self.btn)
        layout16.addWidget(self.mostrar)
        self.setLayout(layout16)
	#función que se usa en la clase editar cadena y barco para dibujar un barco.
    def dibujo(self):
        
        xi=self.barco.pb[0] + self.barco.ls
        yi=self.barco.pb[1] + self.barco.ls
        
        self.b.figure = pl.figure("Dibujo Barco")	
        ax = self.b.figure.add_subplot(111)
        ax.clear()
        ax.plot([-xi,-2,-1,0, 1, 2,xi], [-yi,-2,-1,0,1,2,yi],'+')
        self.canvas = FigureCanvas(self.b.figure)
        
        pl.plot(self.barco.pb[0],self.barco.pb[1],'+')
        patchi=self.barco.dibujar(self.barco.theta,self.barco.pb[0],self.barco.pb[1],'red')

        pl.gca().add_patch(patchi) 
        
        self.canvas.draw()
        pl.show()
        self.close()
	#función que se inicia al cerrar la ventana	
    def closeEvent(self, event):
        
        event.accept()
     #función que se usa para guardar una instancia barco.   
    def guardar(self):
        
        
        if self.barco==None:
          try:
           busca_barco=int(self.tabla1.data['identificador'])
           self.mostrar.setText("")
           if busca_barco >=0: 
            if busca_barco in self.listabarcos:
              self.existe_barco=True
            else:
              self.existe_barco=False
            if self.existe_barco==False:
            
              pb=np.array([float(self.tabla7.data['pb[0]']),float(self.tabla7.data['pb[1]'])])
              vb=np.array([float(self.tabla8.data['vb[0]']),float(self.tabla8.data['vb[1]'])])
              ab=np.array([float(self.tabla6.data['ab[0]']),float(self.tabla6.data['ab[1]'])])
              Fb=np.array([float(self.tabla9.data['Fb[0]']),float(self.tabla9.data['Fb[1]'])])
              mA=np.array([float(self.tabla10.data['mA[0]']),float(self.tabla10.data['mA[1]'])])
              interdist=np.array([float(self.tabla11.data['interdist[0]']),float(self.tabla11.data['interdist[1]'])])
              kvpid=np.array([float(self.tabla5.data['kvpid[0]']),float(self.tabla5.data['kvpid[1]']),float(self.tabla5.data['kvpid[2]'])])
              krpid=np.array([float(self.tabla2.data['krpid[0]']),float(self.tabla2.data['krpid[1]']),float(self.tabla2.data['krpid[2]'])])
              krd=np.array([float(self.tabla3.data['krd[0]']),float(self.tabla3.data['krd[1]']),float(self.tabla3.data['krd[2]'])])
              krl=np.array([float(self.tabla4.data['krl[0]']),float(self.tabla4.data['krl[1]']),float(self.tabla4.data['krl[2]'])])
              barco1=usv.barco(busca_barco,str(self.tabla1.data['tipo']),float(self.tabla1.data['tiempo']),pb,vb,ab,float(self.tabla1.data['alfab']),
		      float(self.tabla1.data['wb']),float(self.tabla1.data['theta'])*np.pi/180,float(self.tabla1.data['Fm']),Fb,float(self.tabla1.data['setl']),float(self.tabla1.data['setw']),
		      float(self.tabla1.data['interrv']),float(self.tabla1.data['interrumbo']),float(self.tabla1.data['thewj']),float(self.tabla1.data['thewjmax'])*np.pi/180,
		      float(self.tabla1.data['Ac']),float(self.tabla1.data['M']),float(self.tabla1.data['mb']),mA,float(self.tabla1.data['mul']),float(self.tabla1.data['mul2']),
		      float(self.tabla1.data['mut']),float(self.tabla1.data['mut2']),float(self.tabla1.data['Ib']),float(self.tabla1.data['mua']),float(self.tabla1.data['mua2']),
		      float(self.tabla1.data['ls']),float(self.tabla1.data['pmax']),float(self.tabla1.data['Tsf']),float(self.tabla1.data['Tsw'])
		      ,float(self.tabla1.data['cfr']),interdist,float(self.tabla1.data['link']),kvpid,krpid,krd,krl)
              app.main.lista_barcos[busca_barco]=barco1
              app.main.btnbyc.setEnabled(True)
              app.main.btnbyc.setStyleSheet("color:white;")
              if len(self.listabarcos)>= 2 and len(app.main.lista_cadenas) > 0:
                 app.main.btncmb.setEnabled(True)
                 app.main.btncmb.setStyleSheet("color:white;")
              self.close()
            
              
            else:
             self.mostrar.setText("no puede haber dos barcos con el mismo identificador ")
             self.mostrar.setOpenExternalLinks(True)
             existe_barco=False
           else:
            self.mostrar.setText("no puede haber identificadores menores que cero ")
            self.mostrar.setOpenExternalLinks(True)
          except ValueError:
           self.mostrar.setText("por favor respeta los tipos")
           self.mostrar.setOpenExternalLinks(True)
        else:
             
              self.close()

#clase que se usa en la interfaz cadena 		
class MyCadena(QWidget):
    
    def __init__(self,listacadenas,cadena,cad,cable_tabla,Anch_tabla,Fd_tabla,Fi_tabla):
        self.cadena=cadena
        self.existe_cadena=False
        self.listacadenas=listacadenas
        tablacad=copy(cad)
        tablacab=copy(cable_tabla)
        tablaanch=copy(Anch_tabla)
        tablafd=copy(Fd_tabla)
        tablafi=copy(Fi_tabla)
        QWidget.__init__(self)
        with open("cadenastyle.css") as f:
          self.setStyleSheet(f.read())
        self.b = QWidget(self)
        self.Etiqueta = QLabel(self.b)
      
        self.Etiqueta.setOpenExternalLinks(True)
        self.Etiqueta.setText("Estructuras Variables Registros")
        self.EtiquetaE = QLabel(self.b)
       
        self.EtiquetaE.setOpenExternalLinks(True)
        self.EtiquetaE.setText("Variables")
        self.mostrar = QLabel()
        self.mostrar.setStyleSheet("font-size: 20px; color: red;")
        self.mostrar.setOpenExternalLinks(False)
        
        
        self.btn2 = QPushButton("Mostrar Dibujo Cadena", self.b)
        self.connect(self.btn2, SIGNAL("clicked()"), self.dibujo_cadena)
        self.btn2.setVisible(False)
        if self.cadena!=None: 
           self.btn2.setVisible(True)
          
        self.btn7 = QPushButton("Guardar cadena", self.b)
       
        self.connect(self.btn7, SIGNAL("clicked()"), self.guardar_cadena)
        self.tabla1=MyTable(tablacad,10,3,cadena)
        self.tabla2=MyTable(tablacab,3,3,cadena)
        self.tabla3=MyTable(tablaanch,2,3,cadena)
        self.tabla4=MyTable(tablafd,2,3,cadena)
        self.tabla5=MyTable(tablafi,2,3,cadena)
        layout = QVBoxLayout()
        layout4 = QVBoxLayout()
        layout7 = QVBoxLayout()
        layout1=QHBoxLayout()
        layout3=QHBoxLayout()
        layout2=QHBoxLayout()
        layout5=QHBoxLayout()
        layout6=QHBoxLayout()
        layout.addWidget(self.Etiqueta)
        layout1.addWidget(self.tabla2)
        layout.addLayout(layout1)
        layout2.addWidget(self.tabla3)
        layout.addLayout(layout2)
        layout5.addWidget(self.tabla4)
        layout.addLayout(layout5)
        layout6.addWidget(self.tabla5)
        layout.addLayout(layout6)
        layout7.addWidget(self.EtiquetaE)
        layout7.addWidget(self.tabla1)
        layout3.addLayout(layout7)
        layout3.addLayout(layout)
        layout4.addLayout(layout3)
        layout4.addWidget(self.btn7)
        layout4.addWidget(self.btn2)
        layout4.addWidget(self.mostrar)
        self.setLayout(layout4)
	
	#función que se usa en la interfaz editar barco y cadena para dibujar una cadena
    def dibujo_cadena(self):
        xi=self.cadena.esl + 10
        yi=5
        self.b.figure = pl.figure("Dibujo Cadena")	
        ax = self.b.figure.add_subplot(111)
        ax.clear()
        ax.plot([-xi,-2,-1,0, 1, 2,xi], [-yi,-2,-1,0,1,2,yi],'+')
        self.canvas = FigureCanvas(self.b.figure)
        color='k'
        ax.plot(self.cadena.cms[0,:],self.cadena.cms[1,:],'o')
        
        barrasi,barrasd=self.cadena.dibujar()
        ax.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],color)
        
        self.canvas.draw()
        pl.show()
        self.close()
		
	#función que se inicia al cerrar la ventana
    def closeEvent(self, event):
                
        event.accept()
	#función que se usa para guardar la instancia cadena
    def guardar_cadena(self):
        
        
        busca_cadena=str(self.tabla1.data["objeto"])
        self.mostrar.setText("")
        if self.cadena==None:
         if busca_cadena.isdigit()==False:
          if busca_cadena in self.listacadenas:
              self.existe_cadena=True
          else:
              self.existe_cadena=False
          
          
          if self.existe_cadena==False:
            
              esl=int(self.tabla1.data["esl"])
              cable=np.array([float(self.tabla2.data["cable[0]"]),float(self.tabla2.data["cable[1]"]),float(self.tabla2.data["cable[2]"])])
              Anch=np.array([float(self.tabla3.data["Anch[0]"]),float(self.tabla3.data["Anch[1]"])])
              Fd=np.array([float(self.tabla4.data["Fd[0]"]),float(self.tabla4.data["Fd[1]"])])
              Fi=np.array([float(self.tabla5.data["Fi[0]"]),float(self.tabla5.data["Fi[1]"])])
              cad1=boom.cadena(esl,cable,Anch,busca_cadena,float(self.tabla1.data["s"]),float(self.tabla1.data["q"]),float(self.tabla1.data["s2"]),float(self.tabla1.data["q2"]),float(self.tabla1.data["A"]),
		                float(self.tabla1.data["L"]),float(self.tabla1.data["m"]),float(self.tabla1.data["mA"]),Fd,Fi)
						
        
              app.main.lista_cadenas[busca_cadena]=cad1
              
              app.main.btnbyc.setEnabled(True)
              app.main.btnbyc.setStyleSheet("color:white;")
              if len(app.main.lista_barcos) >= 2 and len(self.listacadenas) > 0:
                 app.main.btncmb.setEnabled(True)
                 app.main.btncmb.setStyleSheet("color:white;")
              
              self.close()
            
          else:
              self.mostrar.setText("No se puede repetir el nombre de la cadena")
              self.mostrar.setOpenExternalLinks(True)
              existe_cadena=False
              
         else:
              self.mostrar.setText("por favor es un nombre no un numero")
              self.mostrar.setOpenExternalLinks(True)
        else:
              
              
              self.close()	  
			  

#clase que se usa para abrir la interfaz edita barco y cadena.	
class MyEditaByC(QWidget):
    def __init__(self,listabarcos,listacadenas):
     
     QWidget.__init__(self)
     self.barco=None
     self.cadena=None
     self.listabarcos=listabarcos
     self.listacadenas=listacadenas
     self.b = QWidget(self)
     
     with open("stylebyc.css") as f:
          self.setStyleSheet(f.read())
     self.stacked_widget1 = QStackedWidget()
     self.stacked_widget = QStackedWidget()
     self.combobox2 = QComboBox(self.b)
     self.combobox3 = QComboBox(self.b)
     self.combobox = QComboBox(self.b)
     self.combobox.addItem("Elegir Barco Editable")
     self.combobox3.addItem("Elegir Cadena Editable")
     strlist = ['el identificador del barco es {}'.format(barcos) for barcos in listabarcos.keys()]
     strlist2 = ['el nombre de cadena es {}'.format(cadena1) for cadena1 in listacadenas.keys()]		
     self.combobox.addItems(strlist)
     self.combobox3.addItems(strlist2)
     self.combobox.activated.connect(self.editarbarco)
     self.combobox3.activated.connect(self.editarcadena)
     self.combobox2.addItems(['Barco', 'Cadena'])
     self.stacked_widget1.addWidget(self.combobox)
     self.stacked_widget1.addWidget(self.combobox3)
     self.combobox2.activated.connect(self.stacked_widget1.setCurrentIndex)
     self.combobox2.activated.connect(self.stacked_widget.setCurrentIndex)
     self.mostrar = QLabel()
     self.mostrar.setOpenExternalLinks(False)
     self.mostrar1 = QLabel()
     self.mostrar1.setOpenExternalLinks(False)
     self.statusBar = QStatusBar()
     self.layout = QVBoxLayout(self)
     self.layout.addWidget(self.combobox2)
     self.layout.addWidget(self.stacked_widget1)
     self.layout.addWidget(self.stacked_widget)
     self.layout.addWidget(self.mostrar)
     self.layout.addWidget(self.mostrar1)
     self.setLayout(self.layout)
	 
	#función que se usa para editar una cadena.
    def editarcadena(self,y):
        
        if y!=0:
         t=self.combobox3.currentText()
         t=t.replace('el nombre de cadena es','')
         cade1=str(t).strip()
         editaobjeto=self.listacadenas[cade1]
		 #####listas de las variables iniciales de las cadenas#########
         cadena = {'esl':editaobjeto.esl,'objeto':editaobjeto.objeto,'s':editaobjeto.s,'q':editaobjeto.q,
                  's2':editaobjeto.s2,'q2':editaobjeto.q2,'A':editaobjeto.A,'L':editaobjeto.L,'m':editaobjeto.m,'mA':editaobjeto.mA}
         cable_tabla={'cable[0]':editaobjeto.cbl[0], 'cable[1]':editaobjeto.cbl[1], 'cable[2]':editaobjeto.cbl[2]}
         Anch_tabla={'Anch[0]':editaobjeto.Anch[0],'Anch[1]':editaobjeto.Anch[1]}
         Fd_tabla={'Fd[0]':editaobjeto.Fd[0],'Fd[1]':editaobjeto.Fd[1]}
         Fi_tabla={'Fi[0]':editaobjeto.Fi[0],'Fi[1]':editaobjeto.Fi[1]}
		######termina la lista de las variables iniciales de la cadena#####
         self.cadena=editaobjeto
         self.editacadena = MyCadena(self.listacadenas,editaobjeto,cadena,cable_tabla,Anch_tabla,Fd_tabla,Fi_tabla)
         
         self.editacadena.setGeometry(QRect(50, 50, 1200, 400))
         self.editacadena.setWindowIcon(QIcon('images.ico'))
         self.editacadena.setWindowTitle('Editar Cadena')
         self.editacadena.show()
         
         
         
         
         self.combobox3.setCurrentIndex(0)
               
         
       
	#función que se usa para editar una cadena.		   
    def editarbarco(self,x):
        
        if x!=0:
         t=self.combobox.currentText()
         t=t.replace('el identificador del barco es','')
         ide=int(t)
         editaobjeto=self.listabarcos[ide]
		 #####listas de la vista de los barcos (variables iniciales)####
         tablero = {'identificador':editaobjeto.ide,'tipo':editaobjeto.tipo,'tiempo':editaobjeto.tiempo, 'alfab':editaobjeto.alfab, 'wb':editaobjeto.wb,'theta':editaobjeto.theta,'Fm':editaobjeto.Fm,'setl':editaobjeto.setl,'setw':editaobjeto.setw,
                   'interrv':editaobjeto.interrv,'interrumbo':editaobjeto.interrumbo,'thewj':editaobjeto.thewj,'thewjmax':editaobjeto.thewjmax,'Ac':editaobjeto.Ac,'M':editaobjeto.M,'mb':editaobjeto.mb,'mul2':editaobjeto.mul2,'mul':editaobjeto.mul,
                   'mut':editaobjeto.mut,'mut2':editaobjeto.mut,'Ib':editaobjeto.Ib,'mua':editaobjeto.mua,'mua2':editaobjeto.mua2,'ls':editaobjeto.ls,'pmax':editaobjeto.pmax,'Tsf':editaobjeto.Tsf,'Tsw':editaobjeto.Tsw,'cfr':editaobjeto.cfr,'link':editaobjeto.link}
         krpid_tabla={'krpid[0]':editaobjeto.krpid[0] ,'krpid[1]': editaobjeto.krpid[1]  ,'krpid[2]': editaobjeto.krpid[2]}	
         krd_tabla={'krd[0]':editaobjeto.krd[0] ,'krd[1]': editaobjeto.krd[1]  ,'krd[2]': editaobjeto.krd[2]}
         krl_tabla={'krl[0]':editaobjeto.krl[0] ,'krl[1]': editaobjeto.krl[1],'krl[2]': editaobjeto.krl[2]}
         kvpid_tabla={'kvpid[0]':editaobjeto.kvpid[0],'kvpid[1]': editaobjeto.kvpid[1],'kvpid[2]': editaobjeto.kvpid[2]}
         ab_tabla={'ab[0]':editaobjeto.ab[0],'ab[1]':editaobjeto.ab[1]}
         pb_tabla={'pb[0]':editaobjeto.pb[0],'pb[1]':editaobjeto.pb[1]}
         vb_tabla={'vb[0]':editaobjeto.vb[0],'vb[1]':editaobjeto.vb[1]}
         Fb_tabla={'Fb[0]':editaobjeto.Fb[0],'Fb[1]':editaobjeto.Fb[1]}
         mA_tabla={'mA[0]':editaobjeto.mA[0],'mA[1]':editaobjeto.mA[1]}
         interdist_tabla={'interdist[0]':0,'interdist[1]':0}
        #####termina la lista de las variables iniciales de los barcos#####
         self.barco=editaobjeto
         self.editarbarcos = MyBarco(self.listabarcos,editaobjeto,tablero,krpid_tabla,krd_tabla,krl_tabla,kvpid_tabla,
		                          ab_tabla,pb_tabla,vb_tabla,Fb_tabla,mA_tabla,interdist_tabla)
         self.editarbarcos.setGeometry(QRect(100, 100, 1250,550))
         self.editarbarcos.setWindowIcon(QIcon('images.ico'))
         self.editarbarcos.setWindowTitle('Editar Barcos')
         self.editarbarcos.show()
   
         
         
         self.combobox.setCurrentIndex(0)
       

    

        
         
     
#Con esta interfaz se hace el enlace, despues de ello ya no se puede volver atras ya que los barcos y la cadena ya están enlazados.
#		 
class MyCombinar(QWidget):
    def __init__(self,listacadenas,listabarcos,listasimulacion):
     QWidget.__init__(self)
     self.barcoderecho=None
     self.barcoizquierdo=None
     self.listasimulacion=listasimulacion
     self.listabarcos=listabarcos
     self.listacadenas=listacadenas
     self.cadena=None
     self.elonga=0.8
     self.derecho=0
     self.izquierdo=0
     self.centro=0
     self.b = QWidget(self)
     
     with open("style.css") as f:
          self.setStyleSheet(f.read())
     self.btn = QPushButton("Enlazar", self.b)
     self.connect(self.btn, SIGNAL("clicked()"), self.enlazar)
     self.btn.setStyleSheet("color:gray;")
     self.btn.setEnabled(False)
     self.btnelongacion = QPushButton(u"Elongación", self.b)
     self.connect(self.btnelongacion, SIGNAL("clicked()"), self.elongacion)
     self.btnRlink = QPushButton(u"Radio Enlace", self.b)
     self.connect(self.btnRlink, SIGNAL("clicked()"), self.radiolink)
    
     self.combobox2 = QComboBox(self.b)
     self.combobox2.addItem("Elegir Barco Derecho")
	 
     self.combobox3 = QComboBox(self.b)
     self.combobox3.addItem("Elegir Cadena Centro")
	 
     self.combobox = QComboBox(self.b)
     self.combobox.addItem("Elegir Barco a izquierdo")
     self.strlist = ['el identificador del barco es {}'.format(barcos) for barcos in listabarcos.keys()]
     self.strlist2 = ['el nombre de cadena es {}'.format(cad1) for cad1 in listacadenas.keys()]
     self.combobox3.addItems(self.strlist2)
     self.combobox3.activated.connect(self.escojercadenacentro)
     self.combobox.addItems(self.strlist)
     self.combobox.activated.connect(self.escojerbarcoizquierdo)
     self.combobox2.addItems(self.strlist)
     self.combobox2.activated.connect(self.escojerbarcoderecho)
     self.mostrar = QLabel()
     self.mostrar.setOpenExternalLinks(False)
     self.mostrarError = QLabel()
     self.mostrarError.setOpenExternalLinks(False)
     layout = QVBoxLayout(self)
     layout.addWidget(self.combobox)
     layout.addWidget(self.combobox3)
     layout.addWidget(self.combobox2)
     layout.addWidget(self.btnelongacion)
     layout.addWidget(self.btnRlink)
     layout.addWidget(self.btn)
     
     
     layout.addWidget(self.mostrar)
     layout.addWidget(self.mostrarError)
     self.setLayout(layout) 
    #la comunicación se hace a través de esta función que llama a la clase MyRadioLink
    def radiolink(self):
        
        
        self.radiolinkshow = MyRadioLink(self.barcoizquierdo,self.barcoderecho)
        self.radiolinkshow.setGeometry(QRect(100, 100, 400, 50))
        self.radiolinkshow.setWindowIcon(QIcon('images.ico'))
        self.radiolinkshow.setWindowTitle('Escojer Variables a Comunicar')
        self.radiolinkshow.show()
       
    #la elongación se escoge de 0 a 100 
    def elongacion(self):
        
        d,ok = QInputDialog.getInt(self,"Escojer Elongacion ","Elongacion %:",80,0,100,1)
		
        if (ok):
           self.elonga=float(d*0.01)
    #con esta función se enlaza los barcos con la barrera   
    def enlazar(self):
      try:
        if self.barcoizquierdo.ide==self.barcoderecho.ide:
           self.mostrar.setText("no puedes enlazar un mismo barco")
           self.mostrar.setOpenExternalLinks(True)
        else:
            
           self.cadena.normal[:,0] = [1,0]
           self.cadena.normal[:,-1] = [-1,0]
           self.cadena.para[:,0] = [0,1]
           self.cadena.para[:,-1] = [0,-1]
           self.cadena.calcms()
          

#controles de alineacion longitudinal
           self.barcoizquierdo.krl = [2000.,100.,1000]
           self.barcoderecho.krl = [2000.,100.,1000]
#controles de distancia transversal
           self.barcoizquierdo.krd = [1.5,0.01,0.01] #[1,5.,0.1]
           self.barcoderecho.krd = [1.5,0.01,0.01] #[1,5.,0.1]

#controles de rumbo standard  
           self.barcoizquierdo.krpid = [5.,10.,0.1] #[5.,10.,0.1]
           self.barcoderecho.krpid = [5.,10.,0.1] #[5.,10.,0.1]

#controles de velocidad standard
           self.barcoizquierdo.kvpid = [2000,100,500]
           self.barcoderecho.kvpid = [2000,100,500]

###############################################################################

#link-amos las dos radios. 
           self.barcoizquierdo.rlink[1] = self.barcoderecho.rlink[0]
           self.barcoderecho.rlink[1] = self.barcoizquierdo.rlink[0]

		
		
		
		
           self.barcoizquierdo.pb = self.cadena.cms[:,0]  + self.barcoizquierdo.ls/2.*np.array([np.cos(self.barcoizquierdo.theta),np.sin(self.barcoizquierdo.theta)]) + self.cadena.L*self.cadena.para[:,0] + [0, self.elonga * self.cadena.cbl[0]]
           self.barcoderecho.pb = self.cadena.cms[:,-1] + self.barcoderecho.ls/2.*np.array([np.cos(self.barcoderecho.theta),np.sin(self.barcoderecho.theta)]) - self.cadena.L*self.cadena.para[:,-1] + [0,self.elonga * self.cadena.cbl[1]]
           
           if self.derecho > self.izquierdo:
              self.combobox.removeItem(self.derecho)
              self.combobox.removeItem(self.izquierdo)
              self.combobox2.removeItem(self.derecho)
              self.combobox2.removeItem(self.izquierdo)
           else:
              self.combobox.removeItem(self.izquierdo)
              self.combobox.removeItem(self.derecho)
              self.combobox2.removeItem(self.izquierdo)
              self.combobox2.removeItem(self.derecho)
              
           self.combobox3.removeItem(self.centro)
           
           self.mostrar.setText("")
           self.mostrar.setOpenExternalLinks(False)
           self.mostrarError.setText("")
           self.mostrarError.setOpenExternalLinks(False)
           app.main.btnsimula.setEnabled(True)
           app.main.btnsimula.setStyleSheet("color:white;")
           self.close()
      except:
        self.mostrarError.setText("Debes Escojer dos barcos y una cadena")
        self.mostrarError.setOpenExternalLinks(True)
         
           
    #se escoge el barco derecho    
    def escojerbarcoderecho(self,x):
        
        if x!=0:
         self.derecho=x
         t=self.combobox2.currentText()
         t=t.replace('el identificador del barco es','')
         ide=int(t)
         self.barcoderecho=self.listabarcos[ide]
         self.listasimulacion["derecho"]=self.barcoderecho
		 
	#se escoge cadena centro	 
    def escojercadenacentro(self,x):
        
        if x!=0:
         self.centro=x
         t=self.combobox3.currentText()
         t=t.replace('el nombre de cadena es','')
         cade1=str(t).strip()
         self.cadena=self.listacadenas[cade1]
         self.listasimulacion["centro"]=self.cadena
       
	#se escoge el barco izquierdo
    def escojerbarcoizquierdo(self,x):
        
        if x!=0:
         self.izquierdo=x
         t=self.combobox.currentText()
         t=t.replace('el identificador del barco es','')
         ide=int(t)
         self.barcoizquierdo=self.listabarcos[ide]
         self.listasimulacion["izquierdo"]=self.barcoizquierdo
		 

#Con esta clase se crea la gráfica de velocidad, rumbo, distancia entre los barcos, tensión
class MyGraficos(QWidget):
   def __init__(self,listabarcoIsim,listabarcoDsim,listacadenasim,rumbo,step,tam,delta):
     
     QWidget.__init__(self)
     self.rumbo=rumbo
     self.tiempo_paso=delta*tam/step
     self.listabarcoIsim=listabarcoIsim
     self.listabarcoDsim=listabarcoDsim
     self.listacadenasim=listacadenasim
     self.b = QWidget(self)
     with open("style.css") as f:
          self.setStyleSheet(f.read())
     self.stacked_widget1 = QStackedWidget()
     self.stacked_widget = QStackedWidget()
     self.combobox2 = QComboBox(self.b)
     layout = QVBoxLayout(self)
     self.btns=[]
     self.lista=['Modulo de la Velocidad','Tension en los Extremos','Rumbo','Distancia Surge y Sway']
     for i in self.lista:
         self.combobox2.addItem(i)
         btn = QPushButton(' Estado' + i ,self.b)
         btn.clicked.connect(lambda state, x=i: self.button_pushed(x))
         self.btns.append(btn)
         self.stacked_widget.addWidget(btn)
     
         
     
     self.combobox2.activated.connect(self.stacked_widget.setCurrentIndex)
     layout.addWidget(self.combobox2)
     layout.addWidget(self.stacked_widget)
     self.setLayout(layout)
	
	 
   #con esta función se puede escoger la gráfica que queremos ver. 
   def button_pushed(self,num):
        
        self.b.figure = pl.figure("Grafica "+num)
        a=[]
        t=[]
        t=self.tiempo_paso*np.arange(len(self.listabarcoDsim[:-1,2]))
        if num=='Modulo de la Velocidad' :
           for i in zip(self.listabarcoDsim[:-1,2],self.listabarcoDsim[:-1,3]):
               a.append(norm(i))
          # t=self.tiempo_paso*np.arange(size(a))
           
           
           pl.plot(t,a,'g')
           
           
           a=[]
           for i in zip(self.listabarcoIsim[:-1,2],self.listabarcoIsim[:-1,3]):
               a.append(norm(i))
           pl.plot(t,a,'r')
           pl.xlabel('tiempo (s)')
           pl.ylabel('m/s')
        if num=='Tension en los Extremos' :
           for i in zip(self.listabarcoDsim[:-1,-2],self.listabarcoDsim[:-1,-1]):
               a.append(norm(i))
           pl.plot(t,a,'g')
           a=[]
           for i in zip(self.listabarcoIsim[:-1,-2],self.listabarcoIsim[:-1,-1]):
               a.append(norm(i))
           pl.plot(t,a,'r')
           pl.xlabel('tiempo (s)')
           pl.ylabel('N')
        if num=='Rumbo':
           
           pl.plot(t,self.listabarcoDsim[:-1,6],'g')
           pl.plot(t,self.listabarcoIsim[:-1,6],'r')
           pl.xlabel('tiempo (s)')
           pl.ylabel('rad')
        if num=='Distancia Surge y Sway':
           b=[]
           for i in zip(self.listabarcoIsim[:-1,0],self.listabarcoIsim[:-1,1],self.listabarcoDsim[:-1,0],self.listabarcoDsim[:-1,1]):
                a.append((i[2]-i[0])*np.cos(self.rumbo)+(i[3]-i[1])*np.sin(self.rumbo))
                b.append((i[2]-i[0])*np.sin(self.rumbo)-(i[3]-i[1])*np.cos(self.rumbo))
           pl.plot(t,a,'k')
           pl.plot(t,b,'r')
           pl.xlabel('tiempo (s)')
           pl.ylabel('m')
		
        pl.show()
        		



#con ésta clase se muestra en la ventana de los estados de la simulacion por cada iteración	 
class MyEstado(QWidget):
   def __init__(self,rango,listabarcoIsim,listabarcoDsim,listacadenasim):
     
     QWidget.__init__(self)
     self.rango=rango
     self.listabarcoIsim=listabarcoIsim
     self.listabarcoDsim=listabarcoDsim
     self.listacadenasim=listacadenasim
     self.b = QWidget(self)
     with open("style.css") as f:
          self.setStyleSheet(f.read())
     self.stacked_widget1 = QStackedWidget()
     self.stacked_widget = QStackedWidget()
     self.combobox2 = QComboBox(self.b)
     layout = QVBoxLayout(self)
     self.btns=[]
     for i in self.rango:
         self.combobox2.addItem("Estado" + str(i))
         btn = QPushButton('Button {} Estado'.format(i),self.b)
         btn.clicked.connect(lambda state, x=i: self.button_pushed(x))
         self.btns.append(btn)
         self.stacked_widget.addWidget(btn)
     
         
     
     self.combobox2.activated.connect(self.stacked_widget.setCurrentIndex)
     layout.addWidget(self.combobox2)
     layout.addWidget(self.stacked_widget)
     self.setLayout(layout)
	
	 
   #con este botón se muestra la grafica 
   def button_pushed(self,num):
        
        self.b.figure = pl.figure("Simulacion Barco Cadena Barco Estado"+str(num))
        self.b.ax = self.b.figure.add_subplot(111)
        app.main.simula.run_cargar(num,self.b.figure,pl,self.listacadenasim,self.listabarcoIsim,self.listabarcoDsim,self.b.ax)
        pl.axis("equal")
        pl.show()
        		
#con esta clase se realiza la simulación de los procesos de manera concurrente para que se puedan #parar.      	
class Worker(QObject):
    candado = QMutex()
    resultado = pyqtSignal(int)
    resultadoplot=pyqtSignal(int)
    iniciado = pyqtSignal()
    terminado = pyqtSignal(bool)

    def __init__(self,rango,discri):
        super(Worker,  self).__init__()
        self._trabajando = False
        self._abortar = False
        self.rango=rango
        self.discri=discri
        
        

    def iniciar(self):
        self.candado.lock()
        self._trabajando = True
        self._abortar = False
        self.candado.unlock()
        self.iniciado.emit()

    def abortar(self):
        self.candado.lock()
        if self._trabajando:
            self._abortar = True
        self.candado.unlock()

    def procesar(self):
        
        for n in self.rango:
            self.candado.lock()
            abortar = self._abortar
            self.candado.unlock()
            self.resultado.emit(n)
            self.resultadoplot.emit(n)
            if self.discri==1:
                time.sleep(1)
            else:
                time.sleep(0.001)
            if abortar:
                
                break
            
				
        else:
            print ''
           
                

        
        
        self.candado.lock()
        self._trabajando = False
        self.candado.unlock()
        self.terminado.emit(True)		
	 
	 
#Este abre la ventana para hacer una simulación con la configuración deseada.
#	 
class MySimula(QWidget):
    def __init__(self,lista_barcoD_simu,lista_barcoI_simu,lista_cadenaC_simu,listasimulacion,entorno):
	
     QWidget.__init__(self)
     self.simulacion = QWidget(self)
     with open("style.css") as f:
          self.setStyleSheet(f.read())
     self.listabarcoDsim=lista_barcoD_simu
     self.listabarcoIsim=lista_barcoI_simu
     self.listacadenasim=lista_cadenaC_simu
     self.listasimulacion=listasimulacion
     self.tiempo=0
     tablero_simula=None
     if entorno==None:
        tablero_simula={'step':2000,'delta':0.00001,'tam':200000,'d_sway':5,'d_surge':5,'d_rumbo':np.pi/4}
     else:
        tablero_simula={'step':entorno['step'],'delta':entorno['delta'],'tam':entorno['tam'],'d_sway':entorno['d_sway'],'d_surge':entorno['d_surge'],'d_rumbo':entorno['d_rumbo']}
     self.tabla=copy(tablero_simula)
     
     self.Etiqueta = QLabel()
     self.Etiqueta.setText("")
     self.Etiqueta.setOpenExternalLinks(False)
     
     
     
     self.btnguardar = QPushButton("Save", self.simulacion)
     self.connect(self.btnguardar, SIGNAL("clicked()"), self.file_save)
     self.btnguardar.setEnabled(False)
     self.btnguardar.setStyleSheet("color:gray;")
     self.btnrun = QPushButton("run", self.simulacion)
     self.connect(self.btnrun, SIGNAL("clicked()"), self.runa)
     self.btnparar = QPushButton("stop", self.simulacion)
     self.connect(self.btnparar, SIGNAL("clicked()"), self.parar)
     self.btnparar.setEnabled(False)
     self.btnparar.setStyleSheet("color:gray;")
     if entorno==None:
        self.btnrun.setEnabled(False)
        self.btnrun.setStyleSheet("color:gray;")
     else:
        self.btnrun.setEnabled(True)
        self.btnrun.setStyleSheet("color:white;")
     
     self.btnestado = QPushButton("run estado", self.simulacion)
     self.connect(self.btnestado, SIGNAL("clicked()"), self.runa_estado)
	 
	 
     self.btngrafico = QPushButton(u"Gráficos Adicionales", self.simulacion)
     self.connect(self.btngrafico, SIGNAL("clicked()"), self.run_grafico)
     if entorno==None:
        self.btnestado.setEnabled(False)
        self.btnestado.setStyleSheet("color:gray;")
        self.btngrafico.setEnabled(False)
        self.btngrafico.setStyleSheet("color:gray;")
     else:
        self.btnestado.setEnabled(True)
        self.btnestado.setStyleSheet("color:white;")
        self.btngrafico.setEnabled(True)
        self.btngrafico.setStyleSheet("color:white;")
     
     
     self.btncargardatos = QPushButton("Cargar", self.simulacion)
     self.connect(self.btncargardatos, SIGNAL("clicked()"), self.cargardatos)
     self.progress = QProgressBar(self.simulacion)
     self.btncargardatos.setEnabled(True)
     self.btncargardatos.setStyleSheet("color:white;")
     
     
     
     
     self.comboboxpasos= QComboBox(self.simulacion)
     self.comboboxpasos.addItem("Elegir el Modo Pasos")
        
     
        
     
     self.tablaQ=MyTable(self.tabla,6,3,None)
        
     
     strlistpasos = ["normal","intercalado 10 veces","intercalado 100 veces"]
     self.comboboxpasos.addItems(strlistpasos)
     self.comboboxpasos.activated.connect(self.modo)
     self.worker=None
     self.worker2=None	 
     if entorno==None:
        self.comboboxpasos.setEnabled(False)
        
        self.rango=None
     else:
        self.comboboxpasos.setEnabled(True)
        
        self.rango=range(int(self.tabla['tam'])+1)
	 
     
     app.main.btnsimula.setEnabled(False)
     app.main.btnsimula.setStyleSheet("color:gray;")
     app.main.btncmb.setEnabled(False)
     app.main.btncmb.setStyleSheet("color:gray;")
     app.main.btnbarcos.setEnabled(False)
     app.main.btnbarcos.setStyleSheet("color:gray;")
     app.main.btncadena.setEnabled(False)
     app.main.btncadena.setStyleSheet("color:gray;")
     app.main.btnbyc.setEnabled(False)
     app.main.btnbyc.setStyleSheet("color:gray;")
	  
     self.hilo = QThread()
     self.worker = Worker(self.rango,0)
     self.worker.moveToThread(self.hilo)
     self.worker.resultado.connect(self.actualizar_etiqueta)
     self.worker.resultadoplot.connect(self.dibujar_plot)
     self.worker.iniciado.connect(self.hilo.start)
     self.hilo.started.connect(self.worker.procesar)
     self.worker.terminado.connect(self.hilo.quit)
     
	 
	 
     self.hilo2 = QThread()
     self.worker2 = Worker(self.rango,0)
     self.worker2.moveToThread(self.hilo2)
     self.worker2.resultado.connect(self.actualizar_etiqueta)
     self.worker2.resultadoplot.connect(self.cargar_plot)
     self.worker2.iniciado.connect(self.hilo2.start)
     self.hilo2.started.connect(self.worker2.procesar)
     self.worker2.terminado.connect(self.hilo2.quit)

     
     
     layout = QVBoxLayout(self)
     layout.addWidget(self.tablaQ)
     layout.addWidget(self.btncargardatos)
     layout.addWidget(self.comboboxpasos)
     layout.addWidget(self.btnestado)
     layout.addWidget(self.btnrun)
     layout.addWidget(self.btngrafico)
     layout.addWidget(self.btnparar)
     layout.addWidget(self.progress)
     layout.addWidget(self.btnguardar)
     layout.addWidget(self.Etiqueta)
     self.setLayout(layout)
	 
	     
    def file_save(self):
        
        try:
          name = QFileDialog.getSaveFileName(self, 'Save File', filter='(*.npz)')
        
        
          np.savez(str(name),listabarcoDsim=self.listabarcoDsim,listabarcoIsim=self.listabarcoIsim,listacadenasim=self.listacadenasim,listabarcoDsim2=self.listabarcoDsim2,listabarcoIsim2=self.listabarcoIsim2,listacadenasim2=self.listacadenasim2,tam=int(self.tabla['tam']),step=int(self.tabla['step']),delta=float(self.tabla['delta']),d_sway=float(self.tabla['d_sway']),d_surge=float(self.tabla['d_surge']),d_rumbo=float(self.tabla['d_rumbo']))
          self.Etiqueta.setText("")
          self.Etiqueta.setOpenExternalLinks(False)
          self.close()
		
        except:
          self.Etiqueta.setText("Debes escribir donde guardar")
          self.Etiqueta.setOpenExternalLinks(True)
		
    def closeEvent(self, event):
        
        event.accept()
        app.main.btnbarcos.setEnabled(True)
        app.main.btnbarcos.setStyleSheet("color:white;")
        app.main.btncadena.setEnabled(True)
        app.main.btncadena.setStyleSheet("color:white;")
        app.main.lista_barcos=dict()
        app.main.lista_cadenas=dict()
        
        
       
        		
	
    def cargar_plot(self,k):
       
       
       self.cargar(k)	   
       if k==int(self.tabla['tam']):	
             self.listabarcoIsim[-1,0:12] = self.listasimulacion["izquierdo"].thewjmax,self.listasimulacion["izquierdo"].Ac,self.listasimulacion["izquierdo"].mb,self.listasimulacion["izquierdo"].mul,self.listasimulacion["izquierdo"].Ib,self.listasimulacion["izquierdo"].mua,self.listasimulacion["izquierdo"].ls,self.listasimulacion["izquierdo"].mut,self.listasimulacion["izquierdo"].pmax,self.listasimulacion["izquierdo"].pmax,self.listasimulacion["izquierdo"].Ab,self.listasimulacion["izquierdo"].Aw
			 
             self.listabarcoIsim2[-2,0:19] = self.listasimulacion["izquierdo"].ide,self.listasimulacion["izquierdo"].tiempo,self.listasimulacion["izquierdo"].alfab,self.listasimulacion["izquierdo"].wb,self.listasimulacion["izquierdo"].theta,self.listasimulacion["izquierdo"].Fm,self.listasimulacion["izquierdo"].setl,self.listasimulacion["izquierdo"].setw,self.listasimulacion["izquierdo"].interrv,self.listasimulacion["izquierdo"].interrumbo,self.listasimulacion["izquierdo"].thewj,self.listasimulacion["izquierdo"].M,self.listasimulacion["izquierdo"].mul2,self.listasimulacion["izquierdo"].mut2,self.listasimulacion["izquierdo"].mua2,self.listasimulacion["izquierdo"].Tsf,self.listasimulacion["izquierdo"].Tsw,self.listasimulacion["izquierdo"].cfr,self.listasimulacion["izquierdo"].link
			 
			 
             self.listabarcoIsim2[-3,0:2] = self.listasimulacion["izquierdo"].ab
             self.listabarcoIsim2[-4,0:3] = self.listasimulacion["izquierdo"].krpid
             self.listabarcoIsim2[-5,0:3] = self.listasimulacion["izquierdo"].krd
             self.listabarcoIsim2[-6,0:3] = self.listasimulacion["izquierdo"].krl
             self.listabarcoIsim2[-7,0:3] = self.listasimulacion["izquierdo"].kvpid
             self.listabarcoIsim2[-8,0:2] = self.listasimulacion["izquierdo"].ab
             self.listabarcoIsim2[-9,0:2] = self.listasimulacion["izquierdo"].pb
             self.listabarcoIsim2[-10,0:2] = self.listasimulacion["izquierdo"].vb
             self.listabarcoIsim2[-11,0:2] = self.listasimulacion["izquierdo"].Fb
             self.listabarcoIsim2[-12,0:2] = self.listasimulacion["izquierdo"].mA
             self.listabarcoIsim2[-13,0:2] = self.listasimulacion["izquierdo"].interdist
			 
			 
             self.listabarcoDsim[-1,0:12] = self.listasimulacion["derecho"].thewjmax,self.listasimulacion["derecho"].Ac,self.listasimulacion["derecho"].mb,self.listasimulacion["derecho"].mul,self.listasimulacion["derecho"].Ib,self.listasimulacion["derecho"].mua,self.listasimulacion["derecho"].ls,self.listasimulacion["derecho"].mut,self.listasimulacion["derecho"].pmax,self.listasimulacion["derecho"].pmax,self.listasimulacion["derecho"].Ab,self.listasimulacion["derecho"].Aw
			 
             self.listabarcoDsim2[-2,0:19] = self.listasimulacion["derecho"].ide,self.listasimulacion["derecho"].tiempo,self.listasimulacion["derecho"].alfab,self.listasimulacion["derecho"].wb,self.listasimulacion["derecho"].theta,self.listasimulacion["derecho"].Fm,self.listasimulacion["derecho"].setl,self.listasimulacion["derecho"].setw,self.listasimulacion["derecho"].interrv,self.listasimulacion["derecho"].interrumbo,self.listasimulacion["derecho"].thewj,self.listasimulacion["derecho"].M,self.listasimulacion["derecho"].mul2,self.listasimulacion["derecho"].mut2,self.listasimulacion["derecho"].mua2,self.listasimulacion["derecho"].Tsf,self.listasimulacion["derecho"].Tsw,self.listasimulacion["derecho"].cfr,self.listasimulacion["derecho"].link
			 
			 
             self.listabarcoDsim2[-3,0:2] = self.listasimulacion["derecho"].ab
             self.listabarcoDsim2[-4,0:3] = self.listasimulacion["derecho"].krpid
             self.listabarcoDsim2[-5,0:3] = self.listasimulacion["derecho"].krd
             self.listabarcoDsim2[-6,0:3] = self.listasimulacion["derecho"].krl
             self.listabarcoDsim2[-7,0:3] = self.listasimulacion["derecho"].kvpid
             self.listabarcoDsim2[-8,0:2] = self.listasimulacion["derecho"].ab
             self.listabarcoDsim2[-9,0:2] = self.listasimulacion["derecho"].pb
             self.listabarcoDsim2[-10,0:2] = self.listasimulacion["derecho"].vb
             self.listabarcoDsim2[-11,0:2] = self.listasimulacion["derecho"].Fb
             self.listabarcoDsim2[-12,0:2] = self.listasimulacion["derecho"].mA
             self.listabarcoDsim2[-13,0:2] = self.listasimulacion["derecho"].interdist
			 
             self.listacadenasim[-1,0,0:3] = self.listasimulacion["centro"].s,self.listasimulacion["centro"].q,self.listasimulacion["centro"].A
             self.listacadenasim[-1,1,0:3] = self.listasimulacion["centro"].L,self.listasimulacion["centro"].m,self.listasimulacion["centro"].I
             self.listacadenasim[-1,2,0:2] = float(self.tabla['delta']),int(self.tabla['tam'])
             self.listacadenasim[-1,3,0:3] = self.listasimulacion["centro"].cbl
             self.listacadenasim2[-1,4,0:2] = self.listasimulacion["centro"].Anch
             self.listacadenasim2[-1,5,0:2] = self.listasimulacion["centro"].Fd
             self.listacadenasim2[-1,6,0:2] = self.listasimulacion["centro"].Fi
             self.listacadenasim2[-1,7,0:4] = self.listasimulacion["centro"].esl,self.listasimulacion["centro"].s2,self.listasimulacion["centro"].q2,self.listasimulacion["centro"].mA
             self.btnparar.setEnabled(False)
             self.btnparar.setStyleSheet("color:gray;") 
             self.btnrun.setEnabled(True)
             self.btnrun.setStyleSheet("color:white;") 
             self.btnestado.setEnabled(True)
             self.btnestado.setStyleSheet("color:white;") 
             self.comboboxpasos.setEnabled(True)			 
             self.btncargardatos.setEnabled(True)
             self.btncargardatos.setStyleSheet("color:white;") 
             self.btnguardar.setEnabled(True)
             self.btnguardar.setStyleSheet("color:white;") 
             self.tablaQ.setEnabled(True)
             			 
    
       
    def dibujar_plot(self,i):
        self.ax = self.simulacion.figure.add_subplot(111)
        if i==self.rango[len(self.rango)-1]:
              self.btnguardar.setEnabled(True)
              self.btnguardar.setStyleSheet("color:white;")		
              self.btnparar.setEnabled(False)
              self.btnparar.setStyleSheet("color:gray;") 
              self.btnrun.setEnabled(True)
              self.btnrun.setStyleSheet("color:white;") 
              self.btnestado.setEnabled(True)
              self.btnestado.setStyleSheet("color:white;") 
              self.comboboxpasos.setEnabled(True)
               
              self.btncargardatos.setEnabled(True)
              self.btncargardatos.setStyleSheet("color:white;") 
              self.btngrafico.setEnabled(True)
              self.btngrafico.setStyleSheet("color:white;") 
              self.tablaQ.setEnabled(True)
        self.run_cargar(i,self.simulacion.figure,pl,self.listacadenasim,self.listabarcoIsim,self.listabarcoDsim,self.ax)
        pl.show()
        pl.pause(0.1)
        pl.axis("equal")
        pl.hold(True)
	 
    def actualizar_etiqueta(self,l):
       
       self.progress.setValue(l)
	  
    def parar(self):
       
        self.worker.abortar()
        self.btnparar.setEnabled(False)
        self.btnparar.setStyleSheet("color:gray;")
        self.tablaQ.setEnabled(True)
        self.btncargardatos.setEnabled(True)
        self.btncargardatos.setStyleSheet("color:white;") 
        if(int(self.tabla['tam']))!=self.rango[len(self.rango)-1]:
                    self.btnrun.setEnabled(True)
                    self.btnrun.setStyleSheet("color:white;") 
                    
                    
        else:
                    
                    self.btnguardar.setEnabled(True)
                    self.btnguardar.setStyleSheet("color:white;") 
                    
        pl.clf()
        self.worker2.abortar()
        self.Etiqueta.setText("Proceso Parado")
        self.Etiqueta.setOpenExternalLinks(True)
		
		
    
    def run_grafico(self):
       
         
       grafico = MyGraficos(self.listabarcoIsim,self.listabarcoDsim,self.listacadenasim,float(self.tabla['d_rumbo']),
                            float(self.tabla['step']),float(self.tabla['tam']),float(self.tabla['delta']))
       grafico.setGeometry(QRect(100, 100, 400, 50))
       grafico.setWindowIcon(QIcon('images.ico'))
       grafico.setWindowTitle('Esquemas Adicionales')
       grafico.show()
    
    def runa_estado(self):
       
         
       if self.rango==range(int(self.tabla['tam'])+1): self.rango=range(self.listabarcoIsim.shape[0]-1)
       editarcadena = MyEstado(self.rango,self.listabarcoIsim,self.listabarcoDsim,self.listacadenasim)
       editarcadena.setGeometry(QRect(100, 100, 400, 50))
       editarcadena.setWindowIcon(QIcon('images.ico'))
       editarcadena.setWindowTitle('Pintar Estado')
       editarcadena.show()
       
	
    def modo(self,x):
        
        if x!=0:
             if x==1:
              self.rango=range(self.listabarcoIsim.shape[0]-1)
             
             if x==2:
              self.rango=range(0,self.listabarcoIsim.shape[0]-1,10)
              
             if x==3:
              self.rango=range(0,self.listabarcoIsim.shape[0]-1,100)
              
        self.comboboxpasos.setCurrentIndex(0)	  
		
        
    
                 
		
         
    def cargardatos(self):
        
        self.listabarcoIsim = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,16))
        self.listabarcoDsim = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,16))
        self.listacadenasim = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,14,self.listasimulacion["centro"].esl))
        self.listabarcoIsim2 = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,19))
        self.listabarcoDsim2 = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,19))
        self.listacadenasim2 = np.zeros((int(self.tabla['tam'])//int(self.tabla['step'])+2,14,self.listasimulacion["centro"].esl))
        self.rango=range(int(self.tabla['tam'])+1)
        setattr(self.worker2,'rango',self.rango)
        self.progress.setRange(0,int(self.tabla['tam'])) 
        self.btnguardar.setEnabled(False)
        self.btnguardar.setStyleSheet("color:gray;")		
        self.btncargardatos.setEnabled(False)
        self.btncargardatos.setStyleSheet("color:gray;")
        self.btnrun.setEnabled(False)
        self.btnrun.setStyleSheet("color:gray;")
        self.btnestado.setEnabled(False)
        self.btnestado.setStyleSheet("color:gray;")
        self.comboboxpasos.setEnabled(False)
        
        self.tablaQ.setEnabled(False)
        self.btnparar.setEnabled(True)
        self.btnparar.setStyleSheet("color:white;")
        self.btngrafico.setEnabled(False)
        self.worker2.iniciar()
        self.Etiqueta.setText("")
        self.Etiqueta.setOpenExternalLinks(False)  
        
	
    def cargar(self,i):
          
###############################################################################
#EMPIEZA BUCLE CALCULO
#########00#####################################################################
        self.listasimulacion["izquierdo"].propulsion(float(self.tabla['delta']),extf = self.listasimulacion["centro"].T[0:2],dfcm = [-self.listasimulacion["izquierdo"].ls/2.,0])
        self.listasimulacion["derecho"].propulsion(float(self.tabla['delta']),extf = - self.listasimulacion["centro"].T[-2:],dfcm = [-self.listasimulacion["derecho"].ls/2.,0])
        self.listasimulacion["centro"].mtr_s(self.listasimulacion["izquierdo"],self.listasimulacion["derecho"])
        self.listasimulacion["centro"].movimiento(delta = float(self.tabla['delta']))
        self.listasimulacion["izquierdo"].movimiento(delta = float(self.tabla['delta']) )
        self.listasimulacion["derecho"].movimiento(delta = float(self.tabla['delta']))
        self.listasimulacion["izquierdo"].controlador(float(self.tabla['delta']),float(self.tabla['d_rumbo']),1.,1000,np.array([float(self.tabla['d_surge']),-float(self.tabla['d_sway'])]))
        self.listasimulacion["derecho"].controlador(float(self.tabla['delta']),float(self.tabla['d_rumbo']),1.,1000,np.array([-float(self.tabla['d_surge']),float(self.tabla['d_sway'])]))
        self.tiempo= self.tiempo + float(self.tabla['delta'])
        
        ########################Parche para acelerar##################        
        self.listasimulacion["centro"].calcms(self.listasimulacion["centro"].cms[:,0],self.listasimulacion["centro"].cms[:,-1],-1*self.listasimulacion["centro"].ord)
        ##############################################################
        
			 ############GUARDAMOS PASOS DE PINTADO#########
			 ###############################################
        if (i%int(self.tabla['step']) == 0):
            self.listabarcoIsim[i//int(self.tabla['step'])] = self.listasimulacion["izquierdo"].pb[0],self.listasimulacion["izquierdo"].pb[1],self.listasimulacion["izquierdo"].vb[0],self.listasimulacion["izquierdo"].vb[1],self.listasimulacion["izquierdo"].ab[0],self.listasimulacion["izquierdo"].ab[1],self.listasimulacion["izquierdo"].theta,self.listasimulacion["izquierdo"].wb,self.listasimulacion["izquierdo"].alfab,self.listasimulacion["izquierdo"].Fm,self.listasimulacion["izquierdo"].M,self.listasimulacion["izquierdo"].setl,self.listasimulacion["izquierdo"].setw,self.listasimulacion["izquierdo"].thewj,self.listasimulacion["centro"].T[0],self.listasimulacion["centro"].T[1]

            self.listabarcoDsim[i//int(self.tabla['step'])] = self.listasimulacion["derecho"].pb[0],self.listasimulacion["derecho"].pb[1],self.listasimulacion["derecho"].vb[0],self.listasimulacion["derecho"].vb[1],self.listasimulacion["derecho"].ab[0],self.listasimulacion["derecho"].ab[1],self.listasimulacion["derecho"].theta,self.listasimulacion["derecho"].wb,self.listasimulacion["derecho"].alfab,self.listasimulacion["derecho"].Fm,self.listasimulacion["derecho"].M,self.listasimulacion["derecho"].setl,self.listasimulacion["derecho"].setw,self.listasimulacion["derecho"].thewj,self.listasimulacion["centro"].T[-2],self.listasimulacion["centro"].T[-1]

            self.listacadenasim[i//int(self.tabla['step'])] = self.listasimulacion["centro"].cms[0],self.listasimulacion["centro"].cms[1],self.listasimulacion["centro"].v[0],self.listasimulacion["centro"].v[1],self.listasimulacion["centro"].a[0],self.listasimulacion["centro"].a[1],self.listasimulacion["centro"].alfa,self.listasimulacion["centro"].w,self.listasimulacion["centro"].T[2:-1:2],self.listasimulacion["centro"].T[2::2],self.listasimulacion["centro"].normal[0],self.listasimulacion["centro"].normal[1],self.listasimulacion["centro"].para[0],self.listasimulacion["centro"].para[1]
       		
        
         
        
	
	
	
	
	
    def runa(self):
        self.simulacion.figure = pl.figure("Simular Barco Cadena Barco")		
        if self.rango==range(int(self.tabla['tam'])+1) or self.rango==None: self.rango=range(self.listabarcoIsim.shape[0]-1)
        setattr(self.worker,'rango',self.rango)
        setattr(self.worker,'discri',1)
        self.progress.setRange(0,self.listabarcoIsim.shape[0]-2)
        pl.clf()
        self.btncargardatos.setEnabled(False)
        self.btncargardatos.setStyleSheet("color:gray;")
        self.btnestado.setEnabled(False)
        self.btnestado.setStyleSheet("color:gray;")
        self.comboboxpasos.setEnabled(False)
        self.comboboxpasos.setStyleSheet("color:gray;")
        self.btnrun.setEnabled(False)
        self.btnrun.setStyleSheet("color:gray;")
        self.tablaQ.setEnabled(False)
        self.btnparar.setEnabled(True)
        self.btnparar.setStyleSheet("color:white;")
        self.btnguardar.setEnabled(False)
        self.btnguardar.setStyleSheet("color:gray;")		
        self.worker.iniciar()
        self.Etiqueta.setText("")
        self.Etiqueta.setOpenExternalLinks(False)
		
        
        
        
	
    def run_cargar(self,i,figure,pl,listacadenasim,listabarcoIsim,listabarcoDsim,ax):
        
        canvasimula = FigureCanvas(figure)
        #canvasimula.mpl_connect('close_event', handle_close)
        
        #ax.cla()
        ax.hold(True)
        cms = np.array([listacadenasim[i,0,:],listacadenasim[i,1,:]])
        para = np.array([listacadenasim[i,-2,:],listacadenasim[i,-1,:]])
        ax.plot(cms[0,:],cms[1,:],'o')
               
        barrasi = cms + listacadenasim[-1,1,0] * para
        barrasd = cms - listacadenasim[-1,1,0] * para
        ax.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
              
        ax.plot(listabarcoIsim[i,0],listabarcoIsim[i,1],'+r') #r
               
        ax.plot(listabarcoDsim[i,0],listabarcoDsim[i,1],'+g') #g
              
####################Dibujar cable de arrastre################################
        rot = np.array([[np.cos(listabarcoIsim[i,6]),- np.sin(listabarcoIsim[i,6])],[np.sin(listabarcoIsim[i,6]),np.cos(listabarcoIsim[i,6])]])  
        popai =  np.dot(rot, np.array([-listabarcoIsim[-1,6]/2.,0])) + [listabarcoIsim[i,0],listabarcoIsim[i,1]]  
        tipi = para[:,0] * listacadenasim[-1,1,0] + cms[:,0]
        disti = norm(popai - tipi)
        di = disti/listacadenasim[-1,3,0]
        if di > 1: di = 1
        r = bz.bezier4p([[tipi[0]],[tipi[1]]],[[popai[0]],[popai[1]]],1,1,1.5,(1-di) * listabarcoIsim[i,6]+di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),(1-di) * np.arctan2(-para[0,0],-para[0,1])+di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),100)
        ax.plot(r[0][0],r[0][1],'b')
              
               

        
#####################################Dibujar Estado Barco ############################################Izquierdo##########################################        
        vertices = np.array([[-listabarcoIsim[-1,6]/2.,-0.25],[-listabarcoIsim[-1,6]/2.,0.25],[-0.25,0.35],[listabarcoIsim[-1,6]/2.,0],[-0.25,-0.35],[-listabarcoIsim[-1,6]/2.,-0.25]])        
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [listabarcoIsim[i,0],listabarcoIsim[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3]
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'red') #'red'
               
        pl.gca().add_patch(patchi)
               

#######################dibujar cable de arrastre derecha#######################
        rot = np.array([[np.cos(listabarcoDsim[i,6]),- np.sin(listabarcoDsim[i,6])],[np.sin(listabarcoDsim[i,6]),np.cos(listabarcoDsim[i,6])]])  
        popad =  np.dot(rot, np.array([-listabarcoDsim[-1,6]/2.,0])) + [listabarcoDsim[i,0],listabarcoDsim[i,1]]  
        tipd = - para[:,-1] * listacadenasim[-1,1,0] + cms[:,-1]
        distd = norm(popad - tipd)
        dd = distd/listacadenasim[-1,3,0]
        if dd > 1: dd = 1
        r = bz.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,(1-dd) * listabarcoDsim[i,6]+dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),(1-dd) * np.arctan2(-para[0,0],-para[0,1])+dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),100)
        ax.plot(r[0][0],r[0][1],'b')
        #ax.hold(True)
               
#############################Dibujar Estado Barco ########################################Derecho##################################################

        vertices = np.array([[-listabarcoDsim[-1,6]/2.,-0.25],[-listabarcoDsim[-1,6]/2.,0.25],[-0.25,0.35],[listabarcoDsim[-1,6]/2.,0],[-0.25,-0.35],[-listabarcoDsim[-1,6]/2.,-0.25]])              
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [listabarcoDsim[i,0],listabarcoDsim[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3]
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
        
        pl.gca().add_patch(patchd)
        
        ax.plot(listabarcoIsim[i,0],listabarcoIsim[i,1],'+r')
        
        ax.plot(listabarcoDsim[i,0],listabarcoDsim[i,1],'+g')
        ax.set_xlabel('m')
        ax.set_ylabel('m')

        canvasimula.draw()
               


           
        
#clase principal de la Interfaz      	
class MainWindow(QMainWindow):
    def __init__(self, *args):
        QMainWindow.__init__(self, *args)
        with open("style.css") as f:
            self.setStyleSheet(f.read())
        self.lista_barcos=dict()
        self.lista_cadenas=dict()
        self.lista_barcoD_simu=None
        self.lista_barcoI_simu=None
        self.lista_cadenaC_simu=None
        self.lista_simulacion=dict()
        self.principal = QWidget(self)
        self.setGeometry(300, 300, 350, 100)
        layout = QVBoxLayout(self.principal)
        self.setCentralWidget(self.principal)
        self.btnbarcos = QPushButton("Crear Barcos", self.principal)
        self.btnbyc = QPushButton("Editar Barcos y Cadenas", self.principal)
        self.btnbyc.setEnabled(False)
        self.btnbyc.setStyleSheet("color:gray;")
        self.btncadena = QPushButton("Crear Cadena", self.principal)
        self.btncmb = QPushButton("Combinar", self.principal)
        self.btncmb.setEnabled(False)
        self.btncmb.setStyleSheet("color:gray;")
        self.btnsimula = QPushButton("Simular Barco y Cadena", self.principal)
        self.btnsimula.setEnabled(False)
        self.btnsimula.setStyleSheet("color:gray;")
        self.setWindowIcon(QIcon('images.ico'))
        self.setWindowTitle('Interfaz Barcos')
        self.connect(self.btnbarcos, SIGNAL("clicked()"), self.barcocrear)
        self.connect(self.btnbyc, SIGNAL("clicked()"), self.editarbarcoycadena)
        self.connect(self.btncadena, SIGNAL("clicked()"), self.cadenacrear)
        self.connect(self.btncmb, SIGNAL("clicked()"), self.combinar)
        self.connect(self.btnsimula, SIGNAL("clicked()"), self.simular)
		
        self.file_menu = QMenu("File",self.principal)
     
        openFile = QAction(u"&Cargar Datos para una Simulación", self.principal)
        openFile.setStatusTip('Open File')
        openFile.triggered.connect(self.file_open)
        self.menuBar().addMenu(self.file_menu)
        
        
        self.file_menu.addAction(openFile)
		
        self.mostrar = QLabel()
        self.mostrar.setStyleSheet("font-size: 20px; color: red;")
        self.mostrar.setText("")
        self.mostrar.setOpenExternalLinks(False)
		
        
        layout.addWidget(self.btnbyc)
        layout.addWidget(self.btnbarcos)
        layout.addWidget(self.btncadena)
        layout.addWidget(self.btncmb)
        layout.addWidget(self.btnsimula)
        layout.addWidget(self.mostrar)
        self.principal.setLayout(layout)
        self.combina=None
        self.creacadena=None
        self.creabarcos=None
        self.editabyc=None
        self.muestrabyc=None
        self.simula=None
        self.creaobjeto=None
		
    #con ésta se abre la ventana para cargar los ficheros para una simulación.
    def file_open(self):
    
        try:
          name = QFileDialog.getOpenFileName(self, 'Open File', filter='(*.npz)')
          data = np.load(str(name))
          self.lista_barcoD_simu=data['listabarcoDsim']
          self.lista_barcoI_simu=data['listabarcoIsim']
          self.lista_cadenaC_simu=data['listacadenasim']
          self.lista_barcoD_simu2=data['listabarcoDsim2']
          self.lista_barcoI_simu2=data['listabarcoIsim2']
          self.lista_cadenaC_simu2=data['listacadenasim2']
          lista_entorno=dict()
          lista_entorno['step']=data['step']
          lista_entorno['delta']=data['delta']
          lista_entorno['tam']=data['tam']
          lista_entorno['d_sway']=data['d_sway']
          lista_entorno['d_surge']=data['d_surge']
          lista_entorno['d_rumbo']=data['d_rumbo']
          self.btnsimula.setEnabled(True)
        
		
          pb=self.lista_barcoD_simu2[-9,0:2]
          vb=self.lista_barcoD_simu2[-10,0:2]
          ab=self.lista_barcoD_simu2[-8,0:2]
          Fb=self.lista_barcoD_simu2[-11,0:2]
          mA=self.lista_barcoD_simu2[-12,0:2]
          interdist=self.lista_barcoD_simu2[-13,0:2]
          kvpid=self.lista_barcoD_simu2[-7,0:3]
          krpid=self.lista_barcoD_simu2[-4,0:3]
          krd=self.lista_barcoD_simu2[-5,0:3]
          krl=self.lista_barcoD_simu2[-6,0:3]
		
          self.lista_simulacion['derecho'] = usv.barco(int(self.lista_barcoD_simu2[-2,0]),"waterjet",self.lista_barcoD_simu2[-2,1],pb,vb,ab,self.lista_barcoD_simu2[-2,2],self.lista_barcoD_simu2[-2,3],self.lista_barcoD_simu2[-2,4],self.lista_barcoD_simu2[-2,5],Fb,self.lista_barcoD_simu2[-2,6],self.lista_barcoD_simu2[-2,7],self.lista_barcoD_simu2[-2,8],self.lista_barcoD_simu2[-2,9],self.lista_barcoD_simu2[-2,10],self.lista_barcoD_simu[-1,0],self.lista_barcoD_simu[-1,1],self.lista_barcoD_simu2[-2,11],self.lista_barcoD_simu[-1,2],mA,self.lista_barcoD_simu[-1,3],self.lista_barcoD_simu2[-2,12],self.lista_barcoD_simu[-1,7],self.lista_barcoD_simu2[-2,13],self.lista_barcoD_simu[-1,4],self.lista_barcoD_simu[-1,5],self.lista_barcoD_simu2[-2,14],int(self.lista_barcoD_simu[-1,6]),self.lista_barcoD_simu[-1,8],self.lista_barcoD_simu2[-2,15],self.lista_barcoD_simu2[-2,16],self.lista_barcoD_simu2[-2,17],interdist,self.lista_barcoD_simu2[-2,18],kvpid,krpid,krd,krl)
		
        
          pb=self.lista_barcoI_simu2[-9,0:2]
          vb=self.lista_barcoI_simu2[-10,0:2]
          ab=self.lista_barcoI_simu2[-8,0:2]
          Fb=self.lista_barcoI_simu2[-11,0:2]
          mA=self.lista_barcoI_simu2[-12,0:2]
          interdist=self.lista_barcoI_simu2[-13,0:2]
          kvpid=self.lista_barcoI_simu2[-7,0:3]
          krpid=self.lista_barcoI_simu2[-4,0:3]
          krd=self.lista_barcoI_simu2[-5,0:3]
          krl=self.lista_barcoI_simu2[-6,0:3]
		
          self.lista_simulacion['izquierdo']=usv.barco(int(self.lista_barcoI_simu2[-2,0]),"waterjet",self.lista_barcoI_simu2[-2,1],pb,vb,ab,self.lista_barcoI_simu2[-2,2],self.lista_barcoI_simu2[-2,3],self.lista_barcoI_simu2[-2,4],self.lista_barcoI_simu2[-2,5],Fb,self.lista_barcoI_simu2[-2,6],self.lista_barcoI_simu2[-2,7],self.lista_barcoI_simu2[-2,8],self.lista_barcoI_simu2[-2,9],self.lista_barcoI_simu2[-2,10],self.lista_barcoI_simu[-1,0],self.lista_barcoI_simu[-1,1],self.lista_barcoI_simu2[-2,11],self.lista_barcoI_simu[-1,2],mA,self.lista_barcoI_simu[-1,3],self.lista_barcoI_simu2[-2,12],self.lista_barcoI_simu[-1,7],self.lista_barcoI_simu2[-2,13],self.lista_barcoI_simu[-1,4],self.lista_barcoI_simu[-1,5],self.lista_barcoI_simu2[-2,14],int(self.lista_barcoI_simu[-1,6]),self.lista_barcoI_simu[-1,8],self.lista_barcoI_simu2[-2,15],self.lista_barcoI_simu2[-2,16],self.lista_barcoI_simu2[-2,17],interdist,self.lista_barcoI_simu2[-2,18],kvpid,krpid,krd,krl)
		
          self.lista_simulacion['izquierdo'].rlink = [[0,pb,self.lista_barcoI_simu2[-2,4],vb],[]]
          self.lista_simulacion['derecho'].rlink = [[1,pb,self.lista_barcoD_simu2[-2,4],vb],[]]
		
          self.lista_simulacion['izquierdo'].rlink[1] = self.lista_simulacion['derecho'].rlink[0]
          self.lista_simulacion['derecho'].rlink[1] = self.lista_simulacion['derecho'].rlink[0]
        
          cable=self.lista_cadenaC_simu[-1,3,0:3]
          Anch=self.lista_cadenaC_simu2[-1,4,0:2]
          Fd=self.lista_cadenaC_simu2[-1,5,0:2]
          Fi=self.lista_cadenaC_simu2[-1,6,0:2]
		
		
          self.lista_simulacion['centro']=boom.cadena(int(self.lista_cadenaC_simu2[-1,7,0]),cable,Anch,"cadena",self.lista_cadenaC_simu[-1,0,0],self.lista_cadenaC_simu[-1,0,1],self.lista_cadenaC_simu2[-1,7,1],self.lista_cadenaC_simu2[-1,7,2],self.lista_cadenaC_simu[-1,0,2],self.lista_cadenaC_simu[-1,1,0],self.lista_cadenaC_simu[-1,1,1],self.lista_cadenaC_simu2[-1,7,3],Fd,Fi)
			
        		
          self.lista_simulacion['centro'].cms[0]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][0]
          self.lista_simulacion['centro'].cms[1]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][1]
          self.lista_simulacion['centro'].v[0]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][2]
          self.lista_simulacion['centro'].v[1]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][3]
          self.lista_simulacion['centro'].a[0]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][4]
          self.lista_simulacion['centro'].a[1]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][5]
          self.lista_simulacion['centro'].alfa=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][6]
          self.lista_simulacion['centro'].w=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][7]
          self.lista_simulacion['centro'].T[2:-1:2]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][8]
          self.lista_simulacion['centro'].T[2::2]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][9]
          self.lista_simulacion['centro'].normal[0]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][10]
          self.lista_simulacion['centro'].normal[1]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][11]
          self.lista_simulacion['centro'].para[0]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][12]
          self.lista_simulacion['centro'].para[1]=self.lista_cadenaC_simu[lista_entorno['tam']//lista_entorno['step']][13]
		
		
          self.simula = MySimula(self.lista_barcoD_simu,self.lista_barcoI_simu,self.lista_cadenaC_simu,self.lista_simulacion,lista_entorno)
        
          self.simula.setWindowIcon(QIcon('images.ico'))
          self.simula.setWindowTitle('Carga de Simular Barco con Cadena')
          self.simula.show()
        except:
           
            self.mostrar.setText("Debes de escojer un archivo")
            self.mostrar.setOpenExternalLinks(False)
        
    #función que abre la ventana para enlazar los usvs con la barrera. 
    def combinar(self):
        self.combina = MyCombinar(self.lista_cadenas,self.lista_barcos,self.lista_simulacion)
        self.combina.setGeometry(QRect(100, 100, 400, 50))
        self.combina.setWindowIcon(QIcon('images.ico'))
        self.combina.setWindowTitle('Combinar Objetos')
        self.combina.show()	
     
    #esta función abre la ventana para crear cadena	 
    def cadenacrear(self):
        #####listas de las variables iniciales de las cadenas#########
        cadena = {'esl':22,'objeto':'cadena','s':50.,'q':0.,
                  's2':0.,'q2':0.,'A':2.,'L':0.5,'m':2.,'mA':0.}
        cable_tabla={'cable[0]':0.1, 'cable[1]':0.1, 'cable[2]':1000}
        Anch_tabla={'Anch[0]':0.,'Anch[1]':0.}
        Fd_tabla={'Fd[0]':0.0,'Fd[1]':0.0}
        Fi_tabla={'Fi[0]':0.0,'Fi[1]':0.0}
		######termina la lista de las variables iniciales de la cadena#####
        self.es_Cadena=False
        self.es_Barco=False
        self.creacadena = MyCadena(self.lista_cadenas,None,cadena,cable_tabla,Anch_tabla,Fd_tabla,Fi_tabla)
        self.creacadena.setProperty('qwidgetClass', 'creacadena')
        self.creacadena.setGeometry(QRect(50, 50, 1200, 450))
        self.creacadena.setWindowIcon(QIcon('images.ico'))
        self.creacadena.setWindowTitle('Crear Cadena')
        self.creacadena.show()
    #esta función abre la ventana para crear barcos
    def barcocrear(self):
        #####listas de la vista de los barcos (variables iniciales)####
        
        
        tablero = {'identificador':0,'tipo': 'waterjet','tiempo':0.0, 'alfab':0.0, 'wb':0.0,'theta':90,'Fm':0.0,'setl':0.0,'setw':0,
                   'interrv':0.0,'interrumbo':0.0,'thewj':0.0,'thewjmax':30,'Ac':10.0,'M':0.,'mb':350.,'mul2':0.,'mul':1000.0,
                   'mut':10000.0,'mut2':0.,'Ib':1000.0,'mua':1000.0,'mua2':0.,'ls':2.0,'pmax':5000.,'Tsf':0.1,'Tsw':0.1,'cfr':100,'link':0}
        krpid_tabla={'krpid[0]':1 ,'krpid[1]': 1  ,'krpid[2]': 0.1}	
        krd_tabla={'krd[0]':0.01 ,'krd[1]': 0.05  ,'krd[2]': 0.01}
        krl_tabla={'krl[0]':3 ,'krl[1]': 0.,'krl[2]': 0.}
        kvpid_tabla={'kvpid[0]':1000.0,'kvpid[1]': 500,'kvpid[2]': 2}
        ab_tabla={'ab[0]':0.0,'ab[1]':0.0}
        pb_tabla={'pb[0]':0.0,'pb[1]':0.0}
        vb_tabla={'vb[0]':0.0,'vb[1]':0.0}
        Fb_tabla={'Fb[0]':0.0,'Fb[1]':0.0}
        mA_tabla={'mA[0]':0.,'mA[1]':0.}
        interdist_tabla={'interdist[0]':0,'interdist[1]':0}
        #####termina la lista de las variables iniciales de los barcos#####
        self.es_Barco=False
        self.es_Cadena=False
        self.creabarcos = MyBarco(self.lista_barcos,None,tablero,krpid_tabla,krd_tabla,krl_tabla,kvpid_tabla,
		                          ab_tabla,pb_tabla,vb_tabla,Fb_tabla,mA_tabla,interdist_tabla)
        self.creabarcos.setGeometry(QRect(100, 100, 1250,550))
        self.creabarcos.setWindowIcon(QIcon('images.ico'))
        self.creabarcos.setWindowTitle('Crear Barcos')
        self.creabarcos.show()
	#función que abre la ventana para la edición de barcos y cadena
    def editarbarcoycadena(self):
        self.editabyc = MyEditaByC(self.lista_barcos,self.lista_cadenas)
        self.editabyc.setGeometry(QRect(100, 100, 400, 50))
        self.editabyc.setWindowIcon(QIcon('images.ico'))
        self.editabyc.setWindowTitle('Editar Barcos y Cadenas')
        self.editabyc.show()
    
    #función que abre la ventana para la simulación
    def simular(self):
        self.simula = MySimula(self.lista_barcoD_simu,self.lista_barcoI_simu,self.lista_cadenaC_simu,self.lista_simulacion,None)
        self.simula.setProperty('labelClass', 'simula')
        self.simula.setWindowIcon(QIcon('images.ico'))
        self.simula.setWindowTitle('Simular Barco con Cadena')
        self.simula.show()
        
class App(QApplication):
    def __init__(self, *args):
        QApplication.__init__(self, *args)
        self.main = MainWindow()
        self.aboutToQuit.connect(self.byebye )
        self.main.show()

    def byebye( self ):
        self.exit(0)
        
        
        
        

def main(args):
    global app
    app = App(args)
    app.exec_()

if __name__ == "__main__":
    main(sys.argv)

copy_reg.pickle(types.MethodType,MySimula.file_save,MySimula.cargar)
copy_reg.pickle(types.MethodType,MySimula.cargar_plot,MySimula.parar)
copy_reg.pickle(types.MethodType,MySimula.cargardatos,Worker.iniciar)
copy_reg.pickle(types.MethodType,Worker.procesar,Worker.abortar)
copy_reg.pickle(types.MethodType,usv.barco.controlador,usv.tim_prop)
copy_reg.pickle(types.MethodType,boom.cadena.movimiento,usv.barco.movimiento)
copy_reg.pickle(types.MethodType,usv.wj_prop,boom.smatrix.mtrl_s)
copy_reg.pickle(types.MethodType,boom.cadena.calcms,np.fix)
copy_reg.pickle(types.MethodType,sum,np.sqrt)
copy_reg.pickle(types.MethodType,np.abs,np.sum)
copy_reg.pickle(types.MethodType,np.nonzero,np.zeros)
copy_reg.pickle(types.MethodType,np.cos,np.sin)
copy_reg.pickle(types.MethodType,boom.smatrix.mtrr_s,boom.smatrix.mtr_s)
copy_reg.pickle(types.MethodType,boom.smatrix.mtrb_s,boom.smatrix.mtrsp_s)
copy_reg.pickle(types.MethodType,boom.smatrix.mtrlsp_s,boom.smatrix.mtrrsp_s)
copy_reg.pickle(types.MethodType,boom.smatrix.mtrbsp_s,np.sign)
copy_reg.pickle(types.MethodType,usv.barco.csbt,usv.barco.cstb)
copy_reg.pickle(types.MethodType,boom.smatrix.mtrsp_ini,boom.smatrix.mtr_ini)
copy_reg.pickle(types.MethodType,MyTable.tableItemChanged,MyTable.setmydata)
copy_reg.pickle(types.MethodType,MyTable.tableItemExplic)
copy_reg.pickle(types.MethodType,usv.barco.setAtributo,boom.cadena.setAtributo)
copy_reg.pickle(types.MethodType,usv.bezier_cvr.bezier4p,main)
copy_reg.pickle(types.MethodType,MyBarco.guardar,MyCadena.guardar_cadena)
copy_reg.pickle(types.MethodType,MyRadioLink.guardar,CheckableComboBox.handleItemPressed)
copy_reg.pickle(types.MethodType,MyCombinar.enlazar,MyCombinar.escojerbarcoderecho)
copy_reg.pickle(types.MethodType,MyCombinar.radiolink)
copy_reg.pickle(types.MethodType,MyCombinar.escojerbarcoizquierdo,MyCombinar.escojercadenacentro)
copy_reg.pickle(types.MethodType,MyCombinar.elongacion,MySimula.actualizar_etiqueta)
copy_reg.pickle(types.MethodType,MainWindow.barcocrear,MainWindow.cadenacrear)
copy_reg.pickle(types.MethodType,MainWindow.combinar,MainWindow.simular)
copy_reg.pickle(types.MethodType,MainWindow.editarbarcoycadena)

