# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 10:39:43 2016
This module contents the definition of the class cadena. Objets belonging to
this class are models of boom for oil spill recovery.
The module contains also different version of the function employed to
calculate the strains on the elements (links) of the boom. Each version is
apropiated for a different escenario. See below for more details
 
Descriptions for class and class procedures and functions
are included beneath its declaration
       
checked: 06.03.2017
reviewed: 11.05.2017 (mayor change, parameter Anch (anchored to dock included)) 
@author: juan

"""

import numpy as np
from matplotlib import pyplot as pl
from scipy.sparse import csr_matrix


import smatrix

class cadena:
    '''
    Constructor: 
    
    1. __init__(self,n=10,cable = np.zeros(3),Anch=[0,0]): 
        It builds an instance of the class. It can be called  without any
        parameter at all.
        n ->    is the number of links of the boom. If n is omitted during the 
                call, its default values is 10. The size of the computational
                 problem is related with the number of links of the boom, so 
                large booms(n>25) employs diferent procedures to calculate  
                the strains on the boon (see below)
        Cable ->  is an nparray with three values which represent the lenght in
                  meters of the towing cables attached to the boom tips and
                  their elasticity that it is suppose equal for both cables
                  [left, right,k]. If this parameter is omitted, it will take 
                  a 0 value for both tips. This means that the towing forces 
                  are directly applied to boom tips or/and that the towing 
                  ships are directly hitched to the boom. Please notice also 
                  that this will have and influence in the procedure use to 
                  calculate strains on the boom. (See below)
        Anch -> Contains values to anchor one of the tips of the boom to the 
                Dock: Anch =[1,0] means left tip in anchored Anch = [0,1]
                means right tip is anchored Anch = [1,1] means both tips are 
                anchores (this last situation can be use to represent a
                static boom and study the effect of streams, etc, but at 
                present this posibility is not implemented...)
            
        The remaining parameters of the object (boom) are defined by the
        constructor and take defaults values, that can be modified once the
        instance has been created. The rationale behind this is to avoid
        loading the call to the constructor with a lare number of parameters.
    
    Procedures and fuctions:
        
    2. resetear(self,primero = np.zeros(2),hard = False):
    3. listar():
    4. calcms(self,primero = np.zeros(2),ultimo = np.zeros(2), orden = 1):
    5. movimiento(self, delta = 5e-4):
    6. dibujar(self,color = 'black'):
    
    7. mtr_s Defines the method to calculate the strain matrix, needed to
       simulate the dynamic of the boom. The definition of this function varies
       according to the values contained in the variable cable, which are, in 
       turn related with the scenario the boom is deployed. See smatrix.py
       module for details, and comments belows, inside the constructor function
    '''    
    
    def __init__(self,esl=10,cable = np.zeros(3), Anch = np.zeros(2),objeto='cadena',s=2.0,q=0.1,s2=0,q2=0,
	             A=10.0,L=0.5,m=1.0,mA=0,Fd=np.zeros(2),Fi=np.zeros(2)):
        
        #self.objeto = 'cadena'
        
        #boom parameters definition:
        
        #numer of elements (links) building up the boom  
        self.objeto=objeto
        self.Anch=Anch
        self.cbl=cable		
        self.L=L		
        self.esl = esl
        #n=esl
        self.s=s
        self.q=q
        self.s2=s2
        self.q2=q2
        self.A=A
        self.m=m
        self.mA=mA
        self.Fd=Fd
        self.Fi=Fi
        #lenght of the towing cables cbl[0] = left,
        #cbl[1] = right, cbl[2] =  cables elasticity
        self.cbl = cable
        
        
        self.setEslabones()
    		
    def setEslabones(self):

        
        self.I = self.m * self.L**2 / 3 # m*(2*L)**2/12
        
        #unity vectors normal to each boom element
        #they are initilly generated pointing to +y
        self.normal = np.array([np.zeros(self.esl),np.ones(self.esl)])
        

        #unity vectors normal to each boom element. They point to the left end
        #of the element and are initilly generated pointing to -x
        self.para = np.array([-self.normal[1],self.normal[0]]) #vector paralelo 
        
        
        # mass center positions of the boom elements
        self.calcms()
        
        #speed of the boom elements. Repose reference system
        self.v = np.zeros(self.normal.shape)

        #speed of the boom elements. Reference to water
        self.vr =self.v
        
        #Speed norms. Inicialised to 1 to avoid division by zero problems during
        #the start up of the sistem 
        self.vmod = np.ones(self.esl)
        self.vmodr = self.vmod
        
        # aceleration  of the boom elements
        self.a = self.v.copy()        
        
        # boom elements, other kind of perturbations   (not implemented)
        self.per = self.v.copy()

        #angular aceleration boom elements 
        self.alfa = np.zeros(self.normal.shape[1])
        
        #angular speed of boom elements
        self.w = self.alfa.copy()
    
        #strain matrix 
        
        #The following function are necessary to calculate the strains in the
        #boom elements. there are alternative definition depending on the boom
        #situation. The constructor implement the simplest one for the simplest
        #scenarios. See module smatrix for details and to implement other
        #scenarios...
               
        if self.esl < 30:
            ################Matrix initialisation for short booms##############
            self.M = np.zeros((2 * self.esl + 2, 2 * self.esl + 2))
            self.M[0,0]= 1
            self.M[1,1]= 1
            self.M[-2,-2]= 1
            self.M[-1,-1]= 1
            #independent terms
            self.B = np.zeros(2 * self.esl + 2)
            #Strain applied to boom element ends 
            self.T = np.zeros(2 * self.esl + 2)
            
            #full matrix
            if sum(self.cbl) == 0:
                #no towing cables present (default sit)                
                self.mtr_s = smatrix.mtr_s.__get__(self,cadena)                
            elif self.cbl[0] == 0 and self.cbl[1] >= 0:
                #by default this situation represents that there are a towing
                #ship on the right of the boom and the left end is attached to
                #the dock. If this is the case, Fi[0] == -inf (the user should
                # set this at inizialization..) 
                #terms related with the strain in the first link are turned
                #to one to do them equal in both ends of the link
                if self.Anch[0] == 1:
                    self.M[0,2] = -1
                    self.M[1,3] = -1
                self.mtr_s = smatrix.mtrl_s.__get__(self,cadena)
            elif self.cbl[0] >= 0 and self.cbl[1] == 0:
                #by default this situation represents that there are a towing
                #ship on the left of the boom and the right end is attached to
                #the dock. If this is the case, Fd[0] == -inf (the user should
                # set this at inizialization..) 
                #terms related with the strain in the first link are turned
                #to one to do them equal in both ends of the link
                if self.Anch[1] == 1:
                    self.M[-2,-4] = 1
                    self.M[-1,-3] = 1
                self.mtr_s = smatrix.mtrr_s.__get__(self,cadena)                
            else:
                #if there are towing cables in both ends of the boom,
                #there are towing ships attached to them...
                #this includes the special situation of 'cables with length
                #0 cdl = [0,0,k] but elasticity k > 0 in this case the towing
                #ships are directly attached to the boom but the joint includes
                #is cosidered non-rigid...
                self.mtr_s = smatrix.mtrb_s.__get__(self,cadena)
        else:
            ################Matrix initialisation for long booms###############
            #an emtpy  sparse (csr) matrix is built. Dimension:  2*number of
            #links + 2. The csr (compressed sparse row format)  is the most
            #suitable format for the structure of the strains matrix when
            #working with spsolve... (see SciPy documentation for details )
   
            self.M = csr_matrix((2 * self.esl +2,2 * self.esl +2)) 
    
            #pointer to the non-null elements of M 
            self.M.indptr[1:3] = [4,8]
            self.M.indptr[3:-2] = np.arange(14,14+6*(self.M.indptr.size-5),6)
            self.M.indptr[-2:] = [self.M.indptr[-3] + 4, self.M.indptr[-3] + 8]
    
            #Indices of the non-null elements
            self.M.indices= np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5]\
            , dtype = 'int32')
            for i in range(4,self.M.indptr.size-4,2):
                self.M.indices=\
                np.concatenate((self.M.indices,i-2 + self.M.indices[8:20]))
            self.M.indices\
            = np.concatenate((self.M.indices, i + self.M.indices[0:8]))
            self.M.data = np.zeros(self.M.indices.size)
            self.M[0,0]= 1
            self.M[1,1]= 1
            self.M[-2,-2]= 1
            self.M[-1,-1]= 1
            self.B = np.zeros(2 * self.esl +2)
            self.T = np.zeros(2 * self.esl +2)
            
            #funtions for calculating the strains (sparse matrix)
            if sum(self.cbl) == 0:
                #no towing cable present (default -np.infsit)                
                self.mtr_s = smatrix.mtrsp_s.__get__(self,cadena)
            elif self.cbl[0] == 0 and self.cbl[1] >= 0:
                #by default this situation represents that there are a towing
                #ship on the right of the boom and the left end is attached to
                #the dock. Fi[0] == -infty.this should be initialised by the user
                if self.Anch[0] == 1:
                    self.M[0,2] = -1
                    self.M[1,3] = -1
                self.mtr_s = smatrix.mtrlsp_s.__get__(self,cadena)
            elif self.cbl[0] >= 0 and self.cbl[1] == 0:
                #by default this situation represents that there are a towing
                #ship on the left of the boom and the right end is attached to
                #the dock.Fi == -infty.this should be initialised by the user
                if self.Anch[1] == 1:
                    self.M[-2,-4] = 1
                    self.M[-1,-3] = 1
                self.mtr_s = smatrix.mtrrsp_s.__get__(self,cadena) 
            else:
                #if there are towing cables in both ends of the boom,
                #there are towing ships attached to them by default...
                self.mtr_s = smatrix.mtrbsp_s.__get__(self,cadena) 

 
    def setAtributo(self,atributo,valor):
          
        if atributo=="cable[0]":
            self.cbl[0]=float(valor)
        if atributo=="cable[1]":
            self.cbl[1]=float(valor)
        if atributo=="cable[2]":
            self.cbl[2]=float(valor)
        if atributo=="Anch[0]":
            self.Anch[0]=float(valor)
        if atributo=="Anch[1]":
            self.Anch[1]=float(valor)
        if atributo=="Fd[0]":
            self.Fd[0]=float(valor)
        if atributo=="Fd[1]":
            self.Fd[1]=float(valor)
        if atributo=="Fi[0]":
            self.Fi[0]=float(valor)
        if atributo=="Fi[0]":
            self.Fi[1]=float(valor)

    def resetear(self,primero = np.zeros(2),hard = False):
        
        '''2. resetear(self,primero = np.zeros(2)): 
        Reset the boom parameter to its initial values, ie, all parameter are 
        reset to its default values,except number of elements and cables lenght 
        that remain the same.
        primero ->  coordinates (x,y) for the center of mass of the first
                    (left-end) element; Default: (0,0).
                    The initial positions of the remainder elements are
                    calculated accordingly with the value of primero, and the
                    present orientation  (normal and parallel unit vectors )
                    of the elements, using function calcms(primero).  
        
        hard -> if hard is not false then, strains on the element-ends
                are set to zero and normal and parallel unit vector are also
                reset to its initial values: normal = [0,1] parallel [-1,0].
                The boom is then layed out parallel to the x axis.        
        
        WARNING, any change performed in the boom paramaters after created the
        boom objet are lost'''
    
        if  hard is not False:
            #Strain applied to boom element ends 
            if self.esl < 30:
                ################Matrix initialisation for short booms##########
                ##funtions for calculating the strains (strain matrix)
                self.M = np.zeros((2 * self.esl + 2, 2 * self.esl + 2))
                #independent terms
                self.B = np.zeros(2 * self.esl + 2)
                #Strain applied to boom element ends 
                self.T = np.zeros(2 * self.esl + 2)
            else:
                ################Matrix initialisation for long booms###############
                #an emtpy  sparse (csr) matrix is built. Dimension:  2*number of
                #links + 2. The csr (compressed sparse row format)  is the most
                #suitable format for the structure of the strains matrix when
                #working with spsolve... (see SciPy documentation for details )
                
                self.M = csr_matrix((2 * self.esl +2,2 * self.esl +2)) 

                #pointer to the non-null elements of M 
                self.M.indptr[1:3] = [4,8]
                self.M.indptr[3:-2] =\
                np.arange(14,14+6*(self.M.indptr.size-5),6)
                self.M.indptr[-2:]=\
                [self.M.indptr[-3] + 4, self.M.indptr[-3] + 8]      

                #Indices of the non-null elements
                self.M.indices=\
                np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5],\
                dtype = 'int32')
                for i in range(4,self.M.indptr.size-4,2):
                    self.M.indices=\
                    np.concatenate((self.M.indices,i-2 + self.M.indices[8:20]))
                    self.M.indices\
                    = np.concatenate((self.M.indices, i + self.M.indices[0:8]))
                    self.M.data = np.zeros(self.M.indices.size)
                    self.B = np.zeros(2 * self.esl +2)
                    self.T = np.zeros(2 * self.esl +2)
                
            #unity vectors normal to each boom element
            #they are initilly generated pointing to +y
            self.normal = np.array([np.zeros(self.esl),np.ones(self.esl)])        
            #unity vectors normal to each boom element. They point to the left
            #end of the element and are initilly generated pointing to -x
            self.para = np.array([-self.normal[1],self.normal[0]])
            
        #surface of a boom element
        self.s = 2.0        
        #thickness of a boom element
        self.q = 0.1 
        #damping parameter for quadratic damping factor (transverse) for a boom
        #element 
        self.s2 = 0
        #damping parameter for quadratic damping (longitudinal) for a boom
        #element
        self.q2 = 0       
        #damping paramater for turning for a boom element 
        self.A = 10.0        
        # boom element half-length
        self.L = 0.5        
        # boom element mass
        self.m = 1.0        
        # element added mass (only transverse)
        self.mA = 0        
        #Forces applied to boom tips. 
        self.Fd = np.zeros(2)
        self.Fi = np.zeros(2)        
        #moment of inertia for a boom element. 
        self.I = self.m * self.L**2 / 12        
        # mass center positions of the boom:   m elements
        self.calcms(primero)        
        #speed of the boom elements. Repose reference system
        self.v = np.zeros(self.normal.shape)
        #speed of the boom elements. Reference to water
        self.vr =self.v        
        #Speed norms. Inicialised to 1 to avoid division by zero problems during
        #the start up of the sistem  
        self.vmod = np.ones(self.esl)
        self.vmodr = self.vmod        
        # aceleration  of the boom elements
        self.a = self.v.copy()                
        # boom elements, other kind of perturbations   (not implemented)
        self.per = self.v.copy()        
        #angular aceleration boom elements 
        self.alfa = np.zeros(self.normal.shape[1])        
        #angular speed of boom elements
        self.w = self.alfa.copy()
        
                                         
        
    def listar(self):
        ''' Shows a brief list of the current values of the boom paramaters'''
        
        print "\n boom parameters"
        print "----------------------------------------------------------"
        print 'esl   :number of elements                       ', self.esl
        print 's     :element surface                          ', self.s
        print 's2    :element thick                            ', self.s2
        print 'q     :damping coeficient (transverse)          ', self.q
        print 'q2    :damping quadratic coeficient (transverse)', self.q2
        print 'A     :damping coeficient (turning)             ', self.A
        print 'Anch  :Tip(s) anchored to the dock              ', self.Anch 
        print 'L     :half lenght of an element                ', self.L
        print 'm     :element mass                             ', self.m
        print 'mA    :element added mass coef. (transverse)    ', self.mA
        print 'Fd    :force applied to the boom right end      ', self.Fd
        print 'Fi    :force applied to the boom left end       ', self.Fi
        print 'I     :moment of inertia                        ', '{0:.5f}'\
        .format(self.I)
        #print 'mtr_i :method to initialize the strain matrix   ',self.mtr_i
        print 'mtr_s :method to calculate the strain matrix    ',self.mtr_s
        print '----------------------------------------------------------\n'
        print ' Variables and their meaning (values are omited)'  
        print '-----------------------------------------------------------'      
        print 'normal:units vectors normal to the elements'
        print 'para  :units vector parallel to the element'        
        print 'cms   :elements center of mass'
        print 'v     :elements speed'
        print 'vr    :element speed ref. water'
        print 'vmod  :element speed norm'
        print 'a     :element acceleration '       
        print 'per   :element external perturbation'
        print 'alfa  :element angular acceleration'
        print 'w     :element angular speed'
        print 'T     :element tensions'
        print 'M     :strain matrix'
        print 'B     :Independet terms to calculate T, from M*T = B'
        print '----------------------------------------------------------\n'
    
    def calcms(self,primero = np.zeros(2),ultimo = np.zeros(2), orden = 1):
        '''This function calculates the center of mass of the elementes of
        boom, departing from the unit vectors normal to elements, and
        the center of mass of the left-end or the right-end element,
        primero -> center of mass of the left-end element
        
        ultimo ->  center of mass of the right-end element
        
        orden ->   orden =  1: the position of the left-end element (primero)
                   is taken as reference to locate the remainig elements,
                   otherwisethe right-end element (ultimo) is taken
                   as reference.
        '''
        self.cms = np.zeros(self.normal.shape)
        
        if orden == 1:
            self.cms[0,0] = primero[0]
            self.cms[1,0] = primero[1]
            for i in range(1,self.cms.shape[1]):
                #from the first one (left) to the last one (right)
                #centers of mass are calculate for any element,
                # departing from their normal vectors, the components  
                #aren the same ones exchanging x e y and  changing the sign of  
                #x in order to obtain a unit vector parallel to the element.
                #parallel direction are taken from left to right
                #then the position of hte previos elenment is added. 
                self.cms[0,i] = self.cms[0,i-1] + self.L * self.normal[1,i-1] \
                + self.L *self.normal[1,i]
                self.cms[1,i] = self.cms[1,i-1] - self.L * self.normal[0,i-1] \
                -self.L * self.normal[0,i]
            self.ord = 1

        else:
            self.cms[0,-1] = ultimo[0]
            self.cms[1,-1] = ultimo[1]
            for i in range(self.cms.shape[1]-2,-1,-1):                
                #from the last one (right) to the first one (left)
                #centers of mass are calculate for any element,
                # departing from their normal vectors, the components  
                #aren the same ones exchanging x e y and  changing the sign of  
                #x in order to obtain a unit vector parallel to the element.
                #parallel direction are taken from left to right
                #then the position of hte previos elenment is added. 
                self.cms[0,i] = self.cms[0,i+1] - self.L * self.normal[1,i+1] \
                - self.L * self.normal[1,i]
                self.cms[1,i] = self.cms[1,i+1] + self.L * self.normal[0,i+1] \
                + self.L * self.normal[0,i]
                
            self.ord = -1
            
    def movimiento(self, delta = 5e-4):    
        '''This function calculates linear and angular acceleration for the
        elements of the boom, departing for the tensions and speeds. These last
        are used to obtain the damping forces.
        Once the accelaration have been estimated, speeds and position and
        orientation are calculated one-step-ahead using a simpe Euler
        esquema.
        delta -> step size for one-step-ahead Euler calculation. 
        '''
        
        #diference among the strains at the ends of the element        
        dT = np.array([self.T[2:-1:2] - self.T[:-3:2],\
        self.T[3::2] - self.T[1:-2:2]])
        sT = np.array([self.T[2:-1:2] + self.T[:-3:2],\
        self.T[3::2] + self.T[1:-2:2]])
        #relative speeds to water are used to take into account  
        #the currents drag
        self.vmodr = np.sqrt(sum(self.vr**2,0))
        #zero speed norms are change to one to avoid overflow errors            
        self.vmodr[np.nonzero(self.vmodr == 0)] = 1
         
        #Then linear accelerations are calculated            
        self.a = (dT - ((np.abs(np.sum(self.vr*self.normal,0)) * self.s\
        + np.abs(np.sum(self.vr*self.para,0)) * self.q) /self.vmodr\
        + np.abs(np.sum(self.vr*self.normal,0)) * self.s2\
        + np.abs(np.sum(self.vr*self.para,0)) * self.q2)
        * self.vr)/ (self.m + self.mA * self.normal)
        
        
        #and also angular acelerations
        self.alfa =(np.sum(sT * self.normal,0) * self.L - self.A * self.w)\
        /self.I             
        
        #speed are recalculated one-step-ahead using an Euler schema
        
        
        #Linear speed of CoM
        self.v += self.a * delta
               
        #and then  angular speed around the CoM...
        self.w += self.alfa * delta 
        
        #Orientation for the element are modifify accordingly
        
        cthsl = self.normal * np.cos(self.w * delta)
        sthsl = self.normal * np.sin(self.w * delta)
        
        self.normal[:] = [cthsl[0,:] - sthsl[1,:],
        cthsl[1,:] + sthsl[0,:]]
        self.para[:] = [-self.normal[1,:], self.normal[0,:]]
        
        #and also the CoM location
        self.cms += self.v * delta 
            
    def dibujar(self,color = 'k'):
        ''' This function draw schematicly the boom in its current pose\n
        ---o---                      ---o--- \n        
              \\                    /        \n 
                o                  o         \n 
                \\                /          \n
                   ---o--- ---o---           \n 
            color -> a color for the edges, black by default  
        '''
           
        #pl.plot(self.cms[0,:],self.cms[1,:],'o')
        #pl.hold(True)
        barrasi = self.cms + self.L * self.para 
        barrasd = self.cms - self.L * self.para
        return barrasi,barrasd
        #pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],color)
        
    
    
        
            
            