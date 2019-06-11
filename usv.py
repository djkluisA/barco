# -*- coding: utf-8 -*-
"""
NUEVO INTENTO DE CONTROL BASADO EN ANGULO Y DISTANCIA
Created on Tue Jul 12 09:14:20 2016
This module contains the class barco, which can generate instances of USV
 models. 
The class includes several methods to simulated the dynamics of an USV and also
to control its movement. The propulsion system is selected when an instance 
is created. There are at present two posibilities: waterjet and propeller. The
propulsion system afects the dynamics of the USV and also its control.
The USV model can be used for boom towing.
Description for the class contructor and methods are defined below their
declaration.
 
class barco

external functions:
1. trazar(barco,delta,rumbo,vel,tf,paso,ffv=0):
2. wj_prop(barco,delta= 5e-4,extf = [0,0],dfcm = [0,0]):
3. tim_prop(barco,delta= 5e-4,extf = [0,0],dfcm = [0,0] ):

@author: juan
"""



import numpy as np
from numpy.linalg.linalg import norm 
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
import bezier_cvr

#import time

class barco:
    '''
    Constructor:
    
    1. __init__(self,ide,tipo = 'waterjet'):
    It builts an instance of the class.
    ide ->  is an identification number. It allows to distiguish different
            instances of the class when they have been created in the same
            simulation scenario. This parameter is mandatory.
    tipo -> defines the propulsion system for the usv. If not defined, the usv
            will be funrnished with an swivel waterject propulsion system.
            Currently, there is only another propulsion system implemented
            which employ a normal propeller for propulsion and a rudder to
            modify the usv heading. To implement this system, tipo has to be
            set to 'propeller'
            See functions wj_prop and tim_prop, at the end of this module, for 
            details.
    The remaining parameters of the usv are defined by the constructor and can
    be modified once the usv instance has been created.
    
            
    Procedures and functions:

    2. cstb(self,rt):
    3. csbt(self,rb):
    4. dibujar(self):
    5. movimiento(self,delta = 5e-4):
    6. controlador(self, delta, c_rumbo=0, c_velocidad=0,
       ffv = 0,dcons = np.zeros(2)):
    7. planificador(self, puntos, rumbos, veloc = None, tiempos = None,
       ndat =10.  ):
    8. resetear(self):
    9. listar(self):
   10. trazar(barco,delta,rumbo,vel,tf,paso,ffv=0):

    '''
            
    #def __init__(self,ide,tipo = 'waterjet'):
    def __init__(self,ide,tipo='waterjet',tiempo=0.0,pb=np.array([0.0,0.0]),vb= np.array([0.0,0.0]),ab=np.array([0.0,0.0]),alfab=0.0,wb=0.0,theta=np.pi/2.0 ,Fm=0.0,Fb=np.array([0.0,0.0]),setl=0.0,setw=0.0,interrv=0.0,interrumbo=0.0,thewj=0.0,thewjmax=np.pi/6.0,Ac=10.0,M=0.0,mb=1000.0,mA=np.array([0.,0.]),mul=1000.0 ,mul2=0.,mut=10000.0,mut2=0.,Ib=1000.0,mua=1000.0,mua2=0.,ls=2.0,pmax=50.0*735.0,Tsf=0.1,Tsw=0.1,cfr=100,interdist=np.zeros(2) ,link=0,kvpid=[1000.0, 500, 2],krpid=[1, 1, 0.1],krd=[0.01,0.05,0.01],krl=[3,0.,0.]):
         
       #USV parameters difinition
       
       #USV identification number
       self.ide = ide
       
       #USV type, (waterjet or propeller) 
       self.tipo = tipo
       
       #total time (simulation elapsed time)
       self.tiempo = tiempo
       
       #position x,y for the USV center of mass (CoM)
       self.pb = pb
       
       #speed vx, vy
       self.vb = vb
       
       # acceleration ax, ay
       self.ab = ab
       
       #relative to water speed vrx, vry
       self.vr = vb 
              
       #usv angular acceleration
       self.alfab = alfab
       
       #usv angular speed
       self.wb = wb 
       
       #usv heading (x axis is taken as the reference: 0 heading)
       self.theta = theta 
       
       #usv propulsion force  
       self.Fm = Fm;
       
       #components of the total force applied to the USV
       self.Fb = Fb
       
       #set point for the USV engine
       self.setl = setl
       
       #set point for the rudder or waterjet swivel
       self.setw = setw
       
       #speed error (integral)
       self.interrv = interrv #integral del error de velocidad
       
       #heading error (integral)
       self.interrumbo = interrumbo
	   #integral del error de rumbo
       
       #actual waterject or rudder pose
       self.thewj = thewj
       
       #maximum turning angle for waterject or rudder
       self.thewjmax = thewjmax
       
       #Distance from the rudder or waterject to the usv mass center
       self.Ac = Ac

       #torque applied to the ship by the propulsion force  
       self.M = M
       
       #USV mass 
       self.mb = mb #masa del barco en kg
       
       #USV added masses on surge and sway direcctions
       self.mA = mA
       
       #damping parameter for linear damping force in surge direction
       self.mul = mul 

       #damping parameter for quadratic damping force in surge direction
       self.mul2 = mul2

       #damping parameter for linear damping force in sway direction
       self.mut = mut 
       
       #damping parameter for quadratic damping force in sway direction
       self.mut2 = mut2        

       #usv moment of inertia around the usv mass center
       self.Ib = Ib

       #angular damping parameter yaw  
       self.mua = mua
       
       #cuadratic damping parameter yaw
       
       self.mua2 = mua2
       #usv lenght   mb,mA,mul,mul2,mut,mut2,Ib,mua,mua2,ls,pmax,Tsf,Tsw,Ab,Aw,cfr,interdist,link,rlink,
       self.ls = ls
       
       #usv engine max. power
       self.pmax = pmax #Fuerza maxima motor
       
       #engine time constant
       #engine behaviour is model as a first order system         
       self.Tsf = Tsf
       
       #water-jet swivel or rudder time constant 
       #modelled as a first order system
       self.Tsw = Tsw 

       #I do not know what the hell are the two following parameters.... 
       self.Ab = 0.1 #
       self.Aw = 0.1 #
       
       
       #rudder shape factor (not applied to usv with waterjet system)
       self.cfr = cfr 
       

       #####following parameters apply only for cooperative towing task########              #
       
       #distance among usv-mates error(integral)
       #(only applied on cooperative boom towing)
       self.interdist = interdist      
       
       #element of the boom the usv is attached to. This parameter is only
       #applied to towing scenaries with more than two usv towing the same boom
       #and only to those usv that are attached to an element different from
       #first one or the last one.
       self.link = link
       
       
       #radio link of the usv when implemented allows to recibe external set 
       #points, data from other usv and  so for. Currently this parameter is 
       #only partially defined....
       #self.rlink = rlink 
       

       ######################propulsion systems################################
        
       #propulsion system and PID controllers default parameters are asigned
       #according to the usv type
        
       if self.tipo == 'waterjet':
	      self.propulsion =  wj_prop.__get__(self,barco)
       else:
          self.propulsion = tim_prop.__get__(self,barco)
           
           #constans for PID controllers (order is always: kp, kd, ki)  
           #constants for the speed controller 
       self.kvpid = kvpid 
          
           #constant for the speed controler related to the distance errors
           #between USVs (only applied to cooperative towing tasks)
          
           
           #constants for heading controller
       self.krpid = krpid 
           
           
        
       self.krd =krd #distance 
           
       self.krl =krl #alignment
       
           #constans for PID controllers (order is always: kp, kd, ki)
           #constants for the speed controller 
       #self.kvpid = [1000.0, 500, 2] 
          
           
           
           
           #constants for heading controller
        #self.krpid = [1, 1, 0.1] 
           
           #constants for heading controler related to the distance errors
           #between USVs (only applied to cooperative towing tasks)
        
        #self.krd =[0.01,0.05,0.01] #distance 
           
        #self.krl =[3,0.,0.] #alignment
               
    def cstb(self,rt):
        '''cstb(self,rt):
        This function transform a vector rt from fixed (earth)
        reference system to
        to usv reference system.
        The usv reference system x+ axis is oriented on usv surge
        direction. y+ on sway (x+ +90 ) 
        rt -> vector coordinates in fixed reference system        
        '''
        tc = np.cos(self.theta)
        ts = np.sin(self.theta)
        rb = np.dot(np.array([[tc,ts],[-ts,tc]]),rt)
        return(rb)
		
		
    def getAtributo(self,atributo):
        
		
        if atributo=="Tsw":
           return self.Tsw
        if atributo=="cfr":
           return self.cfr
        if atributo=="Tsf":
           return self.Tsf
        if atributo=="pmax":
           return self.pmax
        if atributo=="ls":
            return self.ls
        if atributo=="Ib":
            return self.Ib
        if atributo=="mua":
            return self.mua
        if atributo=="mua2":
            return self.mua2
        if atributo=="mut":
            return self.mut
        if atributo=="mut2":
            return self.mut2
        if atributo=="mul":
            return self.mul
        if atributo=="mul2":
            return self.mul2
        if atributo=="thewjmax":
            return self.thewjmax
        if atributo=="mb":
            return self.mb
        if atributo=="M":
            return self.M
        if atributo=="Ac":
            return self.Ac
        if atributo=="interrv":
            return self.interrv
        if atributo=="interrumbo":
            return self.interrumbo
        if atributo=="thewj":
            return self.thewj
        if atributo=="ide":
            return self.ide
        if atributo=="tipo":
            return self.tipo
        if atributo=="tiempo":
            return self.tiempo
        if atributo=="alfab":
            return self.alfab
        if atributo=="wb":
            return self.wb
        if atributo=="theta":
            return self.theta
        if atributo=="Fm":
            return self.Fm
        if atributo=="setl":
            return self.setl
        if atributo=="setw":
            return self.setw
        if atributo=="ab":
            return self.ab
        if atributo=="pb":
            return self.pb
        if atributo=="vb":
            return self.vb
        if atributo=="Fb":
            return self.Fb
        if atributo=="mA":
            return self.mA
        if atributo=="interdist":
            return self.interdist
        if atributo=="kvpid":
            return self.kvpid
        if atributo=="krpid":
            return self.krpid
        if atributo=="krd":
            return self.krd
        if atributo=="krl":
            return self.krl

    def setAtributo(self,atributo,valor):
        if atributo=="ab[0]":
            self.ab[0]=float(valor)
        if atributo=="ab[1]":
            self.ab[1]=float(valor)
        if atributo=="pb[0]":
            self.pb[0]=float(valor)
        if atributo=="pb[1]":
            self.pb[1]=float(valor)
        if atributo=="vb[0]":
            self.vb[0]=float(valor)
        if atributo=="vb[1]":
            self.vb[1]=float(valor)
        if atributo=="Fb[0]":
            self.Fb[0]=float(valor)
        if atributo=="Fb[1]":
            self.Fb[1]=float(valor)
        if atributo=="mA[0]":
            self.mA[0]=float(valor)
        if atributo=="mA[1]":
            self.mA[1]=float(valor)
        if atributo=="interdist[0]":
            self.interdist[0]=float(valor)
        if atributo=="interdist[1]":
            self.interdist[1]=float(valor)
        if atributo=="kvpid[0]":
            self.kvpid[0]=float(valor)
        if atributo=="kvpid[1]":
            self.kvpid[1]=float(valor)
        if atributo=="kvpid[2]":
            self.kvpid[2]=float(valor)
        if atributo=="krpid[0]":
            self.krpid[0]=float(valor)
        if atributo=="krpid[1]":
            self.krpid[1]=float(valor)
        if atributo=="krpid[2]":
            self.krpid[2]=float(valor)
        if atributo=="krd[0]":
            self.krd[0]=float(valor)
        if atributo=="krd[1]":
            self.krd[1]=float(valor)
        if atributo=="krd[2]":
            self.krd[2]=float(valor)
        if atributo=="krl[0]":
            self.krl[0]=float(valor)
        if atributo=="krl[1]":
            self.krl[1]=float(valor)
        if atributo=="krl[2]":
            self.krl[2]=float(valor)
        
    def csbt(self,rb):
        '''csbt(self,rb):
        This function transform a vector rb from usv
        reference system to
        to fixed (earth) reference system.
        The usv reference system x+ axis is oriented on usv surge
        direction. y+ on sway (x+ +90 ) 
        rb -> vector coordinates in usv reference system
        '''
        tc = np.cos(self.theta)
        ts = np.sin(self.theta)
        rt = np.dot(np.array([[tc,-ts],[ts,tc]]),rb)
        return(rt)
        
    def dibujar(self,theta,pb0,pb1,color = 'blue'):
        '''dibujar(self):
        draw the usv at its current pose'''
        #pl.plot(self.pb[0],self.pb[1],'+')
        vertices = np.array([[-self.ls/2.,-0.25*self.ls/2.],[-self.ls/2.,0.25*self.ls/2.],\
        [-0.25*self.ls/2.,0.35*self.ls/2.],\
        [self.ls/2.,0],[-0.25*self.ls/2.,-0.35*self.ls/2.],\
        [-self.ls/2.,-0.25*self.ls/2.]])
        rot = np.array([[np.cos(theta),- np.sin(theta)],[np.sin(theta), np.cos(theta)]])  
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [pb0,pb1]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3]
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = color)
        return patchi
        #pl.gca().add_patch(patchi)
        

    def movimiento(self,delta = 5e-4):
        '''movimiento(self,delta = 5e-4):
        Calculates the usv movement. It performs an integration step
        for acceleration and speed, following a simple Euler's method.
        It uses a crossed method: speed are integrate from acceleration and
        THEN the speed so obtained is integrated in turn to obtain the
        usv position,
        delta -> contains the integration step(time) used to calculate the usv
                 movement
        '''
        
        #usv speed
        self.vb[0] += self.ab[0] * delta
        self.vb[1] += self.ab[1] * delta 
        
        #usv position
        self.pb[0] += self.vb[0] * delta
        self.pb[1] += self.vb[1] * delta
        
        #usv angular speed
        self.wb += self.alfab * delta 
        
        #heading
        self.theta += self.wb * delta
        
        #The heading is got into the interval -pi <angle < pi 
        self.theta = self.theta -np.fix(self.theta/2./np.pi)\
        -2 * np.pi *((self.theta -np.fix(self.theta/2./np.pi)) > np.pi)\
        +2. * np.pi * ((self.theta -np.fix(self.theta/2./np.pi)) <= -np.pi)
        
        #total simulation time is updated
        self.tiempo = self.tiempo + delta
                
        if self.rlink:
            #if a radio link is available, then position, heading and speed are
            #written in the rlink buffer... 
            self.rlink[0][1] = self.pb
            self.rlink[0][2] = self.theta
            self.rlink[0][3] = self.vb
            
            
    def controlador(self, delta, c_rumbo=0, c_velocidad=0,\
    ffv = 0,dcons = np.zeros(2)):
        '''
        This function implements the speed and heading controller for the 
        usv
        There are several terms corresponding to the expecific scenario the 
        usv is included into. For a single usv only heading and speed error
        are relevant for the case of cooperative towing, them errors in the
        distances between the usv and its usv-mate are also taken into
        account.
        
        delta -> simulation step time

        c_rumbo -> heading set point. For two usv towing a boom, this set
                   point should be common to both of them.

        c_velocidad -> usv speed set point. For two usv towing a boom, this set
                       point should be common to both of them.
                       
        ffv -> feedforward term added to the speed controller. 
               It should be based on trial performed with the usv model
               once its parameter has been fit.

        dcons -> this terms stablishes a set point for the steady navigation
                 distance between two USVs towing a boom.
                 dcons[0] represents the set point navigation distance in surge
                 direction, suposing both usv are navigating parallel along
                 the heading set point direction. for the left-end usv this 
                 distance should be negative and for the right-end ship should
                 be positive. 
                 dcons[1] represents the set point distance in sway direction, 
                 supossing both USVs are navigating parallel along the heading
                 set point direcction. This value is taken positive for the
                 usv in leading position and negative for the ship in following
                 position. Zero in case both usv navigate aligned
            
                   
        
        '''
        #heading vector 
        rumv = np.array([np.cos(self.theta),np.sin(self.theta)])
        
        #speed error         
        errorv = c_velocidad - np.dot(rumv,self.vb)
        
         #integral error with antiwindup (speed) 
        self.interrv += errorv * delta \
        * (np.abs(self.setl) < np.abs(self.pmax))
          
        #heading error
        errumbo = c_rumbo - self.theta\
        -2. * np.pi * ((c_rumbo - self.theta) > np.pi)\
        +2. * np.pi * ((c_rumbo - self.theta) <= -np.pi)
        
      
         
        #the following part of the controller is only implemented if there is a 
        # a radio link active and coordination between two USVs. 
        if self.rlink:
            #heading set point transfor matrix
            c_rumv = np.array([[np.cos(c_rumbo),np.sin(c_rumbo)],\
            [-np.sin(c_rumbo),np.cos(c_rumbo)]])
            #set point distance between USVs             
            c_dist = norm(dcons)
            #set point angle between ships Line-of-Sight and c_rumbo 
            c_par = np.arctan2(dcons[1],dcons[0])
            #distance errors
            
            #vector from the usv to its mate...
            vpi = np.array(self.rlink[1][1]- self.rlink[0][1])
            #vpi represented in c_rumv ref system
            vpicr = np.dot(c_rumv,vpi)
            par  = np.arctan2(vpicr[1],vpicr[0])
            #dist between USVs
            dist = norm(vpi)
            den = dist + (norm(vpi) == 0.)
            #matrix to transform vector for earth reference system to 
            #a system with the x axis laying in the direcction of vpi..
            mvpi = np.array([[vpi[0]/den,vpi[1]/den],\
            [-vpi[1]/den,vpi[0]/den]])
            
            #
             
            
            #components of USV heading vector in vpi ref. system
            rel = np.dot(mvpi,rumv)
            
            
            erdist = np.array([(c_dist -dist)*rel[1],(c_par - par)*rel[1]])
            derdist = - np.dot(mvpi,\
            np.array(self.rlink[1][3]  - self.rlink[0][3])*rel[1])
            self.interdist = self.interdist + erdist * delta
#            print self.ide, derdist
            #dists contains the surge and sway current distances between the
            #USVs, referred to the heading set point...
#            mc_rumbo = np.array([[np.cos(c_rumbo),np.sin(c_rumbo)],\
#            [-np.sin(c_rumbo),np.cos(c_rumbo)]])
            
#            dists = np.dot(mc_rumbo,vpi)           
            
            #distances error normalised to dcons
#            erdist = (dcons - dists)/(np.abs(dcons) + (dcons == 0)) 
#            print self.ide, c_par ,par, erdist
            
            #distances error derivatives
#            derdist = - np.dot(mc_rumbo,\
#            np.array(self.rlink[1][3]  - self.rlink[0][3]))
            
            #if both USVs are initially far apart eah other, then it has no
            #to implement a control system... First the USVs manouver to
            #approach each other as fast as possible... cdn controls this 
            #situation the conditions to start up the controller is that the 
            #distance between the USVs be less than 2 times the sway set point
            #distance
             
#            cdn = np.abs(erdist) <= 2 
            
            #integral distances error. cdn is acting as antiwind-up when 
            #the distance between the USVs is too large. 
#            self.interdist = self.interdist + erdist * delta * cdn

            #errumbo is recalculated here:
            #the first term switches on when USVs are near enought each other
             #(see cdn above).
            #the second term switches on when they are far appart and
            #it tries to get them near each other (errumbo is recalculted
            #considering the direction from the USv to its mate ). 
           
            
#            errumbo = errumbo if np.abs(errumbo) < np.abs( phi - self.theta\
#            -2. * np.pi * ((phi - self.theta) > np.pi)\
#            +2. * np.pi * ((phi - self.theta) <= -np.pi)) else\
#             phi - self.theta -2. * np.pi * ((phi - self.theta) > np.pi)\
#            +2. * np.pi * ((phi - self.theta) <= -np.pi)                    
            
           
                        
            
        else:
            #revisar: posiblemente esto es innecesario si las constantes
            #del pid correspondientes a estos terminos son cero... (ajustar por
            # defecto )    
            erdist = np.zeros(2)
            derdist = np.zeros(2)
           
            
        #integral error with antiwindup 
        self.interrumbo = self.interrumbo + errumbo * delta \
        * (np.abs(self.setw) < np.abs(self.thewjmax))
            
        #The heading controller, consider the heading error and the distance
        # between the USVs for the caseof two USVs towing a boom. NOTICE: 
        #the minus sign in the controler imputs is necesary because the rudder
        #and the usv turn in opposit directions
        self.setw = -self.krpid[0] * errumbo + \
        self.krpid[1] * self.wb - \
        self.krpid[2] * self.interrumbo - \
        (self.krd[0] * erdist[0] + \
        - self.krd[1] * derdist[0] + \
        self.krd[2] * self.interdist[0]) 
        
        #Speed Controller The spped controller, has a feedforward term,
        # regulates the speed according to the speed set point and also 
        #the distance between USVs in surge direcction, to achive a final
        #configuration in 'U' or 'J'...
        #The terms to control the distance between both USVs are multiplied
        #by the escalar product of the actual heading and the set point heading
        #this seasons the efford of control specialy y the USVs are not 
        #alined in the desired heading.
        
        self.setl = ffv + self.kvpid[0] * errorv  + \
        self.kvpid[1] * (np.cos(self.theta) * self.ab[0] \
        + np.sin(self.theta)) * self.ab[1] + \
        self.kvpid[2] * self.interrv -( \
        +self.krl[0] * erdist[1] +\
        -self.krl[1] * derdist[1] +\
        +self.krl[2] * self.interdist[1])
        
          
        
    def planificador(self, puntos, rumbos, veloc = None, tiempos = None,
    ndat =10.  ):
        ''' This planning function build links, using Bezier's 3th degree
        polynomials, a set of way-points to be gone over by an USV. 
        the resulting curves describe a trayectory to be followed by an USV.
        it uses some funtions from bezier_cvr module. (see bezier_cvr for 
        details)

        puntos -> array containing the setpoints to be gone over by the USV.
        
        rumbos -> array of USV heading for the set points include in puntos.
        
        veloc  -> when supplied, array of USV speed for the set points included
                  in puntos. If omitted, speeds on the set points are taken
                  equal to 0.5 m/s, supposing also that the USV starts from
                  repose
                  
        tiempos -> when supplied, array of USV arriving time to the set points
                   included in puntos.               

        

        CAVEAT: This plannig function is still under construction. 
                At present, it does not take into account the real dynamic of
                the USV, thus it can supply unfeasible trajectories.
        
                
        '''
        
        if veloc is None:
             #When no speeds for the set points are supplied,
             #they are taken eqult to 0.5 m/s. The initial USV speed is taken
             #as cero.
             #Anyway, these default values could be change here,
             veloc = np.repeat([0.,0.5],[1,puntos.shape[0]-1])
        if tiempos is None:
            #When no set points arriving times are supplied, they are obtained
            #from the mean between to sucesive set points speeds 
            #and the distance beteen them.
                        
            tiempos = 2 * norm(puntos[1:] - puntos[:-1],axis =1).squeeze() / \
            (veloc[1:] + veloc[:-1])
            
            print tiempos
            
        lstplan = []
            
        for i,j,k,l,m,n,o in zip(puntos[:-1],puntos[1:],veloc[:-1],
        veloc[1:],tiempos, rumbos[:-1], rumbos[1:]):
            #these prints are left here on purpose to let inspect the results 
            print 'estoy en el bucle\n'
            print i, '\n'
            print j, '\n'
            print k, '\n'
            print l, '\n'
            print m, '\n'
            print n, '\n'
            print o, '\n'
            print ndat
            
            lstplan.append(bezier_cvr.bezier4p(i,j,k,l,m,n,o,ndat))
        
        return lstplan
            

    def resetear(self):
        '''  resetear(self:)
        Turns the USV variables to their default values '''
        self.pb = np.array([0.0,0.0]) #position x,y
        self.vb = np.array([0.0,0.0]) #speed, vy
        self.vr = np.array([0.0,0.0]) #speed relative to water vrx, vry
        self.ab = np.array([0.0,0.0]) # aceleration ax, ay
        self.alfab = 0.0 #angular aceleration
        self.wb = 0.0 #angular speed
        self.theta = np.pi/2.0 #USV heading (zero is x axis direction)
        self.Fm = 0.0; #force applied to the ship
        self.setl = 0.0 #USV engine setpoint
        self.setw = 0.0 #USV rudder/waterjet orientation setpoint
        self.thewj = 0.0 #rudder/waterjet orientation
        self.thewjmax = np.pi/6.0 # rudder/waterjet  maximum turn angle
        self.tiempo = 0.0 #simulation time
        
    def listar(self):
        ''' listar(self):
        shows a brief list of the current values of the USV parameters
        '''
        print "\n USV parameters:"
        print "-------------------------------------------------------------"
        print 'ide   :identification number                  ', self.ide
        print 'thewjmax :Maximum watertjet/rudder tunr angle ', self.thewjmax
        print 'Ac       :Distance from the rudder axis to CoM', self.Ac                
        print 'mb       :USV mass                            ', self.mb
        print 'mA       :added mass                          ',self.mA
        print 'mul      :damping linear coeficient  (surge)  ', self.mul
        print 'mul2     :damping cuadratic coeficient (surge)',self.mul2
        print 'Ib       :Inertia moment                      ', self.Ib 
        print 'mua      :damping coeficient (yaw)            ', self.mua
        print 'mua2     :damping cuadratic coeficient (yaw)  ', self.mua2
        print 'ls       :lenght                              ', self.ls
        print 'mut      :damping linear coeficient (sway)    ', self.mut
        print 'mut2     :damping cuadratic coeficient (sway) ', self.mul2
        print 'pmax     :maximun engine power                ', self.pmax        
        print 'Ab       :I do not know what this paramter is ', self.Ab 
        print 'Aw       :I do not know what this param is    ', self.Aw
        print 'Tsf      : time constant of the USV engine    ', self.Tsf
        print 'Tsw      : time constant of USV guiding system', self.Tsw
        print 'tipo     : type of USV propulsion system      ', self.tipo
        print 'cfr      : form factor for rudder guided USVs ', self.cfr 
        print '------------------------------------------------------------\n'
        print 'variables and their meaning (values are omited)'
        print '------------------------------------------------------------'
        print 'pb       :position'                             
        print 'vb       :speed'                               
        print 'vr       :speed relative to water'             
        print 'ab       :aceleration                         '
        print 'alfab    :angular aceleration                 '
        print 'wb       :angular speed                       ' 
        print 'theta    :heading                             '
        print 'Fm       :exerted force                       '
        print 'setl     :engine setpoint                     ' 
        print 'setw     :heading setpoint                    '
        print 'thewj    :waterjet/rudder orientation         '
        print 'M        :moment applied to the USV           '
        print 'tiempo   :USV simulation time                 '
        print '------------------------------------------------------------\n'

def trazar(barco,delta,rumbo,vel,tf,paso,ffv=0):
    '''trazar(barco,delta,rumbo,vel,tf,paso,ffv=0): 
    This function allows to draw the trajectory described by an USV.
    
    barco -> an instance of barco class
    
    delta -> time step for integration
    
    tf -> final simulation time
    
    paso ->  interval for drawing (paso should be greater or equal to delta) 
    ffv -> feedforward parameter for speed controller
    '''
    
    #array of time to calculate the integration results
    tp = barco.tiempo   
    t = np.arange(tp,tf,delta)
    
    #figures to draw the results
    #trajectory
    pl.figure(1)
    pl.hold(True)
    pl.axis('equal')
    pl.xlabel('x')
    pl.ylabel('y')

    #power
    pl.figure(2)
    pl.hold(True)
    pl.xlabel('t (s)')
    pl.ylabel('Power (Kw -> 0.735 CV)')


    #speed
    pl.figure(3)
    pl.hold(True)
    pl.xlabel('t (s)')
    pl.ylabel('speed (m/s)')
    
    #comtrol?
    pl.figure(4)
    pl.hold(True)
    pl.xlabel('t')
    pl.ylabel('setl/1000(s)y setw/0.5(t)')
    
    thp = barco.theta
    vp = barco.vb
    bslp = barco.setl
    bsvp = barco.setw
    pwp = 0
    for i in range(np.size(t)):
        #run the model...
        barco.controlador(delta,rumbo,vel,ffv)
        barco.propulsion(delta)
        barco.movimiento(delta = delta)
        
        if i%paso == 0:
            #draw the results
            #usv trajectory 
            pl.figure(1)            
            barco.dibujar()
            
            #?
            pl.figure(2)
            pl.plot([tp,t[i]],[pwp,norm(barco.vb)*barco.Fm/1000.],'b')
            
            #usv speed
            pl.figure(3)
            pl.plot([tp,t[i]],[np.cos(thp) * vp[0] +
            np.sin(thp) * vp[1],np.cos(barco.theta) * barco.vb[0] +
            np.sin(barco.theta) * barco.vb[1]],'b')
            
            
            pl.figure(4)
            pl.plot([tp,t[i]],[bslp/10000,barco.setl/10000],'k')
            pl.plot([tp,t[i]],[bsvp/0.5,barco.setw/0.5],'g')
                       
                     
            #advance the varible one step ahead
            tp = t[i]    
            thp = barco.theta
            vp = barco.vb.copy()
            bslp = barco.setl
            bsvp = barco.setw 
            pwp = norm(barco.vb)*barco.Fm/1000.


def wj_prop(barco,delta= 5e-4,extf = [0,0],dfcm = [0,0]):
        '''
        This function models a waterjet or other outboard propulsion system.
        Its main feature is the fully accopling between propulsion and
        guiding systems. Departing from the forces a torques exerted on the
        USV, the function turno out the USV acceleration. 
        Whenever a instance of the class barco is created with the option 
        tipo = 'waterjet', the parameter propulsion of the instance points to 
        this function. Making the USV to adopt a waterjet-type propulsion
        and guiding system.
        
        barco -> an instance of the class barco. (Rarely use because all 
        instance of barco has already its own propulsion system defined at
        construction time. Thus the function is mainly use as 
        barco_instance.propulsion(reamining parameters))
        
        delta -> step time for integration
        extf -> external force, (usually the last boom element strain)
        dfcm -> aplication point of the external force... notice that only
                one external force can be applied
        '''
        
        
        #First, the function takes the current engine setpoint value and
        #and assures the setpoint value is within its limits...
        barco.setl = ((barco.setl >= 0) & (barco.setl <= barco.pmax))\
        * barco.setl + barco.pmax * (barco.setl >barco.pmax )        
        
        #the engine is modelled as a first order system. Tsf is the model
        #raising time
        barco.Fm = (barco.setl - barco.Fm) * delta/barco.Tsf + barco.Fm
        
        #The waterjet system is considered perfectly symetrical....
         
        #Moment generated by the waterjet around the USV CoM
        #first the setpoint for the waterjet orientation is fit within its
        #limits
        barco.setw = \
        ((barco.setw <= barco.thewjmax) & (barco.setw >= -barco.thewjmax))\
        * barco.setw\
        + barco.thewjmax * (barco.setw >barco.thewjmax ) \
        - barco.thewjmax * (barco.setw < -barco.thewjmax )
        
        #The waterjet turnrate is also model as a first order system. Tsw is 
        #model rising time   
        barco.thewj = (barco.setw - barco.thewj)* delta/barco.Tsw + barco.thewj
        
        #propulsion force components are calculate accordingly 
        barco.Fb[0] = barco.Fm * np.cos(barco.theta + barco.thewj)
        barco.Fb[1] = barco.Fm * np.sin(barco.theta + barco.thewj)

        #and the moment exerted       
        barco.M = - barco.Ac * barco.Fm * np.sin(barco.thewj)

        #then acceleration are calculated, including the damping forces
        #and other external forces...
        #x component
        barco.ab[0] = (barco.Fb[0]\
        - barco.mul * np.cos(barco.theta) \
        * (barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        - barco.mut * barco.ls * np.sin(barco.theta)\
        * (barco.vr[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))\
        - barco.mul2 * np.cos(barco.theta) \
        * np.sign(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        *(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))**2\
        - barco.mut2 * barco.ls * np.sin(barco.theta)\
        * np.sign(barco.vr[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))\
        * (barco.vr[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))**2\
        + extf[0])\
        / (barco.mb + barco.mA[0] * np.cos(barco.theta)\
        - barco.mA[1] * np.sin(barco.theta))
        
        #y component
        barco.ab[1] = (barco.Fb[1]\
        - barco.mul * np.sin(barco.theta) \
        *(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        - barco.mut * barco.ls * np.cos(barco.theta)\
        * (-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))\
        - barco.mul2 * np.sin(barco.theta)\
        *np.sign(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        * (barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))**2\
        - barco.mut2 * barco.ls * np.cos(barco.theta)\
        * np.sign(-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))\
        *(-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))**2\
        + extf[1])  / (barco.mb + barco.mA[1] * np.cos(barco.theta)\
        + barco.mA[0] * np.sin(barco.theta))
        
        #angular acceleration    
        barco.alfab = (barco.M - barco.mua * barco.ls * barco.wb -\
        np.sign(barco.wb) * barco.mua2 * barco.ls * barco.wb **2 - dfcm[1]\
        * (extf[0] * np.cos(barco.theta) +extf[1] * np.sin(barco.theta)) +\
        dfcm[0]\
        * (extf[1] * np.cos(barco.theta)\
        - extf[0] * np.sin(barco.theta)))/barco.Ib   
    
def tim_prop(barco,delta= 5e-4,extf = [0,0],dfcm = [0,0] ):
        '''
        This function models an inboard engine and rudder as propulsion 
        and guiding system. The USV is propelled by a conventional propeller
        The coupling between propulsion and guide system is 'weak'
        
        Departing from the forces a torques exerted on the
        USV, the function turn out the USV acceleration. 
        Whenever a instance of the class barco is created with the option 
        tipo = 'timon'(in fact watever value except 'waterject),
        the parameter propulsion of the instance points to 
        this function. Making the USV to adopt a propeller-rudder-type
        propulsion and guiding system.
        
        barco -> an instance of the class barco. (Rarely use because all 
        instance of barco has already its own propulsion system defined at
        construction time. Thus the function is mainly use as 
        barco_instance.propulsion(reamining parameters))
        
        delta -> step time for integration
        extf -> external force, (usually the las boom element strain)
        dfcm -> aplication point of the external force..- np.sign(barco.wb) * barco.mua2 * barco.ls * barco.wb **2. notice that only
                one external force can be applied
        '''
        #First, the function takes the current engine setpoint value and
        #and assures the setpoint value is within its limits...
        barco.setl =\
        ((barco.setl >= 0) & (barco.setl <= barco.pmax)) * barco.setl\
        + barco.pmax * (barco.setl >barco.pmax )        
        
        #the engine is modelled as a first order system. Tsf is the model
        #raising time
        barco.Fm = (barco.setl - barco.Fm) * delta/barco.Tsf + barco.Fm

         #The rudder and the USV are considered perfectly symetrical....
        
        #rudder orientation: first the setpoint for rudder orientation is fit
        #within its limits
        barco.setw = \
        ((barco.setw <= barco.thewjmax) & (barco.setw >= -barco.thewjmax))\
        * barco.setw\
        + barco.thewjmax * (barco.setw > barco.thewjmax ) \
        - barco.thewjmax * (barco.setw < -barco.thewjmax )
         
        
        barco.thewj = (barco.setw - barco.thewj)* delta/barco.Tsw + barco.thewj
        
        #propulsion force components. The effect of the rudder is include as
        #a damping forces adding a form factor to enlarge its effect...
        barco.Fb[0] = barco.Fm * np.cos(barco.theta) + barco.cfr * \
        (-barco.vr[0]*np.sin(barco.theta +barco.thewj)
        + barco.vr[1]*np.cos(barco.theta +barco.thewj))\
        * np.sin(barco.theta +barco.thewj)         
        barco.Fb[1] = barco.Fm * np.sin(barco.theta) - barco.cfr * \
        (-barco.vr[0]*np.sin(barco.theta +barco.thewj)
        + barco.vr[1]*np.cos(barco.theta +barco.thewj))\
        * np.cos(barco.theta +barco.thewj)

        #moment exerted
        barco.M = - barco.Ac * barco.cfr * \
        (barco.vr[0]*np.cos(barco.theta +barco.thewj)
        + barco.vr[1]*np.sin(barco.theta +barco.thewj))*np.sin(barco.thewj)

        
        #acceleration
        #x component
        barco.ab[0] = (barco.Fb[0]\
        - barco.mul * np.cos(barco.theta) \
        * (barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        - barco.mut * barco.ls * np.sin(barco.theta)\
        * (barco.vr[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))\
        - barco.mul2 * np.cos(barco.theta) \
        * np.sign(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        *(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))**2\
        - barco.mut2 * barco.ls * np.sin(barco.theta)\
        * np.sign(barco.vb[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))\
        * (barco.vr[0] * np.sin(barco.theta)\
        - barco.vr[1] * np.cos(barco.theta))**2\
        + extf[0])\
        / (barco.mb + barco.mA[0]*np.cos(barco.theta)\
        - barco.mA[1]*np.sin(barco.theta))
        #y component 
        barco.ab[1] = (barco.Fb[1]\
        - barco.mul * np.sin(barco.theta) \
        *(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        - barco.mut * barco.ls * np.cos(barco.theta)\
        * (-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))\
        - barco.mul2 * np.sin(barco.theta) \
        *np.sign(barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))\
        * (barco.vr[1] * np.sin(barco.theta)\
        + barco.vr[0] * np.cos(barco.theta))**2\
        - barco.mut2 * barco.ls * np.cos(barco.theta)\
        * np.sign(-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))\
        *(-barco.vr[0] * np.sin(barco.theta)\
        + barco.vr[1] * np.cos(barco.theta))**2\
        + extf[1])  / (barco.mb + barco.mA[1] * np.cos(barco.theta)\
        + barco.mA[0] * np.sin(barco.theta))
               
        
        #angular acceleration    
        barco.alfab = (barco.M - barco.mua * barco.ls * barco.wb -\
        np.sign(barco.wb) * barco.mua2 * barco.ls * barco.wb **2 - dfcm[1]\
        * (extf[0] * np.cos(barco.theta) +extf[1] * np.sin(barco.theta)) +\
        dfcm[0]\
        * (extf[1] * np.cos(barco.theta)\
        - extf[0] * np.sin(barco.theta)))/barco.Ib