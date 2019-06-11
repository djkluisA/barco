# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:02:05 2016
This module containts several methods to calculate the strains matrix used
to obtain the strains in the links of the booms.

NOTE: boom.py EMPLOY THIS MODULE TO IMPLEMENT THE STRAIN MATRIX of cadena class

There are to kinds of methods. 
The first one allows to initialise the strain matrix, the second one allows
to calculated the matrix and the strain acoording to the specific scenario.

initialisation methods
The methods implemented are:
    short booms (less than 30 links)
            1. mtr_ini(cad):
    long booms (30 links or more)
            1. mtrsp_ini(cad)        

NOTE: cadena class initialisate by default the strain matrix according to the
length of the boom and the presence of towing cables... It does it directly
without making use of these initialisation functions

Matrix and strain calculation:
Towing forces applying to the boom but not towing ships are considered
    Boom with or without towing cables are equivalent from the point of view 
    of strain calculations. The boom class implement by default this functions
    when the lenght of both towing cables are 0.
        short booms (less than 30 links)
            2. mtr_s(cad):
        long booms (30 links or more)
            2. mtrsp_s(cad)

Towing ship at one end of the boom, the other end attached to a deck
In this scenario it is asumed that the ship is attached to the boom using a
towing cable and the other end is directly attached to the deck

    Left end attached to the deck
        short booms (less than 30 links)            
            2. mtrl_s(cad,barcod)
        long booms (30 links ore more)
            2. mtrlsp_s(cad,barcod)
    Right end attached to the deck
        short booms (less than 30 links)
            2. mtrr_s(cad,barcoi)
        long booms (30 links ore more)
            2.mtrrsp_s(cad,barcoi)

Towing ships on both ends of the boom
This is the default situation when lenght of both towing cables is greater
than 0    
        short booms (less than 30 links)            
            2. mtrb_s(cad,barcoi,barcod)
        long booms (30 links ore more)
            2. mtrbsp_s(cad,barcoi,barcod)
            
            
reviewed: 06.03.2017

@author: juan
"""

import numpy as np
import numpy.linalg.linalg as npl
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

###############################################################################
#                        Initialisation fncs                                  #
###############################################################################

#full strain matrix############################################################
def mtr_ini(cad):
    '''
    Strain matrix initialisation for short, less that 30 links, booms.
    
    It creates matrices of zeros to subsequently calculate  
    the strain components (x,y) in the ends of the links.    
    
    The matrix dimension is 2 times plus 2 the number of links of the boom. 
    Notice  that the strain in the right-end of link number i is equal to the
    left-end strain for link number i+1. then, for n links the total number of
    of strains will be n+1, taking into account the strain at the tips of the
    boom. This fuction is pointed by default for (short) instances of the class
    cadena when they are created 
    
    cad -> An instance of the class cadena, i.e, a boom.
    
    M <- Strain Matrixreturn M,B,T 
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B
    '''     
    M = np.zeros((2 * cad.esl + 2, 2 * cad.esl + 2))
    #first and last terms initialisation
    M[0,0]= 1
    M[1,1]= 1
    M[-2,-2]= 1
    M[-1,-1]= 1
    B = np.zeros(2 * cad.esl + 2)
    T = np.zeros(2 * cad.esl + 2)
    return M, B, T

#sparse strain matrix##########################################################
def mtrsp_ini(cad):
    '''
    Strain (sparse) matrix initialisation for large, more than 30 links, booms.  
    
    It creates matrices of zeros to subsequently calculate  
    the strain components (x,y) in the ends of the links.    
    
    The matrix dimension is 2 times plus 2 the number of links of the boom. 
    Notice  that the strain in the right-end of link number i is equal to the
    left-end strain for link number i+1. then, for n links the total number of
    of strains will be n+1, taking into account the strain at the tips of the
    boom. This fuction is pointed by default for (large) instances of the class
    cadena when they are created: (cadena.mtr_i -> smatrix.mtrspc_ini) 
    
    
    cad -> An instance of the class cadena, i.e, a boom.
    
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of nreturn M,B,T  links. M*T = B
    '''              
    #an emtpy  sparse (csr) matrix is built. Dimension:  2*number of links
    #plus 2. The csr (compressed sparse row format)  is the most
    #suitable format for the structure of the strains matrix when working with 
    #spsolve... (see SciPy documentation for details )
   
    M = csr_matrix((2 * cad.esl +2,2 * cad.esl +2)) 
    
    #pointer to the non-null elements of M 
    M.indptr[1:3] = [4,8]
    M.indptr[3:-2] = np.arange(14,14+6*(M.indptr.size-5),6)
    M.indptr[-2:] = [M.indptr[-3] + 4, M.indptr[-3] + 8]      
    
    #Indices of the non-null elements
    M.indices= np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5],\
    dtype = 'int32')
    for i in range(4,M.indptr.size-4,2):
        M.indices= np.concatenate((M.indices,i-2 + M.indices[8:20]))
    M.indices= np.concatenate((M.indices, i + M.indices[0:8]))
    M.data = np.zeros(M.indices.size)
    #first and last terms initialisation    
    M[0,0]= 1
    M[1,1]= 1
    M[-2,-2]= 1
    M[-1,-1]= 1
    B = np.zeros(2 * cad.esl +2)
    T = np.zeros(2 * cad.esl +2)
       
    return M, B, T 
    


###############################################################################
#                  funtions for strains calculation                           #        
###############################################################################        

#external forces no towing ships###############################################

#full strain matrix#########################################################
def mtr_s(cad):
    '''
    This function calculates the components of the strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom.
    
    The calculation are performed departing from the current state of the boom,
    pose, speed and values of cad.Fi and cad.Fd, ie, the forces applied at the
    end of the boom.
    
    This fuction is pointed by default for (short<30 links) instances of the
    class cadena when they are created:(cadena.mtrc_s -> smatrix.mtrc_s) if
    there are not towing cables; both cables length == 0. 
    
    for other scenarios such as towing ships or left end of the boom achored
    see below.
         
    M ->      a previously calculated/inisialisated strain matrix
    B ->      a previously calculaled/inisialisated vector of independent terms
    T ->      a previously calculated/initialisated vector of strain components
    cad ->    An instance of the class cadena, i.e, a boom.
    
              
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal) #inverse of added mass
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping factor
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #rotating damping factor    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
 

        
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]  
        
      
    for i in range(0,cad.normal.shape[1]-1): 
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
        [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]])
    
    #independent terms 
        
    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
        
    #The forces are applied straigthforward to ends of the boom so,
    #their values are  the strains on the left of the first link and on the
    #right of the last one.
        
        
    cad.B[0] = cad.Fi[0]
    cad.B[1] = cad.Fi[1]

    
    cad.B[-2] = cad.Fd[0]
    cad.B[-1] = cad.Fd[1]
       
      
    #solve the linear system M * T = B  
    cad.T = np.linalg.linalg.solve(cad.M,cad.B)
    
#sparse strain matrix##########################################################    
def mtrsp_s(cad):
   '''
    This function calculates the components of the strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom, using sparse matrices and sparse algebra.
    
    The calculation are performed departing from the current state of the boom,
    pose, speed and values of cad.Fi and cad.Fd, ie,the forces applied at the
    end of the boom.
    
    This fuction is pointed by default for (long>=30 links) instances of the
    class cadena when they are created (cadena.mtrc_s -> smatrix.mtrc_s)
    if there are not towing cables; both cables lenght == 0. 
    
    for other scenarios such as towing ships or left end of the boom achored
    see below.
         
    M ->      a previously calculated/inisialisated strain matrix
    B ->      a previously calculaled/inisialisated vector of independent terms
    T ->      a previously calculated/initialisated vector of strain components
    cad ->    An instance of the class cadena, i.e, a boom.
    
              
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B
    '''     
   
     #some terms frecuently used are calculated first...   
   mcA = 1./(cad.m + cad.mA * cad.normal)
   T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
   T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
   #dumping factors
   factor1 = mcA *\
   ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
   cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
   (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
   cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    
        
   factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
   cad.w * cad.normal /cad.I) 
    
        
    
   T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
   T2ya = T1ya[:-1] + T1ya[1:]
    
   #strain matrix inner non zero terms 
   for i in range(0,cad.normal.shape[1]-1):  
       cad.M.data[cad.M.indptr[2*i+2]:cad.M.indptr[2*i+4]] = \
       [T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
       T1ya[i+1],T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
       T1ya[i+1], T1xa[1,i+1]]

   #Independent terms
   cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
   + factor2[0,1:] +factor2[0,0:-1]\
    
   cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
   + factor2[1,1:] +factor2[1,0:-1]
        
    
   
   cad.B[0] = cad.Fi[0]
   cad.B[1] = cad.Fi[1]

   cad.B[-2] = cad.Fd[0]
   cad.B[-1] = cad.Fd[1] 
      
           
   cad.T = spsolve(cad.M,cad.B)
  
#a towing ship in one end and the other moored to the dock#####################
  
               
#left end moored, right attached to a ship#####################################
               
#full strain matrix############################################################               
def mtrl_s(cad, barcod):
    '''
    This function calculates the components of the strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom. For this scenario, the left end of the boom is attached
    to the deck
    M ->      a previously calculated/inisialisated strain matrix
    B ->      a previously calculaled/inisialisated vector of independent terms
    T ->      a previously calculated/initialisated vector of strain components
    cad ->    An instance of the class cadena, i.e, a boom.
    barcod -> an instance of the class barco that represent a model of a towing
              ship, attached to the right-end of the boom                  
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal) #inverse of added mass
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
        
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
 

        
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]  
        
      
    for i in range(0,cad.normal.shape[1]-1): 
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
        [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]])
    
    #independent terms 
        
    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
    #The following elements of the strain matrix should have been fit as it's
    #shown here, somewhere before callng to this function..
    #cad.M[0:2,0:4] = np.array([[1.,0.,1.,0.],[0.,1.,0.,1.]])    

                        
    
    #vector from the tip of the boom to the ship    
    vcbld =  cad.cms[:,-1] - cad.para[:,-1] * cad.L \
    - barcod.csbt([-barcod.ls/2,0]) - barcod.pb
    #distance between the tip of the boom and the ship
    distd = npl.norm(vcbld)
    dt = distd >= cad.cbl[1]
       
    #if  dt = 1 then, the cable is tight
    #unit vector 
    vcbld = vcbld/(distd + (distd == 0.))
       
  
    
    cad.B[-2] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[0]
    cad.B[-1] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[1]
       
    cad.T = np.linalg.linalg.solve(cad.M,cad.B)    

#sparse strain matrix##########################################################
def mtrlsp_s(cad,barcod):
    '''
    This function calculates the components of the  strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom, using sparse matrices and sparse algebra.
    For this scenario, two ships cooperate to tow the boom
    cad ->    An instance of the class cadena, i.e, a boom.
    cad.M ->  a previously calculated/inisialisated strain matrix
    cad.B ->  a previously calculaled/inisialisated vector of independent terms
    cad.T ->  a previously calculated/initialisated vector of strain components
    barcoi -> an instance of the class barco that represent a model of a towing
              ship, attached to the left-end of the boom 
    barcod -> an instance of the class barco that represent a model of a towing
              ship, attached to the right-end of the boom                  
    cad.M <- Strain Matrix
    cad.B <- Vector of independent terms
    cad.T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal)
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #codiciones de cierre y resistencia al giro    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
    
        
    
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]
    
  
    for i in range(0,cad.normal.shape[1]-1):
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M.data[cad.M.indptr[2*i+2]:cad.M.indptr[2*i+4]] = \
        [T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1],T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]

    #independent terms 

    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
    #The following elements of the strain matrix should have been fit as it's
    #shown here, somewhere before callng to this function..
    #cad.M[0:2,0:4] = np.array([[1.,0.,1.,0.],[0.,1.,0.,1.]]) 
    
    
      
    #vector from the tip of the boom to the ship on the right
    vcbld =  cad.cms[:,-1] - cad.para[:,-1] * cad.L \
    - barcod.csbt([-barcod.ls/2,0]) - barcod.pb
    #distance between the tip of the boom and hte ship on the right
    distd = npl.norm(vcbld)
    dt = distd >= cad.cbl[1]
    #if  dt = 1 then, the cable is tight
    #unit vector    
    vcbld = vcbld/(distd + (distd == 0.))
       
   
    cad.B[-2] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[0]
    cad.B[-1] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[1]
               
    cad.T = spsolve(cad.M,cad.B)

#left end attached to a ship, right end moored#################################

#full strain matrix############################################################
     
def mtrr_s(cad, barcoi):
    '''
    This function calculates the components of the strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom. For this scenario, the right end of the boom is attached
    to the deck
    M ->      a previously calculated/inisialisated strain matrix
    B ->      a previously calculaled/inisialisated vector of independent terms
    T ->      a previously calculated/initialisated vector of strain components
    cad ->    An instance of the class cadena, i.e, a boom.
    barcod -> an instance of the class barco that represent a model of a towing
              ship, attached to the right-end of the boom                  
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal) #inverse of added mass
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
        
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
 

        
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]  
        
      
    for i in range(0,cad.normal.shape[1]-1): 
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
        [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]])
    
    #independent terms 
        
    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
    #The following elements of the strain matrix should have been fit as it's
    #shown here, somewhere before callng to this function..
    #cad.M[-2:,-4:] = np.array([[1.,0.,1.,0.],[0.,1.,0.,1.]])    

                        
    
    #vector from the tip of the boom to the ship on the left  
    vcbli =  barcoi.pb + barcoi.csbt([-barcoi.ls/2,0])\
    -cad.para[:,0] * cad.L - cad.cms[:,0] 
    #distance between the tip of the boom and the ship
    disti = npl.norm(vcbli)
    dt = disti >= cad.cbl[0]
       
    #if  dt = 1 then, the cable is tight
    #unit vector 
    vcbli = vcbli/(disti + (disti == 0.))
       
      
    
    cad.B[0] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[0]
    cad.B[1] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[1]
       
    cad.T = np.linalg.linalg.solve(cad.M,cad.B)

#sparse strain matrix##########################################################
def mtrrsp_s(cad,barcoi):
    '''
    This function calculates the components of the  strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom, using sparse matrices and sparse algebra.
    For this scenario, two ships cooperate to tow the boom
    cad ->    An instance of the class cadena, i.e, a boom.
    cad.M ->  a previously calculated/inisialisated strain matrix
    cad.B ->  a previously calculaled/inisialisated vector of independent terms
    cad.T ->  a previously calculated/initialisated vector of strain components
    barcoi -> an instance of the class barco that represent a model of a towing
              ship, attached to the left-end of the boom 
    barcod -> an instance of the class barco that represent a model of a towing
              ship, attached to the right-end of the boom                  
    cad.M <- Strain Matrix
    cad.B <- Vector of independent terms
    cad.T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal)
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #codiciones de cierre y resistencia al giro    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
    
        
    
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]
    
  
    for i in range(0,cad.normal.shape[1]-1):
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M.data[cad.M.indptr[2*i+2]:cad.M.indptr[2*i+4]] = \
        [T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1],T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]

    #independent terms 

    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
    #The following elements of the strain matrix should have been fit as it's
    #shown here, somewhere before callng to this function..
    #cad.M[-2:,-4:] = np.array([[1.,0.,1.,0.],[0.,1.,0.,1.]])   
   
    
    #vector from the tip of the boom to the ship on the left 
    vcbli =  barcoi.pb + barcoi.csbt([-barcoi.ls/2,0])\
    -cad.para[:,0] * cad.L - cad.cms[:,0]
    #distance between the tip of the boom and hte ship on the left
    disti = npl.norm(vcbli)
    dt = disti >= cad.cbl[0]
        
       
    vcbli = vcbli/(disti + (disti == 0.))  
    
             
    cad.B[0] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[0]
    cad.B[1] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[1]
        
   
       
    
               
    cad.T = spsolve(cad.M,cad.B)
    
    
    
#two ships towing the boom from its tips#######################################

#full strain matrix############################################################    
def mtrb_s(cad, barcoi, barcod):
    '''
    This function calculates the components of the strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom.
    For this scenario, two ships cooperate to tow the boom.
    M ->      a previously calculated/inisialisated strain matrix
    B ->      a previously calculaled/inisialisated vector of independent terms
    T ->      a previously calculated/initialisated vector of strains components
    cad ->    An instance of the class cadena, i.e, a boom.
    barcoi -> an instance of the class barco that represent a model of a towing
              ship, attached to the left-end of the boom    
    barcoi -> an instance of the class barco that represents a model of a towing
              ship, attached to the left-end of the boom
              
    M <- Strain Matrix
    B <- Vector of independent terms
    T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal) #inverse of added mass
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #codiciones de cierre y resistencia al giro    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
 

        
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]  
        
      
    for i in range(0,cad.normal.shape[1]-1): 
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
        [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]])
    
    #independent terms 
        
    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
        

    #vector from the tip of the boom to the ship on the left  
    vcbli =  barcoi.pb + barcoi.csbt([-barcoi.ls/2,0])\
    -cad.para[:,0] * cad.L - cad.cms[:,0]
    #distance between the tip of the boom and hte ship on the left
    disti = npl.norm(vcbli)
    dt = disti >= cad.cbl[0]
        
        
    vcbli = vcbli/(disti + (disti == 0.))  
        
             
    cad.B[0] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[0]
    cad.B[1] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[1]
    
    #vector from the tip of the boom to the ship on the right                                     
    vcbld =  cad.cms[:,-1] - cad.para[:,-1] * cad.L \
    - barcod.csbt([-barcod.ls/2,0]) - barcod.pb 
    #distance between the tip of the boom and hte ship on the right 
    distd = npl.norm(vcbld)
    dt = distd >= cad.cbl[1]
       
    vcbld = vcbld/(distd + (distd == 0.))
       
    
       
    cad.B[-2] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[0]
    cad.B[-1] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[1]
      
    #print 'izq', disti, 'dcha', distd ,'\n'      
    cad.T = np.linalg.linalg.solve(cad.M,cad.B)
 
#sparse strain matrix##########################################################    
def mtrbsp_s(cad,barcoi,barcod):
    '''
    This function calculates the components of the  strain matrix M, 
    the components of the vector B of independent terms  and solve the system
    M*T = B to obtain a vector T of components the strains at the ends of the
    links of the boom, using sparse matrices and sparse algebra.
    For this scenario, two ships cooperate to tow the boom
    cad ->    An instance of the class cadena, i.e, a boom.
    cad.M ->  a previously calculated/inisialisated strain matrix
    cad.B ->  a previously calculaled/inisialisated vector of independent terms
    cad.T ->  a previously calculated/initialisated vector of strain components
    barcoi -> an instance of the class barco that represent a model of a towing
              ship, attached to the left-end of the boom 
    barcod -> an instance of the class barco that represent a model of a towing
              ship, attached to the right-end of the boom                  
    cad.M <- Strain Matrix
    cad.B <- Vector of independent terms
    cad.T <- Vector of strains componets, the order is: [Tx0,Ty0,Tx1,Ty1,...
         Txi,Tyi,...Txn,Tyn] for a boom of n links. M*T = B       
    '''
    #some terms frecuently used are calculated first...
    mcA = 1./(cad.m + cad.mA * cad.normal)
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #damping forces
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #codiciones de cierre y resistencia al giro    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
    
        
    
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]
    
  
    for i in range(0,cad.normal.shape[1]-1):
    #terms of the strain matrix corresponding to inner links. the terms
    #corresponding to the links located at the end of the boom are calculated
    #later
        cad.M.data[cad.M.indptr[2*i+2]:cad.M.indptr[2*i+4]] = \
        [T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1],T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]

    #independent terms 

    cad.B[range(2,cad.B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    cad.B[range(3,cad.B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
      
   
    
    #vector from the tip of the boom to the ship on the left 
    vcbli =  barcoi.pb + barcoi.csbt([-barcoi.ls/2,0])\
    -cad.para[:,0] * cad.L - cad.cms[:,0]
    #distance between the tip of the boom and hte ship on the left
    disti = npl.norm(vcbli)
    dt = disti >= cad.cbl[0]
        
       
    vcbli = vcbli/(disti + (disti == 0.))  
    
   
             
    cad.B[0] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[0]
    cad.B[1] = -cad.cbl[2] * (disti -cad.cbl[0]) * dt * vcbli[1]
        
    #vector from the tip of the boom to the ship on the right
    vcbld =  cad.cms[:,-1] - cad.para[:,-1] * cad.L \
    - barcod.csbt([-barcod.ls/2,0]) - barcod.pb
    #distance between the tip of the boom and hte ship on the right
    distd = npl.norm(vcbld)
    dt = distd >= cad.cbl[1]
       
    vcbld = vcbld/(distd + (distd == 0.))
       
    
    
    cad.B[-2] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[0]
    cad.B[-1] = -cad.cbl[2] * (distd -cad.cbl[1]) * dt * vcbld[1]
               
    cad.T = spsolve(cad.M,cad.B)
       