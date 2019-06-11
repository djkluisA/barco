# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 20:10:35 2015
Functions included in this module perform different calculations using 
bezier's curves.

functions:
1. bezier4p(p1,p2,v1,v2,tf,theta1,theta2,n):
2. pintar_bezier(r,v = None,color = 'r'):
3. longitud_curva(p1,p2,v1,v2,tf,theta1,theta2,dt):
See details for function use under the function definition
checked 06.03.2017
@author: juan
"""

import numpy as np
from scipy.misc import comb
import matplotlib.pyplot as pl

def bezier4p(p1,p2,v1,v2,tf,theta1,theta2,n):
    '''
    This function calculates points of a 3th degree Bezier´s curve, using
    four control points. Only the first and last one are supplied, the
    intermediate points are calculated to fit prescribed values of derivatives
    in the first and last control points.         
    
    p1 -> first control point, the curve starts at this point
    p2 -> last control point, the curve ends at this point
    p1 and p2 are defined as column vector, using a list with the following
    structure: p1 = [[x1],[x2]] and p2 = [[x2],[y2]]      
    
    v1 -> Derivative at point p1 (module)
    V2 -> Derivative at point p2 (module)    
    tf -> Final value of the parameter to cover the curve. Usually 
          Bezier's curves are defined using a parameter that varies between 0 
          and 1, Here the parameter is redefined and it varies between 0 and tf
    
    theta1 -> 
    '''
    
    t = np.linspace(0,1,num = n)*tf
        
    #definir los puntos de control 
    p = np.array([p1, p1 + np.array([[v1*tf/3*np.cos(theta1)],[v1*tf/3*np.sin(theta1)]]),\
    p2 - np.array([[v2*tf/3*np.cos(theta2)],[v2*tf/3*np.sin(theta2)]]),p2])
    
    #obtencion puntos de control de la odografa
    d = (p[1:]-p[:-1])        
    
    r = np.zeros((2,n))
    v = np.zeros((2,n))    
            
    for i in range(4):    
        r += tf**-3 * comb(3,i) * t ** i * (tf - t) ** (3-i) * p[i]
    
    
    for i in range(3):
        v += tf**-3 * comb(2,i) * t ** i * (tf - t) ** (2-i) * 3*d[i]
        
    return r,v,t

def pintar_bezier(r,v = None,color = 'r'):
    
    pl.plot(r[0],r[1],color)
    pl.hold(True)
    if v is not None:
        for i,j,k,l in zip(r[0],r[1],v[0],v[1]):
            pl.arrow(i,j,k,l)
            
def longitud_curva(p1,p2,v1,v2,tf,theta1,theta2,dt):
    #es un calculo aproximado integrando. El paso de integracion lo da dt
    p = np.array([p1, p1 + np.array([[v1*tf/3*np.cos(theta1)],\
    [v1*tf/3*np.sin(theta1)]]), p2 - np.array([[v2*tf/3*np.cos(theta2)],\
    [v2*tf/3*np.sin(theta2)]]),p2])    
    d = (p[1:]-p[:-1])
    s = 0
    for t in np.arange(0,tf,dt):
        v = 0        
        for i in range(3):
            v += tf**-3 * comb(2,i) * t ** i * (tf - t) ** (2-i) * 3*d[i]        
        s = s + np.sqrt(v[0]**2 + v[1]**2)*dt
    return s