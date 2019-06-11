# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:17:15 2016

@author: juan
Modelito de corrientes para el simulador de los arrastres
Es una funcion que asocia a cada punto del espacio un vector
velocidad vc este vector debe sumarse a la velocidad  del barco o de los
eslabones para simular el efecto que haría sobre ellos la corriente
Una segunda función permite dibujar los vectores velocidad de la corriente en
una zona determinada a la resolucion (dx,dy que se pida)
"""
import numpy as np
import matplotlib.pyplot as mpl

def gencor(pos,tiempo = 0,fun = lambda x,t=0: np.zeros(x.shape)):
    #x e y deben ser dos vectores (arrays de numpy) de datos las corrientes
    #se calcularan en todos lso pares de puntos x,y
    #pos debe ser un array de datos x,y (cada columna un par de datos)
    #vc = np.zeros(pos.shape)   
    vc = fun(pos,tiempo)
        
    return vc

def vercor(limites,nnodos,tiempo = 0,fun = lambda x,t=0: np.zeros(x.shape)):
    #calcula los valores de la corriente en un reticula (x,y) regular
    #limites debe ser un vector limites =[min_x,max_x,min_y,max_y]    
    #mola mil, mola mil, mola mil
    #nnodes= [puntos_x,puntosy] En los que se dividira la reticula
    escalas = np.array([limites[1] - limites[0],limites[3] - limites[2]])\
    /(np.array(nnodos) - 1.)   
    pc = np.indices([nnodos[0],nnodos[1]])* 1. 
    pc[0] = pc[0] * escalas[0] + limites[0]
    pc[1] = pc[1] * escalas[1] + limites[2]
    vc = fun(pc,tiempo)       
    mpl.quiver(pc[0],pc[1],vc[0],vc[1], width=0.002,units='width')
    
#here start a list of current field generators. Obviously the list could be
#extended adding more functions. So, please, please if you have whim add new
#functions and leave alone those already writen    

def camp0(pos, tiempo = 0):
    #no hay corriente    
    vc = np.zeros(pos.shape)
    return vc
    
def camp1(pos, tiempo = 0):
    #función de prueba para calcular campillos de velocidades
    #la funcion devuelve lo que traga pero espera arrays de numpy
    vc = pos.copy()

    vc[0] = 0.3 * (1-np.exp(-np.abs(pos[1])))*(pos[1]>=0)
    vc[1] = 0.*pos[1]    
    return vc
    
def camp2(pos, tiempo = 0):
    #función de prueba para calcular campillos de velocidades
    #la funcion devuelve lo que traga pero espera arrays de numpy
    vc = pos.copy()

    vc[0] = -0.3 * (1-np.exp(-np.abs(pos[1])))*(pos[1]>=0)
    vc[1] = 0.*pos[1]    
    return vc

campos =[camp0,camp1,camp2]    