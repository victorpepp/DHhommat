# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:13:55 2017

@author: pepp

This code aims to evaluate the Homogenous Forward Kinematics Matrix for a
4-degrees of freedom robotic arm, based on an anthropomorphic manipulator.
"""

from sympy.interactive import printing
printing.init_printing(use_latex=True)
import numpy as np
from sympy import *

"""
Global declarations
"""
length = [315.53, 210.33, 454.67, 707.00, 894.67] # assessed physically
dof = 4

def sine(angle):
    if type(angle) == float:
        return(round(sin(angle), 2))
    return sin(angle)

def cosine(angle):
    if type(angle) == float:
        return(round(cos(angle), 2))
    return cos(angle)

"""
Denavit-Hartenberg parameters
"""

ad = [0., length[2], length[3], length[4]]
alpha = [0.5*round(np.pi, 2), 0., 0., 0.]
dd = [length[0] + length[1], 0., 0., 0.]
theta = list(symbols('theta_1 theta_2 theta_3 theta_4'))
theta[1] = theta[1] + round(0.5*np.pi, 2)

def hommat(i, ad, alpha, dd, theta):
    i-=1
    hommat = Matrix([[cosine(theta[i]), -1*sine(theta[i])*cosine(alpha[i]), sine(theta[i])*sine(alpha[i]), ad[i]*cosine(theta[i])],
           [sine(theta[i]), cosine(theta[i])*cosine(alpha[i]), (-1)*cosine(theta[i])*sine(alpha[i]), ad[i]*sine(theta[i])],
           [0., sine(alpha[i]), cosine(alpha[i]), dd[i]],
           [0., 0., 0., 1.]])
    return hommat

"""
Evaluation
"""

hom = hommat(1, ad, alpha, dd, theta)

for i in range(2, dof):
    hom = hom*hommat(i, ad, alpha, dd, theta)