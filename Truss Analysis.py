# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:42:58 2022

@author: hwlee
"""

from sympy import Symbol, symbols, solve
import numpy as np

# Basic Information of the Structural System
system = '2D Truss'

# Nodes
node1 = [0, 0]
node2 = [3, 0]
node3 = [0, 4]
node4 = [3, 4]
node5 = [6, 4]
node6 = [3, 7]

node_list = [i for i in dir() if ('node' in i) & ('_' not in i)]


# Elements [node_i, node_j, node_i와 x축과의 각도]
element1 = [1, 2, 0]
element2 = [1, 3, 90]
element3 = [2, 3, 135]
element4 = [2, 4, 90]
element5 = [3, 4, 0]
element6 = [4, 5, 0]
element7 = [3, 6, 45]
element8 = [4, 6, 90]
element9 = [5, 6, 135]

element_list = [i for i in dir() if ('element' in i) & ('_' not in i)]

# Degree of Freedom
if system == '2D Truss':
    dof = len(node_list) * 2
elif system == '3D Truss':
    dof = len(node_list) * 3

# Unknown Variables
E = Symbol('E')
A = Symbol('A')
K = np.zeros((dof,dof))

# Force
# F = np.zeros((len(node_list),2))
globals()['F{}1'.format(len(node_list))] = Symbol('F{}1'.format(len(node_list)))
globals()['F{}2'.format(len(node_list))] = Symbol('F{}2'.format(len(node_list)))

F = np.zeros((len(element_list), 2))

for i in range(len(element_list)):
    F[i][0] = globals()['F{}1'.format(i+1)]
    F[i][1] = globals()['F{}2'.format(i+1)]

F[1][1] = 0
F[2][1] = 0
F[5][2] = 300

# Displacement
globals()['u{}1'.format(len(node_list))] = Symbol('u{}1'.format(len(node_list)))
globals()['u{}2'.format(len(node_list))] = Symbol('u{}2'.format(len(node_list)))

u11 = 0
u12 = 0
u21 = 0
u22 = 0


# Length
for i in range(len(element_list)):
    node_i = globals()['element{}'.format(i+1)][0]
    node_j = globals()['element{}'.format(i+1)][1]
    
    x_i = globals()['node{}'.format(node_i)][0]
    y_i = globals()['node{}'.format(node_i)][1]
    x_j = globals()['node{}'.format(node_j)][0]
    y_j = globals()['node{}'.format(node_j)][1]

    globals()['l{}'.format(i+1)] = np.sqrt((x_i-x_j)**2 + (y_i-y_j)**2)

# Angles between local and global axes
for i in range(len(element_list)):
    globals()['lambda{}x'.format(i+1)] = np.cos(np.deg2rad(globals()['element{}'.format(i+1)][2]))
    globals()['lambda{}y'.format(i+1)] = np.cos(np.deg2rad(90 - globals()['element{}'.format(i+1)][2]))


# Stiffness  Matrices(function)
def create_element_k(lambda_x, lambda_y, l):
    
    k = 1/l*np.array([[lambda_x**2, lambda_x*lambda_y, -lambda_x**2, -lambda_x*lambda_y],
                        [lambda_x*lambda_y, lambda_y**2, -lambda_x*lambda_y, -lambda_y**2],
                        [-lambda_x**2, -lambda_x*lambda_y, lambda_x**2, lambda_x*lambda_y],
                        [-lambda_x*lambda_y, -lambda_y**2, lambda_x*lambda_y, lambda_y**2]])
    return k

for i in range(len(element_list)):
    lambda_x = globals()['lambda{}x'.format(i+1)]
    lambda_y = globals()['lambda{}y'.format(i+1)]
    l = globals()['l{}'.format(i+1)]
    globals()['k{}'.format(i+1)] = create_element_k(lambda_x, lambda_y, l)

# Member Stiffness Matrix -> Truss Stiffness Matrix
for i in range(len(element_list)):
    k2K_idx_1 = globals()['element{}'.format(i+1)][0] *2 - 2
    k2K_idx_2 = globals()['element{}'.format(i+1)][0] *2 - 1
    k2K_idx_3 = globals()['element{}'.format(i+1)][1] *2 - 2
    k2K_idx_4 = globals()['element{}'.format(i+1)][1] *2 - 1
    
    K[k2K_idx_1][k2K_idx_1] += globals()['k{}'.format(i+1)][0][0]
    K[k2K_idx_1][k2K_idx_2] += globals()['k{}'.format(i+1)][0][1]
    K[k2K_idx_1][k2K_idx_3] += globals()['k{}'.format(i+1)][0][2]
    K[k2K_idx_1][k2K_idx_4] += globals()['k{}'.format(i+1)][0][3]
    K[k2K_idx_2][k2K_idx_1] += globals()['k{}'.format(i+1)][1][0]
    K[k2K_idx_2][k2K_idx_2] += globals()['k{}'.format(i+1)][1][1]
    K[k2K_idx_2][k2K_idx_3] += globals()['k{}'.format(i+1)][1][2]
    K[k2K_idx_2][k2K_idx_4] += globals()['k{}'.format(i+1)][1][3]
    K[k2K_idx_3][k2K_idx_1] += globals()['k{}'.format(i+1)][2][0]
    K[k2K_idx_3][k2K_idx_2] += globals()['k{}'.format(i+1)][2][1]
    K[k2K_idx_3][k2K_idx_3] += globals()['k{}'.format(i+1)][2][2]
    K[k2K_idx_3][k2K_idx_4] += globals()['k{}'.format(i+1)][2][3]
    K[k2K_idx_4][k2K_idx_1] += globals()['k{}'.format(i+1)][3][0]
    K[k2K_idx_4][k2K_idx_2] += globals()['k{}'.format(i+1)][3][1]
    K[k2K_idx_4][k2K_idx_3] += globals()['k{}'.format(i+1)][3][2]
    K[k2K_idx_4][k2K_idx_4] += globals()['k{}'.format(i+1)][3][3]
    
# Solve!
    


