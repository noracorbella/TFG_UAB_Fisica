#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:23:26 2025

@author: nora
"""

import numpy as np


positions_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.txt"


delta_t = 100 #we will calculate D every 100MCS


TCells = {}
Dynabeads = {}



with open(positions_file,'r') as f:
    next(f)
    for l in f:
        mcs, cell_id, cell_type, x, y = map(float, l.split())
        cell_id = int(cell_id)
        
        if cell_type == 1:
            if cell_id not in TCells:
                
                TCells[cell_id] = {'mcs': [], 'x': [], 'y': []}
            TCells[cell_id]['mcs'].append(mcs)
            TCells[cell_id]['x'].append(x)
            TCells[cell_id]['y'].append(y)
            
        elif cell_type == 2:
            if cell_id not in Dynabeads:
                
                Dynabeads[cell_id] = {'mcs': [], 'x': [], 'y': []}
            Dynabeads[cell_id]['mcs'].append(mcs)
            Dynabeads[cell_id]['x'].append(x)
            Dynabeads[cell_id]['y'].append(y)





def calculate_diffusion(Dictionary):
    
    all_D = []
        
    for Id in Dictionary:
        Data = Dictionary[Id]
        mcs = Data['mcs']
        x = Data['x']
        y = Data['y']
        
        SD = []
        
        for t_0 in range(len(mcs) - delta_t):
            dx = x[t_0 - delta_t] - x[t_0]
            dy = y[t_0 - delta_t] - y[t_0]
            
            SD.append(dx**2 + dy**2)
            
        if SD:
            MSD = np.mean(SD)
            D = MSD / (4 * delta_t)
            all_D.append(D)
    
    if all_D:
        D = np.mean(all_D)
        return D
    else:
        return 0



print("D =", calculate_diffusion(Dynabeads))
    


        























            