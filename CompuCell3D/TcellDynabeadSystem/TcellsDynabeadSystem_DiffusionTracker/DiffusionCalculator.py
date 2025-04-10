#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:23:26 2025

@author: nora
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD

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





def universe_convert(cell_dictionary):
    
    n_mcs = len(cell_dictionary['x'])
    
    coords = np.array([cell_dictionary['x'], cell_dictionary['y']]) #[[x1, x2, ...], [y1, y2, ...]]
    coords_T = coords.T #[[x1, y1], [x2, y2], ...]    
    
    n_mcs = len(coords_T)
    if n_mcs < 2:
        raise ValueError("Not enough frames for MSD calculation!")
    
    coords = np.zeros((n_mcs, 1, 3))  # 1 particle, 3D coords
    coords[:, 0, :2] = coords_T  # Assign X and Y, Z stays 0

    
    u = mda.Universe.empty(1, trajectory = True)
    u.add_TopologyAttr('name', ['particle'])
    u.load_new(coords, order = 'fac') #frames, atoms, components
    
    return u


def msd_calculate(cell_dictionary):
    try:
        u = universe_convert(cell_dictionary)
        msd_analysis = EinsteinMSD(u, select='all', msd_type='xy', fft=True)
        msd_analysis.run()

        # Ensure MSD result is 2D before indexing
        if msd_analysis.results.timeseries.ndim != 2:
            raise ValueError("MSD result is not 2D as expected!")

        times = msd_analysis.results.timeseries[:, 0]
        msd = msd_analysis.results.timeseries[:, 1]
        return times, msd

    except Exception as e:
        # Catch and return the error
        return None, str(e)

print("\nDynabeads MSDs:")
for cell_id, data in Dynabeads.items():
    try:
        times, msd = msd_calculate(data)
        print(f"Bead {cell_id}: MSD at different times: {msd}")
    except Exception as e:
        print(f"Error calculating MSD for Bead {cell_id}: {str(e)}")







'''
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

'''
'''
subroutine dif(switch, nsamp)
if (switch.eq.0) then
    ntel=0
    dtime=dt*nsamp
    do i=1, tmax
        ntime(i)=0
        vacf(i)=0
        r2t(i)=0
    enddo    
  

else if (switch.eq.1) then
    ntel=ntel+1
    if (mod(ntel,it0).eq.0) then
        t0 = t0 + 1
        tt0 = mod(t0-1,t0max)+1
        time0(tt0) = ntel
        do i = , npart
            x0(i,tt0) = x(i)
            vx0(i,tt0) = vx(i)
        enddo
    endif

    do t = 1, min(t0,t0max)
        delt = ntel - time0(t) + 1 
        if (delt.lt.tmax) then:
            ntime(delt) = ntime(delt) + 1
            do i = 1, npart
                vacf(delt) = vacf(delt) + vx(i) * vx0(i, t)
                r2t(delt) = r2t(delt) + (x(i) - x0(i, t)**2)
            enddo
        endif
    enddo

else if (switch.eq.2) then
    do i = 1, tmax
        time = dtime * (i + 0.5)
        vacf(i) = vacf(i) / (npart*ntime(i))
        r2t(i) = r2t(i) / (npart*ntime(i))
    enddo
endif
return
end
'''      




    

        























            