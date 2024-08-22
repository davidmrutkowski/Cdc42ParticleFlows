"""
 * Copyright (C) 2024 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
"""

""" 
    This code calculates a heatmap of the widths of the Cdc42-GTP
    patches from the 3d sphere simulations as in Fig 3C.
    Directory and file location of csv's must match that which is 
    expected based on base_file_list and folder_list.
    Original data can be provided upon request.
""" 

import math
import tkinter as tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

def calc_dist(list1, list2):
    temp_x = list1[0] - list2[0]
    temp_y = list1[1] - list2[1]
    temp_z = list1[2] - list2[2]
    
    return math.sqrt(temp_x*temp_x + temp_y*temp_y + temp_z*temp_z)
    
def calc_sphericalSA(arclength, radius):
    theta = arclength / radius
    
    h = radius - radius*math.cos(theta)
    
    return 2.0*math.pi*radius*h

def sortBy(e):
    start_string = 'D_'
    end_string = '_endo-rate_'
    pos = e.rfind(start_string)
    posExt = e.rfind(end_string)
    name = e[pos+len(start_string):posExt]
    return name
    
 
def calc_mag(vector):
    magnitude = 0.0
    for i in range(0, len(vector)):
        magnitude += vector[i]*vector[i]
        
    return math.sqrt(magnitude)
    
def get_unit_vec(vector):
    curr_mag = calc_mag(vector)
    
    return [element / curr_mag for element in vector]
    
def calcArcLengthDistance(posi, posj, sphereRadius):
    currDistance = calc_dist(posi, posj);
    
    val_asin = 0.5 * currDistance / sphereRadius
    if val_asin >= 1.0 and val_asin < 1.01:
        val_asin = 1.0
        
    arclength = 2.0*sphereRadius * math.asin(val_asin);
    
    return arclength
    

def getValueFromString(s, start_str, end_str):
    start_pos = s.find(start_str)
    
    #tmp_str = s[start_pos + len(start_str):]
    
    end_pos = s[start_pos + len(start_str):].find(end_str) + start_pos + len(start_str)
    
    if start_pos == -1 and end_pos == -1:
        return s
    elif end_pos == -1:
        return s[start_pos+len(start_str):]
    else:
        return s[start_pos+len(start_str):end_pos]
        

def getEdgeList(org_list):
    edge_list = []
    
    for i in range(1, len(org_list)):
        left_val = org_list[i-1]
        right_val = org_list[i]
        
        avg_val = 0.5*(left_val + right_val)
        
        edge_list.append(avg_val)
        
    # now need to find additional two edges beyond org_list on each side
    dist_left = math.log10(edge_list[0]) - math.log10(org_list[0])
    
    print(edge_list[0], org_list[0], math.log10(org_list[0]), dist_left)
    edge_list = [pow(10.0, math.log10(org_list[0])-dist_left)] + edge_list
    
    dist_right = math.log10(org_list[-1]) - math.log10(edge_list[-1])
    edge_list.append(pow(10.0, math.log10(org_list[-1]) + dist_right))
    
    return edge_list
    
def getMeshLists(xlist, ylist):
    x_mesh_list = []
    y_mesh_list = []
    
    for j in range(0, len(xlist)):
        tmp_x_list = []
        tmp_y_list = []
        
        for i in range(0, len(ylist)):
            if len(x_mesh_list) <= i:
                x_mesh_list.append([xlist[j]])
            else:
                x_mesh_list[i].append(xlist[j])
                
            if len(y_mesh_list) <= i:
                y_mesh_list.append([ylist[i]])
            else:
                y_mesh_list[i].append(ylist[i])
        
    return x_mesh_list, y_mesh_list
    
   
initial_radius = 1.8

exclusion_center = [0.0, 0.0, initial_radius]

folder_list = []



# need to change DD value for each of these!
base_file_list = ['RDist_0.00020_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.00045_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.00100_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.00224_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.00500_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.01118_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.02500_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.05590_DD_0.00030__diffFont_frameCutoff0.csv', 'RDist_0.12500_DD_0.00030__diffFont_frameCutoff0.csv']


folder_list = ['D9.0E-5', 'D3.0E-4', 'D1.0E-3', 'D3.3E-3', 'D1.1E-2', 'D3.7E-2']

DD_values = []
rD_values = []

for i in range(0, len(base_file_list)):
    curr_rD = float(getValueFromString(base_file_list[i], 'RDist_', '_DD'))
    
    rD_values.append(curr_rD)

rD_values = sorted(rD_values)  
    
for i in range(0, len(folder_list)):
    curr_DD = float(folder_list[i][1:])
    
    DD_values.append(curr_DD)

DD_values = sorted(DD_values)


rD_edge_values = getEdgeList(rD_values)
DD_edge_values = getEdgeList(DD_values)

DD_mesh, rD_mesh = np.meshgrid(DD_edge_values, rD_edge_values)

DD_mesh_noedge, rD_mesh_noedge = np.meshgrid(DD_values, rD_values)

C = np.zeros((len(rD_values), len(DD_values)))


file_count = 0


for folder in folder_list:
    curr_DD_folder = float(folder[1:])
    
    curr_DD_file = curr_DD_folder
    #correction between foldername and file name
    if curr_DD_folder == 1.1E-2:
        curr_DD_file = 1.11E-2
    elif curr_DD_folder == 3.3E-3:
        curr_DD_file = 3.33E-3
    
    for filename in base_file_list:
        curr_pos_DD = filename.find("DD_")
        curr_pos_Run3 = filename.find("__diffFont")
        

        real_filename = filename[:curr_pos_DD + len("DD_")] + "{:.5f}".format(curr_DD_file) + filename[curr_pos_Run3:]
        
        
        filepath = folder + "/" + real_filename
        
        frame_count = -1
        frame_cutoff = 0

        try:
            with open(filepath) as fp:
                particle_avg_list = []
                
                line = fp.readline()
                
                line = fp.readline()
                
                
                
                arclength_list = []
                GAP_list = []
                
                while line:
                    split_string = line.split(',')
                    
                    curr_arclength = float(split_string[0])
                    
                    curr_Cdc42_GTP = float(split_string[2])
                    
                    curr_GAP = float(split_string[5])
                    curr_Cdc42_GDP = float(split_string[1])
                    
                    
                    arclength_list.append(curr_arclength)
                    GAP_list.append(curr_GAP)
                    
                    line = fp.readline()
                
                cutpoints = 20
                # cut off last 10 points because noisy
                GAP_list = GAP_list[:-cutpoints]
                arclength_list = arclength_list[:-cutpoints]
                
                arclength_width = 0.0
                left_GAP = np.mean(GAP_list[0:5])
                
                right_GAP = np.mean(GAP_list[len(GAP_list)-50:])
                
                avg_GAP = 0.5*(left_GAP + right_GAP)
                
                """for i in range(0, len(arclength_list)):
                    if GAP_list[i] >= avg_GAP:
                        arclength_width = arclength_list[i]
                        
                        break"""
                
                avg_window = 11
                
                for i in range(len(arclength_list)-1, -1, -1):
                    left_edge = i - avg_window
                    right_edge = i + avg_window
                    
                    if right_edge >= len(arclength_list):
                        right_edge = len(arclength_list) - 1
                        
                    if left_edge < 0:
                        left_edge = 0
                        
                    avg_GAP_value = np.mean(GAP_list[left_edge:right_edge])
                    
                    if avg_GAP_value <= right_GAP*0.90:
                    #if avg_GAP_value <= right_GAP*0.90 or avg_GAP_value >= right_GAP*1.10:
                        arclength_width = arclength_list[i]
                        
                        break
                
                
                first_pos = filepath.find("RDist_")
                
                tmp_str = filepath[first_pos+6:]
                
                second_pos = tmp_str.find("_")
                
                curr_rD = float(tmp_str[:second_pos])

                
                curr_pos_rD = rD_values.index(curr_rD)
                curr_pos_DD = DD_values.index(curr_DD_folder)
                
                print(curr_rD, curr_DD_folder, left_GAP, right_GAP, avg_GAP, arclength_width, curr_pos_rD, curr_pos_DD)
                #print(curr_rD, curr_DD_folder, arclength_width)
                
                #if arclength_width > 1.1:
                if arclength_width > 1.6:
                    C[curr_pos_rD, curr_pos_DD] = arclength_width
                else:
                    C[curr_pos_rD, curr_pos_DD] = np.nan
                #print(curr_rD, avg_num_GTP, avg_num_GDP_membrane, avg_num_Scd2GTP_membrane)
                
                #x_values.append(curr_rD)
                #y_values.append(avg_num_GTP)
                
                #break
                
                #file_count += 1
        except (FileNotFoundError):
            print("Couldn't find: ", filepath, " assuming zero")

fig, ax = plt.subplots(figsize=(6,6))

ax.set_xscale('log')
ax.set_yscale('log')

im = plt.pcolormesh(DD_mesh, rD_mesh, C, cmap=plt.cm.Blues, vmin=1.1, vmax=4.0)
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=16) 
cbar.ax.set_title(r'Patch radius ($\mu$m)', size=18, x=0.75, y=1.01)


plt.rcParams['hatch.linewidth'] = 2.0  # previous svg hatch linewidth

#now add hatches where the values are np.nan
for rD in rD_values:
    for DD in DD_values:
        curr_pos_rD = rD_values.index(rD)
        curr_pos_DD = DD_values.index(DD)
        
        if np.isnan(C[curr_pos_rD, curr_pos_DD]):
            bottom = rD_edge_values[curr_pos_rD]
            top = rD_edge_values[curr_pos_rD+1]
            
            left = DD_edge_values[curr_pos_DD]
            right = DD_edge_values[curr_pos_DD+1]
            
            print("Rectangle: ", left, bottom, right, top)
            ax.add_patch(Rectangle((left, bottom),(right-left), (top-bottom), facecolor='gray'))
        
        
color_list = [ (124.0/255.0, 124.0/255.0, 124.0/255.0), 'black', "#c882c8", "#a05abe", "#6e00be", '#003cc8', '#c8003c', '#82c8c8', 'purple', 'green', 'black', 'orange', 'green', 'red']

filepath = 'DepletionBoundary_R2Ridge.csv'

df = pd.read_csv(filepath)

column_count = 0
column_list = list(df)

for i in range(0, len(column_list)):
    column_split = column_list[i].split(',')

    if column_split[0] == 'D':
        curr_linestyle = '-'
        
        x_data = df[column_list[i]]
        y_data = df[column_list[i+1]]
        
        if column_split[1][0:6] == 'growth':
            if column_split[1][0:10] == 'growth-80%' or column_split[1][0:10] == 'growth-50%':
                curr_linestyle = '--'
                if column_split[1][0:10] == 'growth-80%' or column_split[1][0:10] == 'growth-50%':
                    plt.plot(x_data, y_data, label=column_split[1][7:] + " tip occupation", linestyle='-',color=color_list[column_count], linewidth=3.0)
                column_count += 1
        else:
            column_count += 1
        
        
        if column_count >= len(color_list):
            break

D_D_values = [0.00009, 0.0003, 0.001, 0.003333333, 0.011111111, 0.037037037]
koff_D_values = [0.0002, 0.000447214, 0.001, 0.002236068, 0.005, 0.01118034, 0.025, 0.055901699, 0.125]

color_array = []
color_array.append([0, 0, 0, 1, 1, 1, 1, 1, 1])
color_array.append([0, 0, 0, 1, 1, 1, 1, 1, 1])
color_array.append([0, 0, 1, 1, 1, 1, 1, 1, 1])
color_array.append([2, 1, 1, 1, 1, 1, 1, 1, 1])
color_array.append([2, 1, 1, 1, 1, 1, 1, 1, 1])
color_array.append([1, 1, 1, 1, 1, 1, 1, 1, 1])


#color_list = ['blue', 'red', 'green']


for i in range(0, len(koff_D_values)):
    for j in range(0, len(D_D_values)):
        curr_D_value = D_D_values[j]
        curr_koff_value = koff_D_values[i]
        
        curr_color_index = color_array[j][i]
        
        curr_color = color_list[curr_color_index]
        

file_list = ['YSM2468_averagedR2Grid_normalized.csv', 'averagedR2Grid_2xRitC_wYSM768_normalized.csv', 'averagedR2Grid_3xRitC_wYSM767_normalized.csv']
dashed_line_color = color_list[:]

color_order = ['Cdc42-1xRitC', 'Cdc42-2xRitC', 'Cdc42-3xRitC', 'Cdc42-WT (side)', 'Cdc42-WT (tip)']

color_count = 2

for filepath in file_list:
    print(filepath)
    
    X = []
    Y = []
    Z = []
    
    if False:
        pass
    else:
        with open(filepath) as fp:
            # skip first line
            line = fp.readline()
            line = fp.readline()

            tempX_list = []
            tempY_list = []
            tempZ_list = []
            
            priorY = -1.0
            
            while(line):
                split_string = line.split(',')

                if len(split_string) == 3:
                    tempX = float(split_string[0])
                    tempY = float(split_string[1])
                    tempZ = float(split_string[2])
                    
                    # if a new row add to X,Y,or Z and clear temp lists
                    if tempY != priorY and priorY >= 0.0:
                        X.append(tempX_list)
                        Y.append(tempY_list)
                        Z.append(tempZ_list)
                        
                        tempX_list = []
                        tempY_list = []
                        tempZ_list = []
                        
                    tempX_list.append(tempX)
                    tempY_list.append(tempY)
                    tempZ_list.append(tempZ)
                    
                    priorY = tempY
                line = fp.readline()
                
        intermediate_z_list = list(map(max,Z))
        max_val = max(intermediate_z_list)
        min_val = min(intermediate_z_list)
        
        for i in range(0, len(Z)):
            for j in range(0, len(Z[i])):
                if Z[i][j] == max_val:
                    print("MAX: ", filepath, X[i][j], Y[i][j])

                    break
        
        print("colors: ", color_count, color_list[color_count])
        CS = ax.contour(X, Y, Z,[max_val*0.96],linestyles="dashed", colors=[color_list[color_count]], vmin=min_val, vmax=max_val, zorder=10.0, linewidths=3.0)
    
    color_count += 1
      
                

ax.set_xlabel(r'$D_D$ $(\mu m^2/s)$', size=22)
ax.set_ylabel(r'$r_{D}$ $(s^{-1})$', size=22)

ax.tick_params(axis='both', which='major', labelsize=16)

plt.xlim(DD_edge_values[0], DD_edge_values[-1])
plt.ylim(rD_edge_values[1], rD_edge_values[-1])



plt.tight_layout(rect=[0, 0.0, 0.90, 1.0])

plt.show()