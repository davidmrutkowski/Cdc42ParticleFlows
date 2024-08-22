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
    This code calculates the surface concentration profile from the largest Cdc42-cluster
    for the spherocylinder as in e.g. Fig S2C,D for all xyz files in a
    directory given certain criteria of the file names.
    Select cancel in directory dialogue box to end adding folders to calculate the profiles on.
    Output: saves images in both svg and png format as well as a csv of the profiles
    in the folder selected.
""" 

import math
import tkinter as tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
import numpy as np

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
    

def calcDistanceSphere(posi, posj, sphereRadius):
    currDistance = calc_dist(posi, posj)
    
    asin_argument = 0.5* currDistance / sphereRadius;
    
    # this needs to be here for rounding error issues, if it is not, then some points are included as exo/endo that shouldnt be
    if asin_argument > 1.0 and asin_argument < 1.01:
        asin_argument = 1.0
    
    arclength = 2.0*sphereRadius * math.asin(asin_argument);
    
    return arclength

def calcDistanceCylinder(posi, posj, radius):
    posi_z = posi[2]
    posj_z = posj[2]
    
    posi_zeroz = [posi[0], posi[1], 0.0]
    posj_zeroz = [posj[0], posj[1], 0.0]
    
    currDistance_xy = calc_dist(posi_zeroz, posj_zeroz)
    
    asin_argument = 0.5* currDistance_xy / radius
    
    # this needs to be here for rounding error issues, if it is not, then some points are included as exo/endo that shouldnt be
    if asin_argument > 1.0 and asin_argument < 1.01:
        asin_argument = 1.0
    
    arclength_xy = 2.0 * radius * math.asin(asin_argument)
    
    length_z = abs(posi_z - posj_z)
    
    return math.sqrt(arclength_xy*arclength_xy + length_z*length_z)
    
def calcArcLengthDistance(posi, posj, sphereRadius):
    if posi[2] > 0.0 and posj[2] > 0.0:
        # both points are on the sphere
        return calcDistanceSphere(posi, posj, sphereRadius)
        
    elif posi[2] <= 0.0 and posj[2] <= 0.0:
        # both points are on the cylinder
        return calcDistanceCylinder(posi, posj, sphereRadius)
    else:
        # one of the points is on the sphere and one is on the cylinder
        if posi[2] > 0.0:
            # project this point onto the cylinder
            initialPhi = math.atan2(posi[1], posi[0])
            
            posi = np.array([sphereRadius*math.cos(initialPhi), sphereRadius*math.sin(initialPhi), posi[2]]);
        elif posj[2] > 0.0:
            # project this point onto the cylinder
            initialPhi = math.atan2(posj[1], posj[0])
            
            posj = np.array([sphereRadius*math.cos(initialPhi), sphereRadius*math.sin(initialPhi), posj[2]]);
        
        return calcDistanceCylinder(posi, posj, sphereRadius)   




    
def calcArcLengthFromNP(pos, sphereRadius):
    if pos[2] > 0.0:
        # then this point is on the sphere
        theta = math.acos(pos[2] / sphereRadius);
        
        return theta*sphereRadius
    else:
        # then this point is on the cylinder
        return math.pi*0.5*sphereRadius - pos[2]
        

plt.rcParams["font.family"] = "Arial"
    
initial_radius = 1.8
cylinder_length = 3.0*2.0*2.0

exclusion_center = [0.0, 0.0, initial_radius]

root = tk.Tk()
root.withdraw()

folder_selected = filedialog.askdirectory()

folderList = []

run_name = ""

fileList = []
while folder_selected != '':
    if len(fileList) == 0:
        #https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
        directory = os.fsencode(folder_selected)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(run_name + ".xyz") and "2d" not in filename and "frame" not in filename and "Run12" in filename:
                fileList.append(filename)
            else:
                 continue

        print(fileList)
   
    folderList.append(folder_selected)
        
    folder_selected = filedialog.askdirectory()

root.destroy()


combineFolderNames = ''

for tmp_folder in folderList:
    pos_Run_folder = tmp_folder.rfind(run_name)
    
    combineFolderNames += str(tmp_folder[pos_Run_folder:pos_Run_folder+4])
combineFolderNames += "_diffFont_"


for filename in fileList:
    plot_lable_list = ['Cdc42-GDP', 'Cdc42-GTP', 'Cdc42-GDP(cyto)', 'Scd2/1', 'Scd2/1(cyto)', 'Scd2/1/Cdc42-GTP', 'GAP', 'GAP(cyto)', 'All Cdc42', 'WT Cdc42-GTP', 'WT Cdc42-GDP(cyto)', 'Exocytosis', 'Endocytosis']

    colorList = ['#0072b2ff', '#d55e00ff', 'gray', '#cc79a7ff', 'gray', '#f0e442ff', '#000000ff', 'gray', '#009E73ff', '#d55e0070', 'gray', 'red', 'blue']
    
    linestyleList = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
    
    outputList = []


    hist_width = 0.01
    
    num_bins = int(math.ceil((math.pi * initial_radius*0.5 + cylinder_length) / hist_width))
    hist = []
    

    frame_count = -1
    frame_cutoff = 400
    max_frame = 1000

    avgNumBeads = []
    
    
    for folder_count in range(0, len(folderList)):
        tmp_folder = folderList[folder_count]
        
        realfilename = filename[:]
        if folder_count > 0:
            # need to rename filename
            pos_last_underscore = filename.rfind("_")
            
            pos_Run_folder = tmp_folder.rfind("Run")
            
            appendString = tmp_folder[pos_Run_folder:pos_Run_folder+4]
            
            realfilename = filename[:pos_last_underscore+1] + appendString + ".xyz"
            
        print("filename: ", realfilename)
        
        
        
        filepath = tmp_folder + "/" + realfilename
        
        print('Opening: ' + filepath)
        
        with open(filepath) as fp:
            maxCdc42Pos_list = []
            
            line = fp.readline()
            
            while line:
                split_string = line.split()
                
                if len(split_string) == 1:
                    frame_count += 1
                    
                    num_particles = 0
                    try:
                        num_particles = int(split_string[0])
                    except ValueError:
                        print("Value error: ", line)
                        break

                    
                    line = fp.readline()
                    
                    particles = []
                    
                    for n in range(0, num_particles):
                        line = fp.readline()
                        
                        if line == False:
                            print("no line")
                            break
                            
                        split_string = line.split()
                        
                        if len(split_string) >= 4:
                            try:
                                temp_pos = [float(split_string[2]), float(split_string[3]), float(split_string[4])]
                                
                                if(not math.isnan(temp_pos[0])):
                                    temp_type = int(split_string[0])
                                    
                                    while len(particles) <= temp_type:
                                        particles.append([])
                                    
                                    while len(hist) <= temp_type:
                                        hist.append([0] * num_bins)
                                        avgNumBeads.append(0)
                                        
                                    particles[temp_type].append(temp_pos)
                            
                            except ValueError:
                                print("Skipping line: ", line)
                                print("Exiting file")
                                break
                            
                    if frame_count > max_frame:
                        break
                        
                    #1 is RitC Cdc42-GTP, 9 is WT Cdc42-GTP        
                    if frame_count > frame_cutoff and ((len(particles) > 1 and len(particles[1]) >= 1) or (len(particles) > 9 and len(particles[9]) >= 1)):
                        GTP_particle_list = particles[1] + particles[9]
                        
                        # find GTP cluster com
                        cluster_dist = 0.30
                        cluster_id = [-1]*len(GTP_particle_list)
                        current_cluster_ids = []
                        
                        potential_exo_sites = []
                        cluster_sizes = []
                        
                        max_cluster_size = 0
                        
                        curr_cluster = 0
                        
                        for p in range(0, len(GTP_particle_list)):
                            # cluster id of this particle not yet determined
                            if cluster_id[p] == -1:
                                current_cluster_ids.append(p);
                                
                                curr_cluster_size = 0;
                                cluster_centerofmass = [0.0, 0.0, 0.0];
                                
                                cluster_id[p] = curr_cluster
                                
                                while len(current_cluster_ids) > 0:
                                    curr_index = current_cluster_ids.pop()
                                    
                                    
                                    pos_p = GTP_particle_list[curr_index]
                                    
                                    cluster_centerofmass[0] = cluster_centerofmass[0] + pos_p[0]
                                    cluster_centerofmass[1] = cluster_centerofmass[1] + pos_p[1]
                                    cluster_centerofmass[2] = cluster_centerofmass[2] + pos_p[2]
                                    
                                    curr_cluster_size += 1
                                    
                                    
                                    
                                    for q in range(0, len(GTP_particle_list)):
                                        if cluster_id[q] == -1:
                                            pos_q = GTP_particle_list[q]
                                            
                                            arc_distance = calcArcLengthDistance(pos_p, pos_q, initial_radius)
                                            
                                            if arc_distance < cluster_dist:
                                                # curr_index and q are part of the same cluster
                                                cluster_id[q] = curr_cluster
                                                current_cluster_ids.append(q)
                                        
                                
                                if curr_cluster_size > max_cluster_size:
                                    #print(curr_cluster_size, max_cluster_size)
                                    max_cluster_size = curr_cluster_size
                                    cluster_centerofmass[0] = cluster_centerofmass[0] / float(curr_cluster_size)
                                    cluster_centerofmass[1] = cluster_centerofmass[1] / float(curr_cluster_size)
                                    cluster_centerofmass[2] = cluster_centerofmass[2] / float(curr_cluster_size)
                                    
                                    # project this position back onto sphere surface (center of mass will be off sphere in most cases)
                                    cluster_com_unit = get_unit_vec(cluster_centerofmass)
                                    
                                    potential_exo_sites.append([val * initial_radius for val in cluster_com_unit])
                                    cluster_sizes.append(curr_cluster_size)
                                
                                curr_cluster += 1
                        
                        # need to find max Cdc42-GTP concentration point and rotate exclusionCenter based around this point (north pole -> max Cdc42-GTP point)            
                        maxCdc42Pos = [0.0, 0.0, initial_radius]
                        if len(potential_exo_sites) > 0:
                            maxCdc42Pos = potential_exo_sites[-1]
                            
                        maxCdc42Pos_list.append(maxCdc42Pos)
                        
                        
                        #this is not based on the same amount of time unless the time between frames is the same!
                        if len(maxCdc42Pos_list) >= 5:
                            
                            avgCdc42Pos = [np.average([maxCdc42Pos_list[i][j] for i in range(0, len(maxCdc42Pos_list))]) for j in range(0, 3)]
                            
                            #need to project average position onto sphere surface
                            curr_avg_mag = calc_mag(avgCdc42Pos)
                            avgCdc42Pos = [val / curr_avg_mag * initial_radius for val in avgCdc42Pos]
                            
                            #print(maxCdc42Pos[0], maxCdc42Pos[1], maxCdc42Pos[2])
                            
                            for h in range(0, len(particles)):
                                for p in range(0, len(particles[h])):                                   
                                    temp_arclength = calcArcLengthFromNP(particles[h][p], initial_radius)

                                    hist[h][int((temp_arclength / hist_width))] += 1
                                    
                                    avgNumBeads[h] += 1
                                    
                            maxCdc42Pos_list.pop(0)
                            
                        
                line = fp.readline()
            
       
                        
    for h in range(0, len(hist)):
        if avgNumBeads[h] > 0:
            print("avgNumH: ", h, avgNumBeads[h])
            curr_avgNumBeads = avgNumBeads[h] / float(frame_count - frame_cutoff)
            #avgNumBeads = 30000.0
            idealNumDensity = curr_avgNumBeads / (4.0*math.pi*initial_radius**2)
            
            tempX = []
            tempY = []
            for i in range(0, len(hist[h])):
                inner_arclength = i*hist_width
                outer_arclength = (i+1)*hist_width
                
                inner_area = 0.0
                outer_area = 0.0
                
                if inner_arclength < math.pi*initial_radius*0.5 and outer_arclength < math.pi*initial_radius*0.5:
                    #both are on hemisphere
                    inner_area = calc_sphericalSA(inner_arclength, initial_radius)
                    outer_area = calc_sphericalSA(outer_arclength, initial_radius)
                elif inner_arclength >= math.pi*initial_radius*0.5 and outer_arclength >= math.pi*initial_radius*0.5:
                    #both are on cylinder
                    inner_area = 2.0*math.pi*initial_radius*initial_radius + (inner_arclength - math.pi*initial_radius*0.5)*(2.0*math.pi*initial_radius)
                    outer_area = 2.0*math.pi*initial_radius*initial_radius + (outer_arclength - math.pi*initial_radius*0.5)*(2.0*math.pi*initial_radius)
                else:
                    #inner is on sphere, outer is on cylinder
                    inner_area = calc_sphericalSA(inner_arclength, initial_radius)
                    outer_area = 2.0*math.pi*initial_radius*initial_radius + (outer_arclength - math.pi*initial_radius*0.5)*(2.0*math.pi*initial_radius)
               
                
                idealNumInRegion = idealNumDensity * (outer_area - inner_area)
                area_ribbon = outer_area - inner_area
                
                #print(i, area_ribbon)
                
                tempX.append(hist_width*(i+0.5))
                tempY.append(hist[h][i] / float(frame_count - frame_cutoff) / area_ribbon)
                #print('{0} {1} {2} {3}'.format(hist_width*(i+0.5), hist[i] / float(frame_count - 1 - frame_cutoff) / idealNumInRegion, idealNumInRegion, outer_arclength - inner_arclength))
                
            if h != 10:
                if h == 6 or h == 7:
                    outputList.append([sortBy(filename), tempX, [1*elem for elem in tempY], h])
                else:
                    outputList.append([sortBy(filename), tempX, tempY, h])
     outputList.append([sortBy(filename), outputList[0][1], [outputList[0][2][i]+outputList[1][2][i]+outputList[3][2][i] for i in range(0, len(outputList[0][2]))], 8])

    styleList = ['solid', 'solid', 'dotted']
    alphaList = [1.0, 0.5]

    maxY = 3.6

    first_str = 'rD_'
    posLastUnderscore = filename.rfind(first_str)
    posLastDot = filename.rfind('.xyz')


    output_txt_name = folderList[-1] + '/RDist_' + filename[posLastUnderscore+len(first_str):posLastDot] + combineFolderNames + 'frameCutoff' + str(frame_cutoff) + '.csv'
    print("outputting to: ", output_txt_name)
    
    with open(output_txt_name, 'w') as fp:
        #output to csv file
        fp.write('{0},'.format('Arclength from centroid'))
        for i in range(0, len(outputList)):
            fp.write('{0},'.format(plot_lable_list[int(outputList[i][3])]))
        fp.write('\n')
        
        for j in range(0, len(outputList[0][1])):
            fp.write('{0},'.format(outputList[0][1][j]))
            
            for i in range(0, len(outputList)):
                tmp_index = int(outputList[i][3])
                
                #print(i, tmp_index)
                #if tmp_index < 10:
                fp.write('{0},'.format(outputList[i][2][j]))
                    
            fp.write('\n')
                
            
    
    for i in range(0, len(outputList)):
        tmp_index = int(outputList[i][3])

        tmp_label = plot_lable_list[tmp_index]
        
        if i//len(colorList) == 0:
            print(i, tmp_index, colorList[tmp_index])
            
            plt.plot(outputList[i][1], outputList[i][2], color=colorList[tmp_index], linestyle = linestyleList[tmp_index], linewidth=2.0, label=tmp_label)
        else:
            plt.plot(outputList[i][1], outputList[i][2], color=colorList[tmp_index], linestyle = linestyleList[tmp_index], linewidth=2.0, label='')
        
    plt.xlim(0, (initial_radius*math.pi)*0.99)
    plt.ylim(0, 200)
    plt.xlabel(r'Arclength from centroid [$\mu$m]', fontsize=26)
    plt.ylabel(r'Surface concentration [#/$\mu$m$^2$]', fontsize=26)

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=22)
    plt.xticks(np.arange(0, 6, 1.0))
    plt.legend()
    plt.tight_layout()

    print("saving images: ", folderList[-1] + '/RDist_' + filename[posLastUnderscore+len(first_str):posLastDot] + combineFolderNames + 'frameCutoff' + str(frame_cutoff) + '.png')
    
    plt.savefig(folderList[-1] + '/RDist_' + filename[posLastUnderscore+len(first_str):posLastDot] + combineFolderNames + '_avgCentroid_frameCutoff' + str(frame_cutoff) + '.png')
    plt.savefig(folderList[-1] + '/RDist_' + filename[posLastUnderscore+len(first_str):posLastDot] + combineFolderNames + '_avgCentroid_frameCutoff' + str(frame_cutoff) + '.svg')

    plt.clf()