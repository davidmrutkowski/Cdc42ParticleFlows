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
    Code that takes 3d points on sphere surface and projects them to 2d
    using Mollweide projection. Automatically places the largest Cdc42-GTP
    patch of the last frame at the center of the projection.
"""

import os
import tkinter as tk
from tkinter import filedialog
import math
import numpy as np

#https://stackoverflow.com/questions/2301789/how-to-read-a-file-in-reverse-order
def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename, 'rb') as fh:
        segment = None
        offset = 0
        file_size = fh.seek(0, os.SEEK_END)
        
        fh.seek(file_size - offset)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size)).decode(encoding='utf-8')
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

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

def calc_dist(list1, list2):
    temp_x = list1[0] - list2[0]
    temp_y = list1[1] - list2[1]
    temp_z = list1[2] - list2[2]
    
    return math.sqrt(temp_x*temp_x + temp_y*temp_y + temp_z*temp_z)
    
            
            
path_list = []

root = tk.Tk()
root.withdraw()

path = filedialog.askopenfilename(title = "Select xyz file",filetypes = (("XYZ Files","*.xyz"),))

root.destroy()

output_file_name = path[:-4] + '_2dProj_alt.xyz'
output_file = open(output_file_name, 'w')

count = 0
lambda_0 = 0.0

rotation_angle = 0.0

initial_radius = 1.8


GTP_particle_list = []

#find last position of nucleator bead
file_gen_reverse = reverse_readline(path)
for row in file_gen_reverse:
    split_string = row.split()
        
    if(len(split_string) >= 5):
        particle_type = int(split_string[0])
        particle_tag = int(split_string[1])
        pos_x = float(split_string[2])
        pos_y = float(split_string[3])
        pos_z = float(split_string[4])
        
        if particle_type == 1 or particle_type == 9:
            GTP_particle_list.append([pos_x, pos_y, pos_z])
            
    elif len(split_string) == 1:
        break




#GTP_particle_list = particles[1] + particles[9]

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
    
                            
                            


a_vec = np.array(maxCdc42Pos)
a_vec = a_vec / np.linalg.norm(a_vec)


# this is vector that a is rotated to
b_vec = np.array([1.0, 0.0, 0.0])

v_vec = np.cross(a_vec, b_vec)

c = np.dot(a_vec, b_vec)

# now find rotation matrix according to https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
vx_matrix = np.array([[0.0, -v_vec[2], v_vec[1]], [v_vec[2], 0.0, -v_vec[0]], [-v_vec[1], v_vec[0], 0.0]])

R = np.identity(3) + vx_matrix + np.dot(vx_matrix, vx_matrix) * (1.0 / (1.0 + c))



print(R)


with open(path) as fp:
    line = fp.readline()
    
    while line:
        split_string = line.split()
        #print(line)
        
        if len(split_string) >= 5 and split_string[0] != '10':
            x = float(split_string[2])
            y = float(split_string[3])
            z = float(split_string[4])
            
            rotated_pos = np.dot(R, np.array([x, y, z]))

            
            x = rotated_pos[0]
            y = rotated_pos[1]
            z = rotated_pos[2]
            
            #https://gis.stackexchange.com/questions/120679/equations-to-convert-from-global-cartesian-coordinates-to-geographic-coordinates
            r = math.sqrt(x*x + y*y + z*z)
            #print(line)
            
            latitude = math.asin(z / r)
            longitude = math.atan2(y, x) + rotation_angle
            
            #print(longitude, longitude + math.pi/2.0)
            while longitude > math.pi:
                #print(longitude*180.0/math.pi)
                longitude = -math.pi + (longitude - math.pi)
                #print(longitude*180.0/math.pi)
            while longitude < -math.pi:
                print('hello')
                longitude = math.pi - (longitude + math.pi)
            
            
            #https://en.wikipedia.org/wiki/Transverse_Mercator_projection
            #a = 1.8
            # x_2d = a * longitude
            #y_2d = a*0.5 * math.log10((1.0 + math.sin(latitude)) / (1.0 - math.sin(latitude)))
            
            #https://en.wikipedia.org/wiki/Sinusoidal_projection
            #x_2d = (longitude - lambda_0) * math.cos(latitude)
            #y_2d = latitude
            
            #https://en.wikipedia.org/wiki/Mollweide_projection
            theta = latitude
            if abs(latitude) != math.pi / 2.0:
                old_theta = theta
                theta = theta - (2.0 * theta + math.sin(2.0*theta) - math.pi*math.sin(latitude)) / (2.0 + 2.0*math.cos(2.0*theta))
                
                count = 0
                while abs(theta - old_theta) > 0.01:
                    old_theta = theta
                    theta = theta - (2.0 * theta + math.sin(2.0*theta) - math.pi*math.sin(latitude)) / (2.0 + 2.0*math.cos(2.0*theta))
                    count += 1
                    
                #print(count)
                    
            x_2d = 1.8 * 2.0 * math.sqrt(2.0) / math.pi * (longitude - lambda_0) * math.cos(theta)
            y_2d = 1.8 * math.sqrt(2.0) * math.sin(theta)
            
            #https://en.wikipedia.org/wiki/Mercator_projection
            #x_2d = 1.8 * (longitude - lambda_0)
            #y_2d = 1.8 * math.log(math.tan(math.pi/4.0 + latitude/2.0))
            
            #print(x,y,z, latitude, longitude)
            
            #if count >= 20:
            #    exit()
                
            count += 1
            
            output_file.write('{0} {1} {2} {3} {4}\n'.format(split_string[0], split_string[1], round(x_2d,3), round(y_2d,3), 0.0))
            
            
        else:
            output_file.write(line)
            
        line = fp.readline()