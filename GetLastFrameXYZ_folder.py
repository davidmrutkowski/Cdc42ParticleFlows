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

#this code gets the last frame from every xyz in the folder foldername

import os
import tkinter as tk
from tkinter import filedialog

def periodicWrap(dist, boxl):
    return dist - boxl*round(dist / boxl)

root = tk.Tk()
root.withdraw()

foldername = filedialog.askdirectory()

root.destroy()

file_list = []

for file in os.listdir(foldername):
    filename = os.fsdecode(file)
    if filename.endswith(".xyz") and "2d" not in filename and "frame" not in filename:
        file_list.append(foldername + "/" + filename)
     


for filepath in file_list:
    with open(filepath) as fp:
        print(filepath)
        frameCount = -1

        numAtoms = 0

        output_str = ""
        
        new_output_str = ""

        line = fp.readline()
        
        while line:
            splitString = line.split()
            
            if(len(splitString) == 1):            
                
                
                numBeads = 0
                try:
                    numBeads = int(splitString[0])
                except ValueError:
                    break
                    
                frameCount += 1
                    
                splitString = line.split()
                num_particles = int(splitString[0])
                
                output_str = ""
                
                output_str += str(num_particles-1) + "\n"
                #skip blank line
                line = fp.readline()
                
                output_str += line
                #outputFile.write(line)
                
                abnormal_break = False
                for i in range(0, numBeads):
                    line = fp.readline()
                    splitString = line.split()
                    
                    if len(splitString) >= 5:
                        type = int(splitString[0])
                        
                        if type != 10:
                            outputLine = str(type) + " " + splitString[1] + " " + splitString[2] + " " + splitString[3] + " " + splitString[4] + "\n"
                            
                            output_str += outputLine
                            #outputFile.write(outputLine)
                    else:
                        abnormal_break = True
                        break
                
                if abnormal_break == False:
                    new_output_str = output_str[:]
            else:
                line = fp.readline()
    
    print(frameCount)
    outputFile = open(filepath[:-4] + '-frame' + str(frameCount) + '.xyz', 'w')
    outputFile.write(new_output_str)