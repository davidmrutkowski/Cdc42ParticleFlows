import math
import numpy as np
import matplotlib.pyplot as plt

initial_radius = 1.8

filename = "WT_rD_0.17000_DD_0.10000_singleNPPatch.xyz"
#filename = "1ritC_rD_0.02500_DD_0.00330_singleNPPatch.xyz"

arclength_hist = []
arclength_width = 0.1
num_bins = int(math.ceil(math.pi*initial_radius / arclength_width))

for i in range(0, num_bins):
    arclength_hist.append([])

with open(filename) as fp:
    line = fp.readline()
    
    split_string = line.split()
    
    num_particles = int(split_string[0])
    
    line = fp.readline()
    
    for i in range(0, num_particles):
        line = fp.readline()
        
        split_string = line.split()
        
        particle_type = int(split_string[0])
        
        if particle_type == 6:
            # this is a flow particle
            
            pos = np.array([float(split_string[2]), float(split_string[3]), float(split_string[4])])
            
            magnitude = np.linalg.norm(pos)
            
            initial_phi = math.atan2(pos[1], pos[0]);
            initial_theta = math.acos(pos[2] / magnitude);
                    
            flow_vector = np.array([float(split_string[5]), float(split_string[6]), float(split_string[7])])
            
            del_theta = 0.001
            
            unit_normal_z_m = np.array([initial_radius*math.cos(initial_phi)*math.sin(initial_theta-del_theta), initial_radius*math.sin(initial_phi)*math.sin(initial_theta-del_theta), initial_radius*math.cos(initial_theta-del_theta)])
            
            unit_normal_z_p = np.array([initial_radius*math.cos(initial_phi)*math.sin(initial_theta+del_theta), initial_radius*math.sin(initial_phi)*math.sin(initial_theta+del_theta), initial_radius*math.cos(initial_theta+del_theta)])
            
            unit_normal_z = (unit_normal_z_p - unit_normal_z_m) / np.linalg.norm(unit_normal_z_p - unit_normal_z_m)
            
            dot_prod = np.dot(unit_normal_z, flow_vector)
            
            #print(unit_normal_z, flow_vector)
            
            arclength = initial_theta * initial_radius
            
            tmp_bin = int(arclength/ arclength_width)
            
            arclength_hist[tmp_bin].append(dot_prod)
            
            #plt.plot(arclength, dot_prod, marker='o')
 
# these are per event in their raw state, need to calculate average value of del_t from kinetic monte carlo
q = 40.88014155/60.0 + 146.1504417/60.0
avg_t = 1.0 / q * (1.0)
print("avg_del_t: ", avg_t)

ax = plt.gca()

x_values = [0.0] + [(i+0.5)*arclength_width for i in range(0, num_bins)]
y_values = [0.0] + [np.mean(elem) / avg_t for elem in arclength_hist]

plt.errorbar(x_values, y_values, yerr=[0.0] + [np.std(elem) for elem in arclength_hist], color='red')
plt.plot([0.0, math.pi*initial_radius], [0.0, 0.0], color='black', linestyle='--')

#https://stackoverflow.com/questions/42995027/piecewise-linear-interpolation-function-in-python
linear_fit_values = [np.polyfit(x_values[i:i+2], y_values[i:i+2],1) for i in range(len(x_values)-1)]

for i in range(0, len(x_values)-1):
    if i == 0:
        print("if(x >= {0} && x < {1})\n{{\n\tdisp = {2}*x + {3};\n}}\n".format(x_values[i], x_values[i+1], linear_fit_values[i][0], linear_fit_values[i][1]), end="")
    elif i < len(x_values)-2:
        print("else if(x >= {0} && x < {1})\n{{\n\tdisp = {2}*x + {3};\n}}\n".format(x_values[i], x_values[i+1], linear_fit_values[i][0], linear_fit_values[i][1]), end="")
print("else\n{\n\tdisp = 0.0;\n}\n", end="")

print()

# or switch case statement
print("double xval = (x - {0})/{1};\n".format(0.5*arclength_width, arclength_width), end="")
print("int xswitch = (int)xval;\n", end="")
print("if(xval < 0.0)\n{\n\t xswitch = -1;\n}\n", end="")


for i in range(0, len(x_values)-1):
    if i == 0:
        print("switch(xswitch) {{\n\tcase {0}:\n\t\tdisp = {1}*x + {2};\n\t\tbreak;\n".format(i-1, linear_fit_values[i][0], linear_fit_values[i][1]), end="")
    elif i < len(x_values)-2:
        print("\tcase {0}:\n\t\tdisp = {1}*x + {2};\n\t\tbreak;\n".format(i-1, linear_fit_values[i][0], linear_fit_values[i][1]), end="")
print("\tdefault\n\t\tdisp = 0.0;\n}", end="")

ax.tick_params(axis='both', which='major', labelsize=16)

plt.xlim(0.0, math.pi*initial_radius)
plt.ylim(-0.005, 0.005)

plt.xlabel('arclength ($\mu m$)', size=18)
plt.ylabel('flow velocity ($\mu m / s$)', size=18)
plt.tight_layout()
plt.show()