# Cdc42ParticleFlows
Brownian dynamics code that simulates the diffusion and membrane binding/unbinding of proteins involved in S. pombe Cdc42 polarization under the effect of membrane flow. Particles are restricted to be in the cytoplasm or on the surface of a sphere or half spherocylinder. This code was developed for the simulations in Rutkowski, Vincenzetti, Vavylonis, Martin, bioRxiv, https://doi.org/10.1101/2023.07.21.550042

Compiled with the C++17 standard using g++: g++ -std=c++17 -O3 -o Cdc42Diffusion Cdc42Diffusion.cpp

Command line arguments ./Cdc42Diffusion r_D D_D PriorParticleState.xyz

Can get last frame of output xyz for use as PriorParticleState.xyz using python script GetLastFrameXYZ_folder.py

Sphere version of the code is included in the Sphere folder.

Spherocylinder version of the code is included in the Spherocylinder folder, along with a subfolder which contains code and results that were used to determine hard-coded values for the function avgDisplacementAtArclength, a function that calculates the displacement of a particle at a given arclength due to the average flow of exocyotsis/endocytosis determined from sphere simulations.

Visualization of the xyz files can be done with the program Ovito using the attached Ovito state files. These state files will initially attempt to load a specific xyz, but after loading the state file, the correct xyz can be loaded.
