# Cdc42ParticleFlows
Brownian dynamics code that simulates the diffusion and membrane binding/unbinding of proteins involved in S. pombe Cdc42 polarization under the effect of membrane flow. Particles are restricted to be in the cytoplasm or on the surface of a sphere or half spherocylinder. This code was developed for the simulations in Rutkowski, Vincenzetti, Vavylonis, Martin, Nature Communications, https://doi.org/10.1038/s41467-024-52655-1

Compiled on RHEL 8 with the C++17 standard using g++: g++ -std=c++17 -O3 -o Cdc42Diffusion Cdc42Diffusion.cpp

Expected compilation time on a typical Desktop computer is on the order of a minute.

Command line arguments ./Cdc42Diffusion r_D D_D PriorParticleState.xyz

Output from code includes an xyz file of particle positions every snapshotTime (default of 10 s), positions of exo/endocytosis events, number of particles of each state, and ratio of particles in the tip region of the sphere/spherocylinder shape (z >= 0.0) to the particles in the rest of the shape.

Can get last frame of output xyz for use as PriorParticleState.xyz using python script GetLastFrameXYZ_folder.py

Sphere version of the code is included in the sphere folder.

Spherocylinder version of the code is included in the spherocylinder folder, along with a subfolder which contains code and results that were used to determine hard-coded values for the function avgDisplacementAtArclength, a function that calculates the displacement of a particle at a given arclength due to the average flow of exocyotsis/endocytosis determined from sphere simulations.

Exepected run time depends on simulation time desired, but for default values can take around 16.7 hours for the sphere code, while the spherocylinder code takes around 78 hours on a 14600k.

Visualization of the xyz files can be done with the program Ovito using the attached Ovito state files. These state files will initially attempt to load a specific xyz, but after loading the state file, the correct xyz can be loaded.
