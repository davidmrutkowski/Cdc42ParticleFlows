# Cdc42ParticleFlows
Brownian dynamics code that simulates the diffusion and membrane binding/unbinding of proteins involved in S. pombe Cdc42 polarization under the effect of membrane flow. Particles are restricted to be in the cytoplasm or on the surface of a sphere or half spherocylinder. This code was developed for the simulations in Rutkowski, Vincenzetti, Vavylonis, Martin, bioRxiv, https://doi.org/10.1101/2023.07.21.550042

Compilation using g++: g++ -O3 -o Cdc42Diffusion Cdc42Diffusion.cpp

Command line arguments ./Cdc42Diffusion r_D D_D PriorParticleState.xyz

Can get last frame of output xyz for use as riorParticleState.xyz using python script GetLastFrameXYZ_folder.py
