# Gibbs Sampling for Sound Source Separation
This is the repo for the gibbs sampling on a MRF for sound source separation. This code will simulate the sound source separation running on both CPU and the hardware acclerator. For CPU mode, it will directly process the entire MRF without other preprocessing steps. For acclerator mode, it will pratition the entire MRF into different tiles. Each tile will be 24 nodes wide, which is supported for the 12 samplers in the hardware acclerator. Also it would apply an auto normalization on top of the original input to prevent any kind of overflow or underflow situations for the acclerator. 

# Three methods to run the gibbs sampling:
0: Sequential Gibbs Sampling for CPU.

1: Chromatic Gibbs Sampling for hardware acclerator, with fix point numbers.

2: Advanced Chromatic Gibbs Sampling for hardware acclerator, with fix point numbers and runtime normalization.    

# Standard commands to run the code:
For CPU: 
./gibbs 513 125 2 4 50 4 0.5 125 0 1 0

For hardware acclerator: 
./gibbs 512 125 2 4 50 4 0.5 24 1 1 0 

More details for each parameter are shown in the Makefile.

