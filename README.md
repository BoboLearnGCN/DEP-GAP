# DEP-GAP: Depth effects correction and lithology type prediction for geophysical inversion data

Geophysical inversion is the mathematical process of predicting underground geophysical properties at different depths from wave signals detected on the ground, such as seismic waves. Three-dimensional reconstruction of lithology type based on geophysical inversion data enables determining the optimal drilling location and depth. However, the normalization procedures inherent in geophysical inversion often leads to substantial discrepancies between the inferred geophysical properties and those directly measured from drilled rock samples. As a result, geologists usually have to annotate inversion results manually based on expert knowledge. These discrepancies, which can vary with depth even for the same lithology type, will further hamper manual annotation. 

To address these challenges, we develop an unsupervised hierarchical Bayesian model, Depth-corrected Gaussian-Potts Mixture model, abbreviated as DEP-GAP. DEP-GAP integrates a location-and-scale adjustment model with a Potts spatial prior to simultaneously identify lithology types, adjust for depth-related effects, and encourage spatial smoothness.

# System Requirements

## OS Requirements

The C++ source codes support *Linux* operating systems. It has been tested on the following systems:

Linux: Ubuntu 18.04 and CentOS 7

# Installation

1. Installing RngStreams: Please download [RngStreams](http://statmath.wu.ac.at/software/RngStreams/) and install it by
```
wget http://statmath.wu.ac.at/software/RngStreams/rngstreams-1.0.1.tar.gz
tar zxvf rngstreams-1.0.1.tar.gz
cd rngstreams-1.0.1
./configure --prefix=<prefix_path> #prefix_path is your home directory/any directory you want to install at
make
make install
```
2. Compiling DEP-GAP: Please download all source code in *DEP-GAP/cpp* into your working directory and compile by
```
make
```

The compiling of C++ codes will only take a few seconds.


# Run DEP-GAP: 

For general cases, you can change the working directory and run the following two commands in the terminal to conduct MCMC sampling and posterior inference directly:
```
cd /your/working/directory/
/your/DEP-GAP/source/code/directory/DEP-GAP -g [the number of your measurements] -k [the number of lithology type] -b [the number of layers] -x [the range of spatial location 1] -y [the range of spatial location 2] -d [your data directory] -l [initial lithology type] -o[the output directory] -t[the number of iterations] -s[seed] -c[the number of cores for parallel]
```
where 
   - *-d* for the directory to save the posterior sampling and inference result;
   - *-r* for the directory of read count data and dimension information created in the step 3;
   - *-p* for the project name, which should be consistent with the name of read count data and dimension file;
   - *-v* for the version number;
   - *-K* for the total number of cell types;
   - *-i* for the number of iterations in the MCMC sampling;
   - *-o* for the number of iterations printed into the hard disk to control memory usage;
   - *-s* for the seed of MCMC sampling to let the results be reproducible;
   - *-c* for the number of cores for parallel computing;
   - *-b* for the number of burn-in iterations in the posterior inference.

After running MCMC algorithm, there are two folders created to store the results. `MCMC_sampling_K[the number of cell types]` stores the posterior sampling of MCMC algorithm, and `Inference_K[the number of cell types]` saves the posterior inference of all parameters. 
