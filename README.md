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
./DEP-GAP -g [the number of your measurements] -k [the number of lithology type] -b [the number of layers] -x [the range of spatial location 1] -y [the range of spatial location 2] -d [your data directory] -l [initial lithology type] -o[the output directory] -t[the number of MCMC iterations] -s[seed] -c[the number of cores for parallel]
```
where 
   - *-g* for the number of physical measurements;
   - *-k* for the lithology type number;
   - *-b* for the layer number;
   - *-x* for the range of spatial coordinate 1;
   - *-d* for the range of spatial coordinate 2;
   - *-l* for the initial lithology type label;
   - *-o* for the directory of output;
   - *-t* for the number of iterations;
   - *-c* for the number of cores for parallel computing;

After running MCMC algorithm, your output directory stores the posterior sampling of MCMC algorithm. 

# Posterior analysis
## Posterior means
After saving the MCMC records in the output directory, the posterior estimates of parameters can be computed and the first half iterations will be regraded as burn-in.

```
cd /your/working/directory/
./posterior_means.R -i [MCMC record directory] -o [the posterior estimates directory]
```
