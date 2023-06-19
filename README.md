# QuantumBreaking

This code reproduces the data used in 

M. Michel and S. Zell, *The Timescales of Quantum Breaking*, [arXiv::2306.09410](https://arxiv.org/abs/2306.09410)

We refer to this paper for a detailed description. 

The generated data can be found on zenodo [10.5281/zenodo.8039017](https://doi.org/10.5281/zenodo.8039017).

# Project Structure

<pre>
TimeEvolver                     TimeEvolver software package included as submodule
scripts                         bash scripts for running the programm for various parameter choices
CMakeLists.txt                  file needed for cmake
quantumBreaking.cpp             main function
quantumBreakingHamiltonians.*   defining the Hamiltonian of the system
</pre>

# External Requirements

Note that this project requires following external libraries
* Intel Math Kernel Library (MKL) 
* BOOST 
* HDF5 


# Compilation 
```
mkdir build; cd build     # create and use a build directroy
cmake ..                  # configuration
cmake --build .           # compilation and linking (or type "make")
```

Important: the submodule TimeEvolver is not pulled automatically. If it's the first time you need to use ``--init`` first. 
```
git submodule update --init --recursive
```

For problems during compilation, e.g. related to the version of BOOST, we refer to the repository of the [`TimeEvolver`](https://github.com/marco-michel/TimeEvolver). 
