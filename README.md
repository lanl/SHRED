# SHRED

Stochastic and Hybrid Representation Electronic structure by Density functional theory

© 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

# Is this project active?
You may notice that this version of SHRED appears inactive. Note that this version of SHRED constitutes a stable release and is meant primarily for standard use cases and is updated infrequently (next updated expected late 2025). An active development version, maintained by the LANL team, lives a private repository. Friendly users interested in access to latest capabilities and collaboraton should contact Alexander J. White at alwhite@lanl.gov to discuss access. 

# What is SHRED for?

SHRED is a plane-wave DFT code similar to ABINIT, VASP and Quantum Espresso. It utilizes mostly the same conventions as ABINIT and has been tested against ABINIT whenever possible. While VASP and ABINIT were originally designed to target condensed phase solids (elegantly using symmetry to simplify computation), SHRED is designed to tackle numerically difficult problems for large-scale disordered systems and/or high-temperature (dense plasmas) calculations. In addition to standard deterministic finite temperature Kohn-Sham DFT within the Mermin approximation, SHRED can perform stochastic and mixed-stochastic deterministic calculations. Stochastic DFT is a linear scaling in system size and inverse scaling in temperature approach to solving for the finite-temperature electron density within the Kohn-Sahm formalism. Deterministic KS calcualtions scale cubically with system size and temperature, albiet with a much lower prefactor. Mixed DFT combines both approaches to improve precision and/or reduce computational costs.  SHRED can also perform DFT calcualtions utilizing Orbital-Free (OF) DFT with Thomas-Fermi Perrot or Thomas-Fermi Perrot Von Weizsacker Kinetic Free-Energy (TFvW) functionals.

In addition to single point calculation and Born-Oppenheimer molecular dynamics, SHRED can perform linear-electron response (conductivity, optical response, thermal conductivity, thermopower) within the Kubo-Greenwood approximation, real-time TD-DFT and Nonadiabatic Ehrenfest TD-DFT and stopping power calculations. These can be performed for all flavors of DFT available in SHRED with both HGH and PAW pseudopotentials. For TD-OF-DFT SHRED can deploy a nonadiabtic (i.e. current dependent) kinetic energy functional developed by A. J. White, or utilize the TD-TFvW approximation.  

Exchange correlation functionals are handeled through libXC, and are currently limited to LDA and GGA functionals. 

# License
This software is open source under the GNU General Public License; either version 3.0 of the License, or (at your option) any later version. See the included GPLv3.pdf for the full license. 

# Contributors 
Main Author:

    Alexander J. White (Los Alamos National Laboratory, T-1 Group)

Contributing Authors:
    
    Legacy Code- Ondrej Certik (GSI Technology; Former, Los Alamos National Laboratory, T-1 Group)
    
Significant Testing: 

    Vidushi Sharma (Los Alamos National Laboratory, T-1 Group)
    Lee Collins (Los Alamos National Laboratory, T-1 Group)
    Katarina Nichols (Univ. of Rochester; Los Alamos National Laboratory, T-1 Group)


# Build
Some Makefiles are provided, the default is gcc with openmpi. 
These dependencies can be installed with Homebrew on MAC OS and  wget for unix (probably). The code uses some fortran 2018 so you need a gfortran 8 or newer.
Dependencies include gcc (if your gcc/gfortran is out of date), open-mpi, lapack, fftw, and libxc. These can all be installed on MAC OS with Homebrew. 
https://brew.sh/

    brew install gcc
    brew install open-mpi
    brew install lapack
    brew install fftw
    brew install libxc

The homebrew default symlinks work with the main Makefile without needing to set. So at this point it is

    make
or
    make gnu_debug

for productive or debug builds. You can overwrite (for example if you use wget) the default paths with your own paths via:

    make FFTW_DIR= your path XC_DIR= your path LA_DIR= your path

For the HPC example makefiles you need to set the paths either as above or by exporting

    export FFTW_DIR= your path 
    export XC_DIR= your path 
    export LA_DIR= your path
    make

While the primary MAC OS build is provided, SHRED is really intended as a HPC code. The primary Makefile can be used as a template for your specific system. Additional Makefiles for LANL HPC systems are included, providing examples for various compilers and compiler options. It builds and installs directly in the shred folder. At this time, configuration and/or CMAKE builds are not available. 

# Execution
The build creates an executable main.exe, so it runs with N ranks:

    mpirun -n N /YOUR PATH TO REPO/shred/main.exe
    
The code reads a single input file called "input". The tex and PDF input file manual for the input file is located in the Documentation folder. Pseduopotential files also need to be included, with the filenames specified in the input files. Hartwigsen-Goedeker-Hutter (HGH) pseudopotentials are the only Norm-conserving PP's. SHRED utilizes the same format as ABINIT, PP's can be downloaded at: https://sourceforge.net/p/cp2k/code/HEAD/tree/trunk/potentials/Goedecker/abinit/. Projector Augmented Wave in the "XML" format can be read by SHRED through libPAW. This format can be generated by AtomPAW. JTH and GBRV datasets can be downloaded in xml format at https://www.abinit.org/psp-tables and http://www.physics.rutgers.edu/gbrv/ respecively. 

# Libraries
In addition to core mathematical and computer science libraries (lapack, open-mpi, fftw). The code relies on two physics libraries. The first, LibPAW, calcualtes many of the PAW energy, and Hamiltonian terms. We have had to modify LibPAW from ABINIT's version, this modified version is included in the repository: 

    ABINIT: Overview, and focus on selected capabilities J. Chem. Phys. 152, 124102 (2020).

The LibXC library calculates the exhange-correlation terms from provided electron densities and gradients. This library is linked and can be updated using homebrew. 

    Recent developments in Libxc - A comprehensive library of functionals for density functional theory, Software X 7, 1 (2018)

# Readings
We assume a basic understanding of DFT in the plane wave pseudopotential methods. If not, there are too many good sources to choose from, simply Google "DFT planewaves". Otherwise, these are some of the main papers describing the methods and algorithims used in SHRED. 

PAW in LibPAW/SHRED language: 

    Implementation of the Projector Augmented-Wave Method in the ABINIT code. Comput. Mat. Science 42, 337 (2008).
    
Kubo Greenwood approximation for linear response of electrons:
    
    Electronic transport coefficients from ab initio simulations and application to dense liquid hydrogen Phys. Rev. B 83, 235120 (2011).

Orbital-Free DFT:

    Very-high-temperature molecular dynamics Phys. Rev. E 73, 016403 (2006)

Stochastic DFT (what is implemented here):

    Stochastic density functional theory at finite temperatures. Phys. Rev. B 97, 115207 (2018)
    Transition to metallization in warm dense helium-hydrogen mixtures using stochastic density functional theory within the Kubo-Greenwood formalism Phys. Rev. B 100, 195101(2019)

Mixed Stochastic Deterministic DFT (mDFT):

    Fast and Universal Kohn Sham Density Functional Theory Algorithm for Warm Dense Matter to Hot Dense Plasma  Physical Review Letters 125, 055002

Stopping Power Calculations with real-time TD-DFT, TD-mDFT and Time-Dependent Orbital Free DFT :

    Time-dependent orbital-free density functional theory for electronic stopping power: Comparison to the Mermin-Kohn-Sham theory at high temperatures Physical Review B 98 (14), 144302
    Mixed stochastic-deterministic time-dependent density functional theory: application to stopping power of warm dense carbon Journal of Physics: Condensed Matter 34 (17), 174001
    
Eigensolver & Density Mix Algorithims:

    Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set Phys. Rev. B 54, 11169 (1996)
    Parallel eigensolvers in plane-wave density functional theory Computer Physics Communications 187, 98-105 (2015) 

Direct minimization for orbital-free electron density:
  
    Conjugate-gradient optimization method for orbital-free density functional calculations J. Chem. Phys. 121, 2030 (2004);
    
