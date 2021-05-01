# Dissociation Dynamics of H2+ in XUV and IR laser fields

## Methods

The following methods should be used:

* [Finite difference method](https://en.wikipedia.org/wiki/Finite_difference_method)
* [Exact diagonalization](https://en.wikipedia.org/wiki/Exact_diagonalization)
* [Split-step Fourier method](https://en.wikipedia.org/wiki/Split-step_method)
* [Fourier analysis](https://en.wikipedia.org/wiki/Split-step_method)


## Data

The following data is provided in atomic units and stored in csv format:

* dipole_coupling.dat: (R-dependent) dipole matrix between relevant states
* H2nuclwf.dat: H2 ground state wave function (initial state)
* H2p_pot_gerade.dat: Born-Oppenheimer surface of one relevant state (binding)
* H2p_pot_ungerade.dat: Born-Oppenheimer surface of other relevant state (non-binding)

Hints:

1. Not all the data is given on the same spatial grid,
   so you may have to interpolate and extrapolate.
2. Also, the given grid spacing and size may not be optimal for the numerical
   integration that you do. Think carefully about your choice of grid size
   and spacing, also to not waste computational resources!
3. The dipole operator is non-diagonal in the internal states.
   How to handle this issue to still get good computational performance
   is the crux of the problem. 
   
Atomic Units:

1. Length in Bohr radii (r_0 = 0.529177 1e-10 m)
2. Energy given in Hartree energies (E_H = 27.211 eV)
3. For further (e.g. time) see [here](https://de.wikipedia.org/wiki/Atomare_Einheiten) or [here](https://en.wikipedia.org/wiki/Hartree_atomic_units#Units)


## Steps

### 1. Find the vibrational eigenstates in the H2+ ground state potential

1. ground state potential = H2p_pot_gerade.dat 
2. Born Oppenheimer approximation -> get H
3. exact diagonalization
4. verify := eigenstates CAN be similar to HO, dependent on how low-lying the eigenstates are and therefore how close the potential is to HO (For HO computations see experiments.ipynb)

Methods:

* Exact diagonalization
* Finite difference method


### 2. Simulate the wave packet propagation without the IR laser field

(1. compute H(t) using split-step fourier method)
(2. apply to H2nuclwf.dat)

1. decompose H2nuclwf.dat into vibational eigenstates
2. compute time-evolution of eigenstates using split-step fourier method

Methods:

* Split-step fourier
* fourier analysis


### 3. Simulate the dynamics of the time-dependent system
Simulate the dynamics of the time-dependent system with two coupled energy surfaces,
coupled by the time-dependent IR laser field. Visualize the wave packet dynamics.

1. determine time-dependent barrier "height"
2. compute effect on the time-evolution of the vibrational eigenstates (Do they: Escape/ Escape partially/ Nothing at all?)


### 4. Scan the delay time and analyze the momentum distribution
Scan the delay time and analyze the momentum distribution of the dissociation products.

1. See time dependent-escape effects (Where is the connection between escape and break-up?)
2. ?

## References

* [Molecular Dissociative Ionization and Wave-Packet Dynamics Studied Using Two-Color XUV and IR Pump-Probe Spectroscopy, Kelkensberg et al. 2009](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.103.123005)
* [Molecular wave-packet dynamics on laser-controlled transition states, Fischer et al. 2016](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.93.012507)
* [Avoided Crossing for understanding of Dipole-Coupling](https://www.univie.ac.at/columbus/workshops/Tianjin2016/Lectures/PLASSER/FP2_handout.pdf)
* [Project Description](https://uebungen.physik.uni-heidelberg.de/uebungen/download/4233/Programming%20projects.pdf)
