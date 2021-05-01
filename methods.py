"""
A collection of methods employed in the H2+ dynamics project.

The following methods are provided by this module:

* load_data: A simple helper function to load the different data files used in this project.
* hamiltonian: Computes the Hamiltonian for a given potential as a matrix.
* TimeEvolution: The time evolution of a wave function under a time independent Hamiltonian.
* SplitStepMethod: The time evolution of a wave function with an arbitrary Hamiltonian.
* LaserFieldPotential: The time dependent potential in the laser field.

In addition to the different methods a constant for the reduced molecular mass `mu` is defined.
"""

import os
import numpy as np
import scipy.sparse as ss


# reduced molecule mass ( mu = m1*m2/(m1 + m2) = m_p/2 ) in atomic units (m_p in electron masses)
mu = 1836/2


def load_data(filename, filepath="data"):
    """Load data from a tab separated csv file.

    Args:
        filename (str): Name of the csv file.
        filepath (str): The path of the file.

    Returns:
        (np.array) A matrix containing the data from the file.
    """
    return np.loadtxt(os.path.join(filepath, filename), delimiter="\t")


def hamiltonian(potential, dx):
    """Compute the hamiltonian for the given potential.

    Args:
        potential (np.array): The potential as vector containing the diagonal.
        dx (float): step size of hamiltonian.
    Returns:
        (scipy.spare.csr_matrix): The hamiltonian as a sparse matrix.
    """
    size = len(potential)

    # diagonals of kinetic energy matrix
    main_diag = (1.0 / (dx**2 * mu)) * np.ones(size) + potential
    off_diag = (-0.5 / (dx**2 * mu)) * np.ones(size-1)

    return ss.diags([off_diag, main_diag, off_diag], [-1, 0, 1])


class TimeEvolution:
    """Functor for the time evolution of a given initial state and time independent hamiltonian.

    Args:
        initial_state (np.array): The initial state.
        hamiltonian (np.array): The Hamiltonian of the system.
    """
    def __init__(self, initial_state, hamiltonian):
        # TODO: add some assertions to check dims
        eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)

        mask = eigenvals <= 0.0
        
        self.initial_state = np.conjugate(eigenvecs[:, mask].T) @ initial_state
        self.hamiltonian = hamiltonian
        self.eigenvals = eigenvals[mask]
        self.eigenvecs = eigenvecs[:, mask]

    def __call__(self, time):
        """Compute the state at the given time.

        Args:
            time (float): A point in time.

        Returns:
            (np.array) The state at the given time.
        """
        evolution = self.eigenvecs * np.exp(-1j * self.eigenvals * time)
        return np.sum(self.initial_state * evolution, axis=1)


class SplitStepMethod:
    """Implementation of the split-step fourier method.

    Note:
        If the potential is constant provide it on construction.
        Otherwise pass is to the call method.

    Args:
        num_points (int): The number of elements in frequency domain.
        step_size (float): The step size in the frequency domain.
        potential (np.array): The potential energy of the system (optional).
    """
    def __init__(self, num_points, step_size, potential=None):
        self.kvalues_squared = (2 * np.pi * np.fft.fftfreq(num_points, d=step_size))**2
        
        self.num_points = num_points
        self.step_size = step_size
        self.potential = potential

    def __call__(self, state, dt, potential=None):
        """Compute a single step with the split-step method.

        Args:
            state (np.array): The state of the system.
            dt (float): The size of the time step.
            potential (np.array): The potential energy if not already provided in the init method.

        Returns:
            (np.array): The state after the time step has been aplied.
        """
        potential = self.potential if potential is None else potential

        # apply potential
        state = np.exp(-1j * dt * potential) * state

        state = np.fft.fft(state)

        # apply only half (!) of the kinetic part
        state = np.exp(((-1j * dt) / (2 * mu)) * self.kvalues_squared) * state
  
        return np.fft.ifft(state)


class LaserFieldPotential:
    """Represents the time dependent potential inside the laser field.

    Args:
        pot_even (np.array): the binding part of the constant potential
        pot_uneven (np.array): the repellent part of the constant potential
        dipole_coupling (np.array): dipole coupling between states
        dx (float): step size of potential
        tau (float): time delay of one oscillation period
    """
    def __init__(self, pot_even, pot_uneven, dipole_coupling, dx, tau=1075.3):
        # TODO: add more documentation and maybe dx
        # provide constants taken from reference papers (in atomic units)
        self.tau = tau
        self.t_fwhm = 355.7
        self.omega = 0.0599
        self.f_0 = 0.01205

        self.dipole_coupling = dipole_coupling[:, 1]
        self.pot_sum = pot_even[:, 1] + (pot_uneven[:, 1] - self.omega)
        self.pot_diff = pot_even[:, 1] - (pot_uneven[:, 1] - self.omega)

    def __call__(self, t):
        """Compute V(t).

        Args:
            t (float): a point in time

        Returns:
            (np.array) The potential as a vector
        """
        # TODO: add more documentation
        g_t = np.exp(-(2*np.sqrt(np.log(2)) * (t - self.tau) / self.t_fwhm)**2) 
        f_t = self.f_0 * g_t * np.cos(self.omega * (t - self.tau))
        pot_laser = -f_t * self.dipole_coupling

        potential = np.sqrt(self.pot_diff**2 + 4 * np.real(pot_laser**2))
        return 0.5 * self.pot_sum - 0.5 * potential
