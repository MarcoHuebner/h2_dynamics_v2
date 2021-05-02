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
import scipy
import matplotlib
import numpy as np
import scipy.sparse as ss
import scipy.sparse.linalg as ssla

from typing import Tuple
from matplotlib import animation


# reduced molecule mass ( mu = m1*m2/(m1 + m2) = m_p/2 ) in atomic units (m_p in electron masses)
mu = 1836/2


def load_data(filename: str, filepath: str = "data") -> np.ndarray:
    """Load data from a tab separated csv file.

    Args:
        filename (str): Name of the csv file.
        filepath (str): The path of the file.

    Returns:
        (np.array) A matrix containing the data from the file.
    """
    return np.loadtxt(os.path.join(filepath, filename), delimiter="\t")


def preprocess_data(dip_coupling: np.ndarray, wave_fct: np.ndarray, pot_uneven: np.ndarray,
                    pot_even: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate loaded .dat files to hard-coded range of 300 with steps of 0.05 atomic units.

        Args:
            dip_coupling (np.array): Dipole coupling array.
            wave_fct (np.array): Given initial wave function.
            pot_uneven (np.array): Photon-excited, anti-binding potential.
            pot_even (np.array): Ground state potential.

        Returns:
            (Tuple[np.array*4]) The interpolated arrays.
    """
    # global variables
    unit_limit = 300.0  # should be a multiple of dx, defines total grid size
    dx = 0.05           # constant for fineness of grid

    # interpolate and add the first missing dipole coupling values
    interpolation = np.stack((
        np.linspace(dip_coupling[0, 0], dip_coupling[1, 0], 7, endpoint=True),
        np.linspace(dip_coupling[0, 1], dip_coupling[1, 1], 7, endpoint=True)
    ), axis=1)

    dipole_coupling = np.concatenate([interpolation[1:-1], dip_coupling[1:]])

    # interpolate between values at 10 and 1000 Bohr radii for larger grids
    if unit_limit > 10.0:
        units = np.linspace(10.0 + dx, unit_limit, int((unit_limit - 10.0) / dx), endpoint=True)
        interpolation = np.stack((units, 0.494112 * units), axis=1)
        dip_coupling = np.concatenate([dipole_coupling, interpolation])

    # interpolate potentials and wave function for larger grids
    if unit_limit > 30.0:
        units = np.linspace(30.0 + dx, unit_limit, int((unit_limit - 30.0) / dx), endpoint=True)

        wave_fct = np.concatenate([
            wave_fct[wave_fct[:, 0] <= 30.0],
            np.stack((units, np.zeros(len(units))), axis=1)
        ])

        pot_uneven = np.concatenate([
            pot_uneven[pot_uneven[:, 0] <= 30.0],
            np.stack((units, np.zeros(len(units))), axis=1)
        ])

        pot_even = np.concatenate([
            pot_even[pot_even[:, 0] <= 30.0],
            np.stack((units, -6.66667e-06 * np.ones(len(units))), axis=1)
        ])

    # limit to desired range, ensure same grid (discard values at 0.0 due to diverging potential)
    pot_uneven = pot_uneven[pot_uneven[:, 0] <= unit_limit]
    pot_even = pot_even[pot_even[:, 0] <= unit_limit]
    wave_function = wave_fct[wave_fct[:, 0] <= unit_limit][1:]
    dipole_coupling = dip_coupling[dip_coupling[:, 0] <= unit_limit]

    return dipole_coupling, wave_function, pot_uneven, pot_even


def hamiltonian(potential: np.ndarray, dx: float) -> scipy.sparse.csr_matrix:
    """Compute the hamiltonian for the given potential.

    Args:
        potential (np.array): The potential as vector containing the diagonal.
        dx (float): step size of hamiltonian.
    Returns:
        (scipy.sparse.csr_matrix): The hamiltonian as a sparse matrix.
    """
    # diagonals of kinetic energy matrix (second order diff. quot.) + potential
    main_diag = (1.0 / (dx**2 * mu)) * np.ones(len(potential)) + potential
    off_diag = (-0.5 / (dx**2 * mu)) * np.ones(len(potential)-1)

    return ss.diags([off_diag, main_diag, off_diag], [-1, 0, 1])


def save_anim(anim: matplotlib.animation.FuncAnimation, filename: str, fps: int = 1) -> None:
    """Load data from a tab separated csv file.

        Args:
            anim (matplotlib.animation.FuncAnimation): Animation figure.
            filename (str): The name of the file, should be a .gif.
            fps (int): Frames to show per second in .gif.
    """
    writer = animation.PillowWriter(fps=fps)
    anim.save('plots/' + filename, writer=writer, dpi=300)

    return None


class TimeEvolution:
    """Functor for the time evolution of a given initial state and time independent hamiltonian.

    Args:
        initial_state (np.array): The initial state.
        hamiltonian (np.array): The Hamiltonian of the system.
    """
    def __init__(self, initial_state: np.ndarray, hamiltonian: np.ndarray) -> None:
        assert initial_state.shape[0] == hamiltonian.shape[0], "Make sure initial state dim is compatible with " \
                                                               "hamiltonian dim for matrix multiplication"

        # compute eigenvals (-energies) and eigenvecs (-states)
        eigenvals, eigenvecs = ssla.eigsh(hamiltonian, k=hamiltonian.shape[0])

        # select only bound states
        mask = eigenvals <= 0.0
        
        self.initial_state = np.conjugate(eigenvecs[:, mask].T) @ initial_state
        self.hamiltonian = hamiltonian
        self.eigenvals = eigenvals[mask]
        self.eigenvecs = eigenvecs[:, mask]

    def __call__(self, time: float) -> np.ndarray:
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
    def __init__(self, num_points: int, step_size: float, potential: np.ndarray = None) -> None:
        self.kvalues_squared = (2 * np.pi * np.fft.fftfreq(num_points, d=step_size))**2
        
        self.num_points = num_points
        self.step_size = step_size
        self.potential = potential

    def __call__(self, state: np.ndarray, dt: float, potential: np.ndarray = None) -> np.ndarray:
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
        tau (float): time delay of one oscillation period
    """
    def __init__(self, pot_even: np.ndarray, pot_uneven: np.ndarray, dipole_coupling: np.ndarray,
                 tau: float = 1075.3) -> None:
        # provide constants taken from reference papers (in atomic units)
        self.tau = tau
        self.t_fwhm = 355.7
        self.omega = 0.0599
        self.f_0 = 0.01205

        self.dipole_coupling = dipole_coupling[:, 1]
        # define intermediate results for shorter code
        self.pot_sum = pot_even[:, 1] + (pot_uneven[:, 1] - self.omega)
        self.pot_diff = pot_even[:, 1] - (pot_uneven[:, 1] - self.omega)

    def __call__(self, t: float) -> np.ndarray:
        """Compute V(t).

        Args:
            t (float): a point in time

        Returns:
            (np.array) The potential as a vector
        """
        # compute eq. (11) of linked project report in smaller steps
        g_t = np.exp(-(2*np.sqrt(np.log(2)) * (t - self.tau) / self.t_fwhm)**2) 
        f_t = self.f_0 * g_t * np.cos(self.omega * (t - self.tau))
        # assemble above's steps to W(t)
        pot_laser = -f_t * self.dipole_coupling

        # calculate final steps for solved eigenvalue eq. (7) of coupled energy surfaces
        potential = np.sqrt(self.pot_diff**2 + 4 * np.real(pot_laser**2))

        return 0.5 * (self.pot_sum - potential)
