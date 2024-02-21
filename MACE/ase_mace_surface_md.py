from ase import units
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.md.verlet import VelocityVerlet
import ase.io.castep
import ase.md.analysis
import getpass
import os
from mace.calculators import mace_anicc, mace_mp, mace_off
from ase.visualize import view
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, fftfreq
from ase.io.trajectory import Trajectory
# from scipy.ndimage.filters import gaussian_filter1d as gaussian
from scipy.ndimage import gaussian_filter1d as gaussian
from ase import build


# List only the top level folders in a directory
def folder_list(mypath=os.getcwd()):
    """
    This function lists only the top level folders in a directory specified by mypath.
    If no directory is specified, it defaults to the current directory.

    Parameters:
    mypath (str): The path of the directory to be listed. Defaults to the current directory.

    Returns:
    onlyfolders (list): A list of top level folders in the specified directory.
    """
    onlyfolders = [f for f in os.listdir(mypath) if os.path.isdir(os.path.join(mypath, f))]
    return onlyfolders


def file_list(mypath=os.getcwd()):
    """
    This function lists only the files in a directory specified by mypath.
    If no directory is specified, it defaults to the current directory.

    Parameters:
    mypath (str): The path of the directory to be listed. Defaults to the current directory.

    Returns:
    onlyfiles (list): A list of files in the specified directory.
    """
    onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    return onlyfiles


def get_all_files(directory):
    """
    This function walks through a given directory and its subdirectories,
    and returns a list of all files in them.

    Parameters:
    directory (str): The path of the directory to be walked through.

    Returns:
    file_list (list): A list of paths for all files in the given directory and its subdirectories.
    """
    file_list = []
    for root, directories, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_list.append(file_path)
    return file_list


# List only files which contain a substring
def sub_file_list(mypath, sub_str):
    """
    This function lists only the files in a specified directory that contain a given substring.

    Parameters:
    mypath (str): The path of the directory to be listed.
    sub_str (str): The substring to filter the files by.

    Returns:
    list: A list of files in the specified directory that contain the given substring.
    """
    return [i for i in get_all_files(mypath) if sub_str in i]


def calc_vel_com(trajectory):
    """
    This function calculates the centre of mass velocity for a given trajectory.

    Parameters:
    trajectory (list): A list of Atom objects, each representing a state in the trajectory.

    Returns:
    vel_com (numpy.ndarray): A 2D array with the x, y, and z components of the centre of mass velocity for each state in the trajectory.

    The function works by iterating over each state in the trajectory. For each state, it calculates the total mass and the mass-weighted velocity for each dimension (x, y, z). The mass-weighted velocity is then divided by the total mass to get the centre of mass velocity for that dimension and state. The results are stored in a 2D array, with one row for each dimension and one column for each state.
    """
    # Initialize a 2D array to store the centre of mass velocities.
    # The array has one row for each dimension (x, y, z) and one column for each state in the trajectory.
    vel_com = np.zeros([len(trajectory), 3], dtype=float)

    # Iterate over each state in the trajectory.
    for i, atoms in enumerate(trajectory):
        # Get the masses of all atoms in the current state.
        masses = atoms.get_masses()
        # Calculate the total mass of all atoms in the current state.
        total_mass = np.sum(masses)

        # Iterate over each dimension (x, y, z).
        for j in range(3):
            # Calculate the mass-weighted velocity for the current dimension and state.
            # Multiplying the velocities of all atoms by their masses,
            # Summing the results, and dividing by the total mass.
            vel_com[i, j] = np.sum(atoms.get_velocities()[:, j] * masses) / total_mass

    # Return the 2D array with the centre of mass velocities.
    return vel_com


def calc_pdos(V, dt):
    """
    Calculate the phonon density of states from a trajectory of
    velocities (power spectrum of the velocity auto-correlation
    function).

    Parameters
    ----------
    V : :obj:`numpy.ndarray`
        (dims N x T) velocities of N degrees of freedom for
        trajetory of length T
    dt : float
        time between steps in trajectory (fs)

    Returns
    -------
    freq : :obj:`numpy.ndarray`
        (dims T) Frequencies (cm^-1)
    pdos : :obj:`numpy.ndarray`
        (dims T) Density of states (a.u.)
    """

    V = V.reshape(V.shape[0], -1).T

    n_steps = V.shape[1]

    # mean velocity auto-correlation for all degrees of freedom
    vac2 = [np.correlate(v, v, 'full') for v in V]
    vac2 /= np.linalg.norm(vac2, axis=1)[:, None]
    vac2 = np.mean(vac2, axis=0)

    # power spectrum (phonon density of states)
    pdos = np.abs(fft(vac2)) ** 2
    pdos /= np.linalg.norm(pdos) / 2  # spectrum is symmetric

    freq = fftfreq(2 * n_steps - 1, dt) * 33356.4095198152  # Frequency in cm^-1

    return freq[:n_steps], pdos[:n_steps]


dt = 1.0
nt = 100

user = getpass.getuser()
load_dir = r"C:\Users\{}\OneDrive - University of Surrey\Papers\paper_water_surface_diffusion\data".format(user)
directories = folder_list(load_dir)
files = sub_file_list(os.path.join(load_dir, directories[0]), ".md")
file = files[0]
print(file)

# Load the MD file
f = open(file, 'r')
traj = ase.io.castep.read_castep_md(f)
f.close()

# build water dimers
atoms = build.molecule('H2O')
atoms.center(vacuum=2.0)
# atoms = atoms.repeat((2, 1, 1))
atoms.center(vacuum=8.0)
atoms.set_pbc([True, True, True])

# load the md file
atoms = traj[0]
atoms.set_pbc([True, True, True])

# Make the calculator
# calc = mace_mp(model="large", default_dtype="float64", device='cpu', dispersion=True)
calc = mace_mp(model="small", default_dtype="float32", device='cpu', dispersion=False)
# calc = mace_off(model="large", default_dtype="float64",device="cpu", dispersion=False)
atoms.calc = calc

# Set the momenta
MaxwellBoltzmannDistribution(atoms, temperature_K=150.0)
Stationary(atoms)  # zero linear momentum
ZeroRotation(atoms)  # zero angular momentum
# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, timestep=dt * units.fs, trajectory='md.traj', logfile='-', loginterval=1)
dyn.run(1000)

traj = read("md.traj", index=":")
view(traj)

# Calculate the centre of mass velocity
vel_com = calc_vel_com(traj)

# Read velocities
V = np.array([atoms.get_velocities() for atoms in traj])
# Subtract the center of mass velocity
V_sub = V - vel_com[:, None, :]
freq, pdos_sub = calc_pdos(V_sub, dt)

# pdos_sub = gaussian(pdos_sub, sigma=1)

# read velocities
V = np.array([atoms.get_velocities() for atoms in traj])
freq, pdos = calc_pdos(V, dt)
# smoothing of the spectrum removes numerical artifacts due to finite time trunction of the FFT
# pdos = gaussian(pdos, sigma=1)

# plot a histogram of the velocities
plt.hist(V_sub.flatten(), bins=1000, density=True, alpha=0.5)
plt.hist(V.flatten(), bins=1000, density=True, alpha=0.5)
plt.xlim(-1.0, 1.0)
plt.yscale("log")
plt.show()

plt.plot(freq, pdos_sub)
plt.plot(freq, pdos, color="red")
plt.xlim(0, 4000)
plt.yscale("log")
plt.xlabel('Frequency [cm$^{-1}$]')
plt.ylabel('Density of states [a.u.]')
plt.show()
