import rebound                                                                                                                                                           
import numpy as np                  
import matplotlib as mpl
from matplotlib import pyplot as plt    
import glob, os, re

# mpl.rc('font',**{'family':'serif','serif':['Palatino']})
# mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 20

def natural_sort(l): 
    """Sorts list of filenames numerically."""
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def plot_nbody_output(main_directory, directory, filename):
    num_of_files = 0
    for file in glob.iglob(directory + '/' + '*.bin'):
        num_of_files += 1

    a_i = np.zeros(num_of_files)
    e_i = np.zeros(num_of_files)
    P_i = np.zeros(num_of_files)
    a_o = np.zeros(num_of_files)
    e_o = np.zeros(num_of_files)
    P_o = np.zeros(num_of_files)
    omega_i = np.zeros(num_of_files)
    omega_o = np.zeros(num_of_files)
    lambda_i = np.zeros(num_of_files)
    lambda_o = np.zeros(num_of_files)
    t = np.zeros(num_of_files)

    for i, file in enumerate(natural_sort(glob.glob(directory + '/' + '*.bin'))):
        sim = rebound.Simulation.from_file(file)
        orbits = sim.calculate_orbits()
        a_i[i] = orbits[0].a
        e_i[i] = orbits[0].e
        P_i[i] = orbits[0].P/(2*np.pi)
        a_o[i] = orbits[1].a
        e_o[i] = orbits[1].e
        P_o[i] = orbits[1].P/(2*np.pi)
        omega_i[i] = orbits[0].omega*180/np.pi
        omega_o[i] = orbits[1].omega*180/np.pi
        lambda_i[i] = orbits[0].l*180/np.pi
        lambda_o[i] = orbits[1].l*180/np.pi
        t[i] = sim.t/(2*np.pi)
        del sim     
        
    # Make plot
    fig, ax = plt.subplots(3, 2, figsize=(25,25))

    ax[0, 0].plot(t/1e06, a_i, 'C0.')
    ax[0, 1].plot(t/1e06, e_i, 'C1.')
    ax[0, 0].set_ylabel(r'$a_i [au]$')
    ax[0, 1].set_ylabel(r'$e_i$')

    ax[1, 0].plot(t/1e06, a_o, 'C0.')
    ax[1, 1].plot(t/1e06, e_o, 'C1.')
    ax[1, 0].set_ylabel(r'$a_o [au]$')
    ax[1, 1].set_ylabel(r'$e_o$')
    
    ax[2, 0].plot(e_o*np.cos(omega_o - omega_i),e_o*np.sin(omega_o - omega_i) , 'C0.')
    ax[2, 1].plot(t/1e06, omega_o - omega_i, 'C0.')
    ax[2, 0].set_ylabel(r'$P_o/P_i$')
    ax[2, 1].set_ylabel(r'$\omega_o - \omega_i$')

    for a in ax.ravel():
        a.set_xlabel('t [Myr]')
        a.grid(True)

    ax[1,0].set_ylim(0, 3.5)
    ax[0,1].set_ylim(0, 0.35)
    ax[1,1].set_ylim(0, 0.35)

    plt.savefig(main_directory + '/' + filename)

# Iterate over all folders with nbody output, make plots and save them
M1 = np.array([1.2, 1.6, 2.0])
q = [0.3, 0.6, 0.8]
P = [20, 200]
e = [0., 0.1, 0.3]

from itertools import product
for M1_, q_, P_, e_ in product(M1, q, P, e):
    main_directory = 'nbody_sims/' + str(M1_) + '_' + str(q_) + '_'\
        + str(P_) + '_' + str(e_) 
    
    if(e_ > 0.1):
        directory = main_directory + '/' + '1.00_6.90'
        directory_dense = main_directory + '/' + '1.00_6.80_dense'
        n_res = 7
    else:
        directory = main_directory + '/' + '1.00_5.90'
        directory_dense = main_directory + '/' + '1.00_5.80_dense'
        n_res = 6
        
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory_dense):
        os.makedirs(directory_dense)
    
    plot_nbody_output(main_directory, directory, 'plot.png')
    plot_nbody_output(main_directory, directory_dense, 'plot_dense.png')


