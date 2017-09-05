import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
import rebound
from rebound.interruptible_pool import InterruptiblePool
import xml.etree.ElementTree as ET, urllib.request, gzip, io
from query_catalog_CBPs import find_CBPs

# mpl.rc('font',**{'family':'serif','serif':['Palatino']})
# mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 20

def simulation(params):
    """
    Runs a single rebound simulation and calculates the MEGNO value.
    
    This function sets up a single rebound simulation, and calculates the MEGNO value
    based on the approach described in the Rein & Tamayo paper on variational equations.
    
    Parameters
    -----------
    m1 : float
        Mass of the first body in msun.
    m2 : float
        Mass of the second body in msun.
    m3 : float
        Mass of the third body in msun.
    P_i : float
        Period of the 'inner binary' in days.
    e_i : float
        Eccentricity of the inner binary.
    P_o : float
        Period of the 'outer binary' in days.
    e_o : float
        Eccentricity of the outer binary
    
    Returns
    -------
    megno : float
        The MEGNO value for the chosen system parameters.
    """

    m1, m2, m3, P_i, e_i, P_o, e_o, integration_time = params
    
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.ri_whfast.safe_mode = 0
    
    # Convert periods to rebound units
    P_i *= 2*np.pi/365.25 # 2pi*years
    P_o *= 2*np.pi/365.25 
    
    sim.dt = P_i/10
    sim.add(m=m1)
    sim.add(m=m2, P=P_i, e=e_i) 
    sim.add(m=m3, P=P_o, e=e_o) 

    sim.move_to_com()
    sim.init_megno()
    sim.exit_max_distance = 20.
    try:
        sim.integrate(integration_time*2.*np.pi, exact_finish_time=0) # integrate 
        # for 500 years, integrating to the nearest
        #timestep for each output to keep the timestep constant and preserve WHFast's symplectic nature
        megno = sim.calculate_megno()
        return megno
    except rebound.Escape:
        return 10. # At least one particle got ejected, returning large MEGNO.

def resonance_width(m1, m2, m3, e_i, e_o, n):
    """
    Calculates the resonance width of an n-th order resonance in units of the 
    'inner binary' period. 
    """
    m123 = m1 + m2 + m3
    m12 = m1 + m2
    
    xi = np.arccosh(1/e_o) - np.sqrt(1 - e_o**2)
    width = (6*0.71**.5/((2*np.pi)**(1/4.)))*\
            (m3/m123 + n**(2/3.) *(m12/m123)**(2/3.)*(m1*m2/m12**2))**.5 *\
            (e_i**.5/e_o)*(1 - 13/24.*e_i**2.)**.5*(1 - e_o**2)**(3/8.)*\
            n**(3/4.)*np.exp(-n*xi/2.);
    return width

def stability_plot(params):
    """
    Plots a MEGNO map of a given system, together with resonant widths.
    
    For a given set of system parameters, the function plots a heat map 
    of MEGNO values as a function of period and 'outer binary' eccentricity.
    On top of this map it also plots the analytical prediction for resonant widths.
    
    Parameters
    ----------
    ax : matplotlib axes object
    
    m1 : float
        Mass of the first body in msun.
    m2 : float
        Mass of the second body in msun.
    m3 : float
        Mass of the third body in msun.
    P_i : float
        Period of the 'inner binary' in days.
    e_i : float
        Eccentricity of the inner binary."""

    m1, m2, m3, P_i, e_i, integration_time = params

    # Grid resolution
    N = 150 
    
    P_outer = np.linspace(1., 15.2, N)*P_i
    e_outer = np.linspace(0., 0.82, N)

    parameters = []

    for e in e_outer:
        for P in P_outer:
            parameters.append((m1, m2, m3, P_i, e_i, P, e, integration_time))

    # This part uses parallelization 
    pool = InterruptiblePool()
    results = pool.map(simulation, parameters)
    results2d = np.array(results).reshape(N, N)
    x = (P_outer - (P_outer[1] - P_outer[0])/2)/P_i # shift ticks to center of cell
    y = e_outer - (e_outer[1] - e_outer[0])/2
    
    return x, y, results2d

cbp_transit = find_CBPs('transit')

fig, ax = plt.subplots(3,2, figsize=(15,15), sharey=True, sharex=True)
fig.subplots_adjust(wspace=0.05, hspace=0.15)

k = 0

planet_list = ['Kepler-1647 b', 'Kepler-34 (AB) b', 'Kepler-35 (AB) b','Kepler-38 (AB) b',
              'Kepler-413 b', 'KIC 9632895 b','KOI-2939 b','PH-1 A(ab) b']

for i, name in enumerate(planet_list):        
    m1 = cbp_transit.loc[name,'mstar1']
    m2 = cbp_transit.loc[name,'mstar2']
    m3 = cbp_transit.loc[name,'m']*0.0009543
    a_i = cbp_transit.loc[name,'a_binary']
    P_i = cbp_transit.loc[name,'P_binary']
    e_i = cbp_transit.loc[name,'e_binary']
    P_o = cbp_transit.loc[name,'P']
    e_o = cbp_transit.loc[name,'e']
    if (np.isnan(P_i) == True):
        P_i = np.sqrt(a_i**3/(m1 + m2))*365.25
    if (np.isnan(m3) == True):
        m3 = 0.   

    if (P_o/P_i < 20):
        x, y, z = stability_plot([m1,m2,m3,P_i,e_i, 5000])
        im = ax.ravel()[k].pcolormesh(x, y, z, vmin=1.9, vmax=4)
        im.set_edgecolor('face') # remove white lines in some pdf viewers
        ax.ravel()[k].scatter(P_o/P_i, e_o, color='white')
        ax.ravel()[k].set_title(name + ', $e_i={:.2f}$'.format(e_i))

        #Plot orbit crossing line
        ax.ravel()[k].plot(x, 1 - ((m1 + m2)/(m1 + m2 + m3))**(1/3)*x**(-2/3), 'k--')

        # Plot resonance widths
        eo = np.linspace(0, 1, 500) 

        for i, n in enumerate(range(2, 16)):
            width = resonance_width(m1, m2, m3, e_i, eo, n)   
            ax.ravel()[k].fill_betweenx(eo, -width + n, width + n,
                                        facecolor='black', alpha=0.15)
            ax.ravel()[k].set_xlim(1, 14.9)

        k += 1
        
cb = fig.colorbar(im, ax=ax.ravel().tolist(), pad=0.02)
cb.set_label(r"MEGNO $\langle Y \rangle$")   
        
for a in ax.ravel():
    a.set_ylim(0, 0.8)
    a.set_xticks(np.arange(1, 16))
    a.set_yticks(np.linspace(0., 0.8, 9))
    
for a in ax[:, 0]:
    a.set_ylabel(r" $e_o$")
    
for a in ax[-1, :]:
    a.set_xlabel(r"$P_o/P_i$")
        
# cbp_transit.index.values
plt.savefig('megno_kepler.pdf', bbox_inches='tight')
