#!/usr/bin/env python
import rebound
import reboundx
import numpy as np
import sys 
import os

# Get args from terminal, specify P_o in units of P_i
def get_argvs():
    """Gets arguments passed through terminal, returns them as floats.
    M1[msun], q, P_i[days], m_planet[msun*1e-03], P_o/P_i, e_o"""
    name, m1, q, P_i, e_i, m3, P_o = sys.argv 
    return float(m1), float(q), float(P_i), float(e_i),\
        float(m3), float(P_o)

m1, q, P_i, e_i, m3, P_o = get_argvs()

# Construct filenames in order to access previously generated stellar evolution
# data
stars_filename = \
    "nbody_sims/{:.1f}_{:.1f}_{:d}_{:.1f}/{:.1f}_{:.1f}_{:d}_{:.1f}.npy"\
        .format(m1, q, int(P_i), e_i, m1, q, int(P_i), e_i)
stars_directory = "nbody_sims/{:.1f}_{:.1f}_{:d}_{:.1f}/"\
        .format(m1, q, int(P_i), e_i, m1, q, int(P_i), e_i)

# Load stellar evolution data into numpy array
stars = np.load(stars_filename)

# Normalize time to zero, change units to rebound units
t = stars[:, 0]*2*np.pi*1e06
t -= t[0]

M1 = stars[0, 1] # msun
M2 = stars[0, 2]
e = stars[:, 5]
P = 2*np.pi*(stars[:, 6]/365.25) # rebound units
a = stars[:, 7]  # au

def quick_simulation():
    """Run quick simulation without stellar evolution to check system stability,
    raise an error if the system is unstable. Determine the secular forced 
    eccentricity of the planet and use it as an initial condition in the
    main simulation."""
    sim_q = rebound.Simulation()
    m_p = m3*9.5458*1e-04 # to jupiter masses
    P_p = P_o*P[0]
    sim_q.add(m=M1)
    sim_q.add(m=M2, P=P[0], e=e[0]) 
    sim_q.add(m=m_p, P=P_p, e=0.) 
    sim_q.move_to_com() # Moves to the center of momentum frame 
    sim_q.integrator = "ias15"
    ps = sim_q.particles
    
    times = 2*np.pi*np.linspace(0, 500, 2000)
    e_temp = np.zeros(len(times))
    pomega = np.zeros(len(times))
   
    planet_ejected = False

    for i, time in enumerate(times):
        orbits = sim_q.calculate_orbits()
        if(orbits[1].a < 10):
            sim_q.integrate(time)
        else:
            planet_ejected = True 
        pomega[i] = orbits[1].pomega
        e_temp[i] = orbits[1].e
    
    if(planet_ejected == True):
        print("planet ejected")
        sys.exit() 

    # Calculate the centroid in ecos(pomega),esin(pomega) space, this is the
    # secular forced eccentricity

    x_cen = np.mean(e_temp*np.cos(pomega))
    y_cen = np.mean(e_temp*np.sin(pomega))
    e_forced = np.sqrt(x_cen**2 + y_cen**2)

    return e_forced
    
e_forced = quick_simulation()

# Initialize simultion
sim = rebound.Simulation() 

m_p = m3*9.5458*1e-04 # Planet mass in jupiter masses
P_p = P_o*P[0]

sim.add(m=M1, hash='star1') 
sim.add(m=M2, P=P[0], e=e[0], hash='star2') 
sim.add(m=m_p, P=P_p, e=e_forced, hash='planet') 

sim.move_to_com() # Moves to the center of momentum frame 
sim.integrator = "ias15"
ps = sim.particles 

# Initilize reboundx objects 
rebx = reboundx.Extras(sim) 
mod_effect = rebx.add("modify_orbits_forces") 
 
# Create array of t with more points between each binary_c step
def denser_array(t, n):
    """Takes an array of binary_c times, returns new array with n times 
    more linearly spaced points in between each two consecutive points in 
    the original array."""
    t_dense = np.linspace(t[0], t[1], n) 
    for j in range(1, len(t) - 1): 
        array = np.linspace(t[j], t[j + 1], n) 
        t_dense = np.concatenate([t_dense, array]) 
    return t_dense

n = 10000 
t_dense = denser_array(t, n)
Nout = len(t_dense) 

# Create directory for output binary files
output_directory = stars_directory + "{:.2f}_{:.2f}_dense"\
            .format(m3, P_o)

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

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

def adiabatic_criterion(m1, m2, m3, P_i, e_i, e_o, n_res, Pdot):
    """
    Test wheater the resonance crossing is adiabatic.
    
    Given the system parameters, this function calculates the ratio between
    the librational timescale and the resonance width crossing timescale defined
    by Pdot of the binary. The calculations are based on the resonance widths 
    derived in Mardling (2013).
    
    Parameters
    ----------
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
    e_o : float
        Eccentricity of the outer binary.
    n_res: int
        Resonance order.
    Pdot: float
        Time derivative of P_i in units of days/Myr
    """
    # Calculate resonance widths
    widths = resonance_width(m1, m2, m3, e_i, e_o, n_res)
    P_o = n_res*P_i
    nu_o = 2*np.pi/P_o # Outer mean motion
    omega_n = nu_o*widths/2 # Libration frequency
    lib_period = 2*np.pi/omega_n # Libration period
            
    crossing_time = 2*widths*P_i/Pdot
    
    timescale_ratio = np.abs(lib_period/crossing_time)
    
    return timescale_ratio

# Test wheather the adiabatic criterion is satisfied
Pdot = np.zeros(len(t))
Pdot[:-1] = np.diff(P)/np.diff(t) # days/(2pi*yr)

# Find time at which a 6:1 or 7:1 resonance is crossed
if(P_p/P[0] < 6):
    n_res = 6
    for i, p in enumerate(P):
        if(P_p/p > 6):
            index = i
            break
else:
    n_res = 7
    for i, p in enumerate(P):
        if(P_p/p > 7):
            index = i 
            break

# Calculate timescale ratio at resonance crossing
timescale_ratio = adiabatic_criterion(M1, M2, m_p, P, e, e_forced,
        n_res, Pdot)

with open(stars_directory  + "initial_log.txt", "w") as f:
    f.write("Total integration time: {0}\n".format(t[-1]/(2*np.pi*1e06)))
    f.write("Total number of files: {0}\n".format(Nout))
    f.write("Time at resonance crossing {0}\n".format(t[index]/(2*np.pi*1e06)))
    f.write("e_i at resonance crossing {0}\n".format(e[index]))
    f.write("Ratio between libration timescale and the resonance crossing\
    timescale for {0}:1 resonance: {1}\n".format(n_res, timescale_ratio[index]))

for i, time in enumerate(t_dense): 
    # Update e-folding timecales every n steps 
    if(i % n == 0): 
        k = int(i/n) 
        delta_t = t[k + 1] - time 
        
        # Avoid problems with log0 and division by zero
        large_number = 1e30
        log_a = np.log(a[k + 1]/ps[1].a) 

        if(e[k + 1] > 0):
            log_e = np.log(e[k + 1]/ps[1].e) 
        else:
            log_e = 0.

        # Avoid division by zero
        if (log_a == 0):
            tau_a = large_number
        else:
            tau_a = delta_t/log_a

        if (log_e == 0):
            tau_e = large_number
        else:
            tau_e = delta_t/log_e

        ps[1].params["tau_e"] = tau_e 
        ps[1].params["tau_a"] = tau_a 

        # Stop the integration if the planet is ejected
        if (ps[2].a > 10):
            with open("planet_ejected.txt", "w") as f:
                f.write("Planet ejected during main Simulation")
            sys.exit()

    sim.integrate(time)

    # Save data to binary file
    nbody_filename = output_directory + "/{:.5f}.bin"\
            .format(time/(2*np.pi*1e06))
    sim.save(nbody_filename)
