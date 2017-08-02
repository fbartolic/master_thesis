import numpy as np
import pandas as pd
import ctypes

def evolve_binary(tmax, m1, m2, Z, P, e):
    """
    Evolves a binary system with the binary_c code.
    
    This function calls C code which uses the binary_c
    API to evolve a binary system. The output structure
    has to be defined in the C code.
    
    Parameters
    ----------
    tmax : float
        Maximum evolution time.
    M1 : float
        Mass of star 1.
    M2 : float
        Mass of star 2.
    Z : float
        Metallicity of the binary (assumed equal for both stars)
    P : float
        Inital orbital period of the binary.
    e : float
        Initial eccentricity.
        
    Returns
    -------
    DataFrame  
        Pandas DataFrame with the following output:
        t[Myr], m1[msun], m2[msun], R1[rsun], R2[rsun], e,
        P[days], a[au], type1, type2, in_RLOF, coel, merged
        
        type1 and type2 refer to stellar types. Possible stellar
        types are "Low mass MS", "MS", "Hertzsprung gap", "Giant branch",
        "Core Helium burning", "Early AGB", "TP AGB", "Naked Helium MS",
        "Naked Helium HG", "Naked Helium GB", "Helium WD", "CO WD", "ON WD",
        "NS","BH","Massless remnant". With Low mass MS designated by 0,
        MS by 1 and so on.
        
    """
    if not (m1 >= m2):
        raise AssertionError("binary_c requires that m1 > m2")

    # Load shared C library into Python via ctypes
    lib = np.ctypeslib.load_library('libevolve', '.')

    # Specify input and output types for called functions
    lib.evolve_binary.argtypes = [ctypes.c_char_p, ctypes.c_char_p] # input type
    lib.evolve_binary.restype = None # Output type 

    tmax=str(tmax); m1=str(m1); m2=str(m2); Z=str(Z); P=str(P); e=str(e)
    
    # Enter binary parameters
    params = bytes("binary_c M_1 " + m1 + " M_2 " + m2 + " metallicity " + Z\
              + " orbital_period " + P + " max_evolution_time " + tmax\
                   + " eccentricity " + e, 'utf8')
    
    # Call C function for system evolution
    buffer = ctypes.create_string_buffer(4000000)
    lib.evolve_binary(params, buffer)

    # Read byte values of the output into memory
    values = ctypes.cast(buffer, ctypes.c_char_p).value
    string = str(values, 'utf-8')

    # Remove last part of string about element yields
    string = string[:string.rfind('XYIELD')]

    array = np.fromstring(string, sep=",")
    n_columns = 13
    array = array.reshape(int(len(array)/n_columns), n_columns)
    array[:, 6] *= 365.25 # convert P units to days
    array[:, 7] *= 0.00465 # convert a units to au
    
    # Create DataFrame object
    columns = ['t', 'M1', 'M2', 'R1', 'R2', 'e',
        'P', 'a', 'type1', 'type2', 'in_RLOF', 'coel','merged']
    
    return pd.DataFrame(array, columns=columns, index=array[:, 0])
