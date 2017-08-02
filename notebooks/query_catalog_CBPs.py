import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET, urllib.request, gzip, io

def find_CBPs(detection_method):  
    """This function queries the Open Exoplanet Catalogue for circumbinary planets
    and returns all detections as a Pandas DataFrame object."""  

    # Link to OEC
    url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
    # Load XML tree into memory
    oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.request.urlopen(url).read())))
    print(type(oec))

    # Find all circumbinary planets (i.e., children of <binary> tag)
    print("Total number of CBPs: ", len(oec.findall(".//binary/planet")), "\n")

    for planet in oec.findall(".//binary/planet"): # note, "//" searches through all sublevels
        print(planet.findtext("name"), "-", planet.findtext("discoverymethod"))  

    # Find all binary systems with multiple CBPs
    print("\nSystems with multiple CBPs: \n")
    for binary in oec.findall(".//binary"):
        n_p = len(binary.findall("./planet"))
        if n_p > 1:
            print(binary.findtext(".//name"),",", n_p,"planets")


    # Store all CBPs in a pandas DataFrame object
    columns = ['m', 'm_erru', 'm_errorl', 'R', 'a', 'P',
               'e', 'I', 'mstar1', 'mstar2', 'q', 'a_binary', 'P_binary', 'e_binary']
    rows_transit = []
    rows_timing = []

    # Find all transit CBPs, i.e., MS CBPs detected by Kepler
    for planet in oec.findall(".//binary/planet"): # note, "//" searches through all sublevels
        if (planet.findtext("discoverymethod")=="transit"):
            rows_transit.append(planet.findtext("name"))

    # Find all ETV CBPs, i.e., PCE CBPs
    for planet in oec.findall(".//binary/planet"): 
        if (planet.findtext("discoverymethod")=="timing"):
            rows_timing.append(planet.findtext("name"))   

    # Create DataFrame object
    planets_transit = pd.DataFrame(columns=columns, index=rows_transit)
    planets_timing = pd.DataFrame(columns=columns, index=rows_timing)

    # Iterate over each planet in binary systems, fill in the column data
    for binary in oec.findall(".//binary"):

        # get masses of stars in binary
        mstars = [0, 0]
        for i, star in enumerate(binary.findall("./star")):
            # check if mass string is either empty or None
            if star.findtext("mass") is not None and star.findtext("mass"):
                mstars[i] = float(star.findtext("mass"))
        mstars = np.sort(np.array(mstars))
        mstars = mstars[::-1]

        # iterate over planets in binary
        for planet in binary.findall("./planet"):
            mass_errors = [0, 0]
            if planet.find("./mass") is not None:
                mass_errors[0] = planet.find("./mass").get('errorplus')
                mass_errors[1] = planet.find("./mass").get('errorminus')

            data = [
                    planet.findtext("mass"),
                    mass_errors[0],
                    mass_errors[1],
                    planet.findtext("radius"),
                    planet.findtext("semimajoraxis"),
                    planet.findtext("period"),
                    planet.findtext("eccentricity"),
                    planet.findtext("inclination"),
                    mstars[0],
                    mstars[1],
                    mstars[0]/mstars[1],
                    binary.findtext("semimajoraxis"),
                    binary.findtext("period"),
                    binary.findtext("eccentricity")
                ]

            if (planet.findtext("discoverymethod")=="transit"):
                planets_transit.loc[planet.findtext("name")] = data

            if (planet.findtext("discoverymethod")=="timing"):
                planets_timing.loc[planet.findtext("name")] = data

    # Convert all values in DataFrame to floats
    planets_transit = planets_transit.apply(pd.to_numeric, errors='raise')
    planets_timing = planets_timing.apply(pd.to_numeric, errors='raise')
    
    if detection_method=='transit':
        return planets_transit
    if detection_method=='timing':
        return planets_timing
