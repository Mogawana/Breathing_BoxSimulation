##################################################################
#
# Author: Orapeleng Mogawana
#   This file was written by Orapeleng Mogawana (o.mogawana@gmail.com) for Breathing_Box ##

##################################################################
#   The main routine solves the first oder differential equations using euler method
#   To generate numerical solutions for the gas conversion rate into stars and versa-versa


# Assuming the instantaneous mixing approximation
# (IMA), and that infall occurs after outflows in each timestep

def total_gas_mass(sfr, R, alp, bta, t):
    """This function executes star formation and gas return from supernovae.
    It futher calculates the gas outflows and infalls based on IMA

    Parameters
    ----------
    sfr:float
        star-formation rate.
    R  :float
        Super Novae returned fraction
    alp:float
        Outflow efficiency
    bta:float
        Infall rate
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining gas mass """

    # Enrich the gas from SNe (Fractional mass returned = 0.43)
    Mg_e = (-sfr + R*sfr)*t

    # Remove this newly-enriched, fully-mixed gas from the box
    Mg_r = (- alp*sfr)*t

    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mg_i = (bta*sfr)*t
    
    return Mg_e, Mg_r, Mg_i

def total_stellar_mass(sfr, R, t):
    """This function executes star formation and gas return from supernovae.

    Parameters
    ----------
    sfr:float
        star-formation rate.
    R  :float
        Super Novae returned fraction
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining stellar mass """

    # Execute star formation and gas return from supernovae
    return (sfr - R*sfr)*t

def metal_mass_gasphase(zg, sfr, yz, alp, bta, znf, t):
    """This function executes new metal enrichment of gas by Super Novae
    and infall gas from outside the box. It futher accounts the gas outflows
    for IMA

    Parameters
    ----------
    zg :float
        Gas metallicity
    sfr:float
        star-formation rate.
    yz :float
        Metal yield
    alp:float
        Outflow efficiency
    bta:float
        Infall rate
    znf:float
        infall gas metallicity
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining metallicity in gas phase """

    # enrich the gas with  new  metals
    Mgz_e = (-zg*sfr + yz*sfr)*t

    # Remove this newly-enriched, fully-mixed gas from the box
    Mgz_r = (- alp*zg*sfr)*t

    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mgz_i =(bta*znf*sfr)*t
    
    return Mgz_e, Mgz_r, Mgz_i
    
def metal_mass_starphase(zg, sfr, yz, t):
    """"This function executes new metal enrichment of gas by Super Novae
    and infall gas from outside the box.

    Parameters
    ----------
    zg :float
        Gas metallicity
    sfr:float
        star-formation rate.
    yz :float
        Metal yield
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining metallicity in stellar phase """

    # enrich the stars with  new  metals
    return (zg*sfr - yz*sfr)*t

def total_mass_starpop(sfr, R, t):
    """This function calculates the total mass of stellar population 
    in each time step.

    Parameters
    ----------
    sfr:float
        star-formation rate.
    R  :float
        Super Novae returned fraction
    t  :float
        Integration time step

    Return
    ---------
    float:
        Total remaining stellar mass population"""

    # total mass of each stellar population
    return (1 - R)*sfr*t

def total_mass_metalpop(zg, sfr, R, t):
    """"This function calculates the total metal mass of stellar population 
    in each time step.

    Parameters
    ----------
    zg :float
        Gas metallicity
    sfr:float
        star-formation rate.
    R  :float
        Super Novae returned fraction
    t  :float
        Integration time step

    Return
    ---------
    float:
        Total remaining stellar metal mass population """

    # total metal mass of each stellar population
    return zg*(1 - R)*sfr*t