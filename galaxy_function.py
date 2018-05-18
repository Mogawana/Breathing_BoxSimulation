##################################################################
#
# Author: Orapeleng Mogawana
#   This file was written by Orapeleng Mogawana (o.mogawana@gmail.com) for Breathing_Box model##

##################################################################
#   The main routine solves the first oder differential equations using euler method
#   To generate numerical solutions for the gas conversion rate into stars and versa-versa


# Assuming the instantaneous mixing approximation
# (IMA), and that infall occurs after outflows in each timestep

def enrich_gas(sfr, R, t):
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
        Total remaining gas mass """

    # Enrich the gas from SNe
    Mg_e = (-sfr + R*sfr)*t

    return Mg_e


def remove_egas(sfr, alp, t):
    """This function calculates the outflow gas mass driven my supernovae
    and stellar winds

    Parameters
    ----------
    sfr:float
        star-formation rate.
    alp:float
        Outflow efficiency
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining gas mass """

    # Remove this newly-enriched, fully-mixed gas from the box
    Mg_r = (- alp*sfr)*t

    return Mg_r


def gas_infall(sfr, bta, t):
    """This function calculates the infall gas mass from outside the box.

    Parameters
    ----------
    sfr:float
        star-formation rate.
    bta:float
        Infall rate
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining gas mass """

    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mg_i = (bta*sfr)*t

    return Mg_i


def total_stellar_mass(sfr, R, t):
    """This function executes star formation and gas return from 
    supernovae.

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


def enrich_metalgas(zg, sfr, yz, t):
    """This function calculates the new metal enrichment of gas by 
    supernovae
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
        Total remaining metallicity in gas phase 
    """
    # enrich the gas with  new  metals
    Mgz_e = (-zg*sfr + yz*sfr)*t

    return Mgz_e


def remove_metalgas(zg, sfr, alp, t):
    """This function calculates fraction of the new metals removed from the 
    gas by supernovae
    Parameters
    ----------
    zg :float
        Gas metallicity
    sfr:float
        star-formation rate.
    alp:float
        Outflow efficiency
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining metallicity in gas phase 
    """
    # Remove this newly-enriched, fully-mixed gas from the box
    Mgz_r = (- alp*zg*sfr)*t

    return Mgz_r


def infall_metalgas(sfr, bta, znf, t):
    """This function calculates fraction of the new metals removed from the gas by supernovae
    Parameters
    ----------
    sfr:float
        star-formation rate.
    bta:float
        Infall rate
    znf:float
        infall gas metallicity
    t  :float
        Integration time step
    Return
    ---------
    float
        Total remaining metallicity in gas phase 
    """
    # Remove this newly-enriched, fully-mixed gas from the box
    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mgz_i =(bta*znf*sfr)*t
    
    return Mgz_i


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

def star_from_rate(alp_e, mg, mc, tdy):
    """"This function calculates the star fromation rate 
    in each time step.

    Parameters
    ----------
    apl_e :float
            star formation effeciency
    mg    :float
            total gas mass.
    mc    :float
            critical mass threshold
    tdy   :float
            disc dynamical time

    Return
    ---------
    float:
        Star formation rate at each time step """

    sfr = ( alp_e*(mg - mc) ) / tdy

    return sfr


def infall_rate(R, sfr):
    """"This function calculates the star fromation rate 
    in each time step.

    Parameters
    ----------
    R   :float
         fractional gas mass returned by SNe
    sfr :float
         star formation rate

    Return
    ---------
    float:
        infall rate of accreting gas """
    bta = (1 - R)*sfr

    return bta

class simulation:
    def __init__(self, t, Mgas, 
                 Zgas, Mstar, 
                 Zstar, Mzgas, 
                 Mzstar, MstarSP, 
                 ZstarSP):

        """Initial conditions place holder class object to initialize and
        evolve with the simulation"""
        self.t        = t 
        self.Mgas     = Mgas     # total mass of gas
        self.Zgas     = Zgas     # total gas metalicity
        self.Mstar    = Mstar    # total stellar mass
        self.Zstar    = Zstar    # total stellar metallicity
        self.Mzgas    = Mzgas    # total metal mass in gas phase
        self.Mzstar   = Mzstar   # total metal mass in stellar phase
        self.MstarSP  = MstarSP  # total mass of stellar population
        self.ZstarSP  = ZstarSP  # total metal mass in stellar population
