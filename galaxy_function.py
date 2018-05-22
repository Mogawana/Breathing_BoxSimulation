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
    And calculates the gas mass remaining

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
        Total remaining gas mass 

    # Enrich the gas from SNe
    Mgas = (-sfr + R*sfr)*t
    """
    return (-sfr + R*sfr)*t


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
        Total remaining gas mass

    # Remove this newly-enriched, fully-mixed gas from the box
    Mgas = (- alp*sfr)*t
    """
    return (- alp*sfr)*t


def gas_infall(bta, t):
    """This function calculates the total gas mass accreted 
    from outside the box.

    Parameters
    ----------
    
    bta:float
        Infall rate
    t  :float
        Integration time step

    Return
    ---------
    float
        Total remaining gas mass

    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mgas = (bta)*t
    """
    return (bta)*t


def total_stellar_mass(sfr, R, t):
    """This function calculates the mass of stars formed by exercuting
    star formation and gas return.

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
        Total remaining stellar mass

    # Execute star formation and gas return from supernovae
    Mstar = (sfr - R*sfr)*t
    """
    return (sfr - R*sfr)*t


def enrich_metalgas(zg, sfr, yz, t):
    """This function calculates the new metal enrichment of gas by 
    supernovae.

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
    
    # enrich the gas with  new  metals
    Mzgas = (-zg*sfr + yz*sfr)*t
    """
    return (-zg*sfr + yz*sfr)*t


def remove_metalgas(zg, sfr, alp, t):
    """This function calculates fraction of the new metals removed from the 
    gas by supernovae.

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
    
    # Remove this newly-enriched, fully-mixed gas from the box
    Mzgas = (- alp*zg*sfr)*t
    """
    return (- alp*zg*sfr)*t


def infall_metalgas(bta, znf, t):
    """This function calculates fraction of the new metals accreted
    into the gas from outside the box.

    Parameters
    ----------
    
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
    
    # Remove this newly-enriched, fully-mixed gas from the box
    # Infall is assumed to come from outside the box,
    # and is added to the gas component
    # infall occurs after outflows in each timestep
    Mzgas =(bta*znf)*t
    """
    return (bta*znf)*t


def metal_mass_starphase(zg, sfr, R, t):
    """"This function calculates the metal mass content trapped in stars.

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
    float
        Total remaining metallicity in stellar phase 

    # enriched stars with  new  metals
    Mzstar = (zg*sfr - zg*R*sfr)*t
    """
    return (zg*sfr - zg*R*sfr)*t


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
        Total remaining stellar mass population

    # total mass of each stellar population
    MstarSP = (1 - R)*sfr*t
    """
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
        Total remaining stellar metal mass population

    # total metal mass of each stellar population
    MzstarSP = zg*(1 - R)*sfr*t
    """
    return zg*(1 - R)*sfr*t

def star_form_rate(alp_e, mg, mc, tdy):
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
        Star formation rate at each time step

    sfr = ( alp_e*(mg - mc) ) / tdy
    """
    return ( alp_e*(mg - mc) ) / tdy


def infall_rate(R, sfr):
    """This function calculates the gas infall rate in each
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
        infall rate of accreting gas
    beta = bta = (1 - R)*sfr
    """
    return (1 - R)*sfr

class simulation(object):
    """Initial conditions place holder class object to initialize and
    evolve with the simulation

    Parameters
    ----------
    t: float
        integration time
    Mgas: float
        total mass of gas
    Zgas: float
        total gas metalicity
    Mstar: float 
        total stellar mass
    Zstar: float 
        total stellar metallicity
    Mzgas: float
        total metal mass in gas phase
    Mzstar: float
        total metal mass in stellar phase
    MstarSP: float 
        total mass of stellar population
    ZstarSP: float
        metallicity in stellar population
    MzstarSP: float 
        total metal mass in stellar population
    SFR : float 
        star formation rate
    beta: float 
        infall rate

    Return
    ---------
    None
    """

    def __init__(self, t, Mgas, 
                 Zgas, Mstar, 
                 Zstar, Mzgas, 
                 Mzstar, MstarSP, 
                 ZstarSP, MzstarSP, 
                 SFR, beta):

        self.t        = t 
        self.Mgas     = Mgas
        self.Zgas     = Zgas
        self.Mstar    = Mstar
        self.Zstar    = Zstar
        self.Mzgas    = Mzgas
        self.Mzstar   = Mzstar
        self.MstarSP  = MstarSP
        self.ZstarSP  = ZstarSP
        self.MzstarSP = MzstarSP
        self.SFR      = SFR
        self.beta     = beta