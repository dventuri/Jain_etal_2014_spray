from math import pi, sqrt
from CoolProp.CoolProp import PropsSI


def circ_area(Diameter):
    '''Calculates circular area given a diameter'''
    Area = pi*Diameter*Diameter/4
    return Area


def th_flowrate(Area_nozzle, dP_nozzle, rho_fluid):
    '''Calculates the theoretical flow rate based on
       Bernoulli equation'''
    theoretical_flow_rate = Area_nozzle*sqrt(2*dP_nozzle/rho_fluid)
    return theoretical_flow_rate


def act_flowrate(discharge_coefficient, theoretical_flowrate):
    '''Calculates the actual flow rate based on the
       theoretical flow rate and coefficient of discharge'''
    actual_flow_rate = discharge_coefficient*theoretical_flowrate
    return actual_flow_rate


def flowrate2velocity(flow_rate, area):
    '''Calculates velocity through an area given the volumetric
       flow rate'''
    velocity = flow_rate/area
    return velocity


def calc_Weber(rho_fluid, velocity, nozzle_diameter, surface_tension):
    '''Calculates the Weber number'''
    Weber_number = rho_fluid*velocity*nozzle_diameter/surface_tension
    return Weber_number


def calc_Reynolds_p(nozzle_diameter, rho_fluid, pressure, viscosity_fluid):
    '''Calculates the Reynolds number based on the
       atomization pressure'''
    pressure_Re_number = nozzle_diameter*sqrt(rho_fluid*pressure)/viscosity_fluid
    return pressure_Re_number


def calc_spray_cone_angle(nozzle_diameter, swirl_chamber_diameter, Weber_number,
                          pressure_based_Reynolds_number):
    '''Calculates the spray cone angle using Jain et. al. (2014)
       suggested correlation'''
    spray_angle = (
        334.32**(-0.165) *
        (swirl_chamber_diameter/nozzle_diameter)**(-0.484) *
        Weber_number**(0.043) *
        pressure_based_Reynolds_number**(-0.065)
    )/2.0
    return spray_angle


if __name__ == '__main__':

    # atomization pressure
    dP = 6          # [bar]
    dP = dP*1e5     # convert [bar] -> [Pa]
    Patm = 101325   # [Pa]
    Tamb = 298.15   # [K]

    # liquid density (water @ Tamb, Patm)
    rho_l = PropsSI('D', 'T', Tamb, 'P', Patm, 'water')   # [kg/m³]

    # liquid viscosity (water @ Tamb, Patm)
    mu_l = PropsSI('V', 'T', Tamb, 'P', Patm, 'water')    # [Pa s]

    # liquid-gas surface tension (water-air @ Tamb, dP+Patm)
    # sigma = PropsSI('I', 'T', Tamb, 'P', Patm, 'water')
    sigma = 71.97       # [dynes/cm]
    sigma = sigma/1000  # convert [dynes/cm] -> [N/m]

    # coefficient of discharge
    Cd = 0.7213         # adimensional

    # swirl chamber diameter
    Ds = 6.0            # [mm]
    Ds = Ds/1000        # convert [mm] -> [m]

    # nozzle outlet orifice diameter
    Do = 1.65           # [mm]
    Do = Do/1000        # convert [mm] -> [m]

    # nozzle outlet orifice area
    Ao = circ_area(Do)  # [mm²]

    # theoretical flow rate (Bernoulli)
    Qth = th_flowrate(Ao, dP, rho_l)        # [m³/s]

    # actual flow rate
    Qact = act_flowrate(Cd, Qth)            # [m³/s]

    # nozzle mean velocity
    Um = flowrate2velocity(Qact, Ao)

    # Weber number
    We = calc_Weber(rho_l, Um, Do, sigma)       # adimensional

    # pressure-based Reynolds number
    Rep = calc_Reynolds_p(Do, rho_l, dP, mu_l)  # adimensional

    # spray cone angle
    beta = calc_spray_cone_angle(Do, Ds, We, Rep)   # [rad]
    beta = beta*180.0/pi

    # specific heat capacity
    cp_l = PropsSI('C', 'P', Patm, 'T', Tamb, 'water')

    # latent heat
    H_l = PropsSI('H', 'P', Patm, 'Q', 0, 'water')
    H_v = PropsSI('H', 'P', Patm, 'Q', 1, 'water')
    latent_heat = H_v - H_l

    # saturation temperature
    Tsat = PropsSI('T', 'P', Patm, 'Q', 0, 'water')

    # for AIR:
    rho_air = PropsSI('D', 'T', Tamb, 'P', Patm, 'air')
    mu_air = PropsSI('V', 'T', Tamb, 'P', Patm, 'air')
