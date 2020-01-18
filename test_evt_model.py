#!/usr/bin/env python3
'''
Created on 25 Aug 2019

@author: jmht
'''

from evt_model import calc_saturated_vapour_pressure_air
from evt_model import calc_saturated_vapour_pressure_air_FAO
from evt_model import calc_stomatal_resistance
from evt_model import calc_net_radiation
from evt_model import calc_saturated_vapour_concentration_air
from evt_model import calc_vapour_concentration_air
from evt_model import calc_vapour_concentration_deficit
from evt_model import calc_vapour_concentration_surface
from evt_model import calc_sensible_heat_exchange
from evt_model import calc_latent_heat_flux
from evt_model import calc_epsilon


def test_saturated_vapour_pressure_air():
    temp_air = 25.0
    saturated_vapour_pressure = 3158
    assert abs(calc_saturated_vapour_pressure_air(temp_air) - saturated_vapour_pressure) < 1


def test_saturated_vapour_pressure_air_FAO():
    temp_air = 25.0
    saturated_vapour_pressure = 3.1688
    assert abs(calc_saturated_vapour_pressure_air_FAO(temp_air) - saturated_vapour_pressure) < 0.001


def test_saturated_vapour_concentration_air():
    """Values taken from: https://www.engineeringtoolbox.com/relative-humidity-air-d_687.html"""
    temp_air = 20.0
    ref_concentration = 17.2 # they use 17.8 - difference due to their use of 273 rather then 273.15 for zero Kelvin
    assert abs(calc_saturated_vapour_concentration_air(temp_air) - ref_concentration) < 0.1

def test_calc_vapour_concentration_air():
    """Values taken from: https://www.engineeringtoolbox.com/relative-humidity-air-d_687.html"""
    temp_air = 20.0
    relative_humidity = 57.8
    ref_concentration = 10
    assert abs(calc_vapour_concentration_air(temp_air, relative_humidity) - ref_concentration) < 0.1


def test_vapour_concentration_deficit():
    """Data from experiment 3 from table 1 of paper"""
    temp_air = 21
    relative_humidity = 76
    vapour_concentration_deficit = 4.4
    assert abs(calc_vapour_concentration_deficit(temp_air, relative_humidity) - vapour_concentration_deficit) < 0.1


def test_stomatal_resistance():
    """Data from experiment 1(A) from table 1 of paper """
    ppfd = 140
    ref_stomatal_resistance = 289
    assert abs(calc_stomatal_resistance(ppfd) - ref_stomatal_resistance) < 1


def test_stomatal_resistance_null():
    ppfd = 0
    ref_stomatal_resistance = 450
    assert abs(calc_stomatal_resistance(ppfd) - ref_stomatal_resistance) < 1


def test_net_radiation():
    """Data from set A from table 2 of paper"""
    ppfd = 600
    reflection_coefficient = 0.05
    cultivation_area_coverage = 0.95
    net_radiation = calc_net_radiation(ppfd, reflection_coefficient, cultivation_area_coverage)
    assert abs(net_radiation - 108.3) < 0.1

# def test_sensible_heat_exchange():
#     # NEEDS CHECKING
#     temp_air = 21.0
#     temp_surface = 23.0
#     lai = 3.0
#     vapour_resistance = 100
#     ref_value =  0.06
#     assert(abs(calc_sensible_heat_exchange(temp_air, temp_surface, lai, vapour_resistance) - ref_value) < 0.01)

def test_epsilon():
    """
    Data from table 1 (experiment 1) ppfd of 600
    Following data from the paper:

    Air temperature (Ta): 21C
    LAI: 3.0
    Aerodynamic boundary layer resistance (R_a): 100 s m-1
    Stomatal resistance (R_s): 158 s m-1

    Also need the following constants, for which I’ve taken standard SI values:
    Heat capacity of air (Rc): 1003 J  kg-1 C-1

    From figure 3 it looks like the transpiration rate for 600 ppfd is approximately: 0.051g m-2 s

    From the FAO page on Crop evapotranspiration (http://www.fao.org/3/X0490E/x0490e0i.htm) it looks like you multiply by 2450
    to go from g m-2 s to W m-2, so that makes the transpiration rate for 600 ppfd approximately: 125 W m-2 (this seems roughly
    correct to me).

    If I assume that the temperature difference is of the order of a few degrees (say Ts = 23C) and Ts - Ta = 2.0
    (again, it’s just to show an order of magnitude), then:


    Equation 6 gives the latent heat exchange (and therefore the transpiration rate) as:
    La_E = LAI * Lambda * ((X_sts - X_sta) / (R_s + R_a))

    Equation 7 to calculate the vapour concentration at the surface is:
    X_sts = X_sta + (Rc / Lambda) * Epsilon * (Ts - Ta)

    Rearranging to calculate Epsilon gives me:
    Epsilon = (La_E * (R_s + R_a) ) / (LAI * Rc * (Ts - Ta))

    Epsilon = 7.0988361358001234e-06

    """
    temp_air = 23.0
    epsilon = calc_epsilon(temp_air)
    ref = 2.58659e-06
    assert(abs(epsilon - ref) < 1.0e-08 )


def test_latent_heat_flux():
    # NEEDS CHECKING
    pass



if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
