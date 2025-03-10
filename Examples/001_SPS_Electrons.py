import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import xtrack as xt

from synchrotron_integrals import SynchrotronIntegral as synint

import numpy as np

from scipy.constants import c as clight

line = xt.Line.from_json('Example_Data/001_sps.json')
line.particle_ref = xt.Particles(energy0=20e9, mass0=xt.ELECTRON_MASS_EV)

line.discard_tracker()
slicing_strategies = [
    xt.Strategy(slicing=xt.Teapot(1)),  # Default
    xt.Strategy(slicing=xt.Teapot(2), element_type=xt.Bend),
    xt.Strategy(slicing=xt.Teapot(8), element_type=xt.Quadrupole),
]
line.slice_thick_elements(slicing_strategies)

tw_no_rad = line.twiss(strengths=True)

line.configure_radiation(model='mean')

tw_rad = line.twiss(method='6d', strengths=True, eneloss_and_damping=True, radiation_integrals=True)

print(tw_rad.keys())

# Compute radiation integrals
Integrals = synint(line)

# Momentum Compaction Factor
# NOTE: CORRECT
alpha_cI = Integrals.momentum_compaction()
alpha_cX = tw_rad.momentum_compaction_factor
print("\n-----------------------------------------------------------------------------------")
print(f"Momentum Compaction Factors:\nalpha_cIx = {alpha_cI} [-]")
print(f"alpha_xs = {alpha_cX} [-]\n")
print(f"alpha_xs/alpha_cIx - 1 = {alpha_cI/alpha_cX - 1:2e}")

# Energy loss
# NOTE: CORRECT
U_0I = Integrals.energy_loss()
U_0X = tw_rad.eneloss_turn

print("\n-----------------------------------------------------------------------------------")
print(f"Energy Loss:\nU_0I = {U_0I} [eV]")
print(f"U_0X = {U_0X} [eV]\n")
print(f"U_0I/U_0X - 1 = {U_0I/U_0X - 1:2e}")

# Radiation damping [s⁻¹]
# NOTE: CORRECT
alpha_xIs = tw_rad.rad_int_damping_constant_x_s
alpha_yIs = tw_rad.rad_int_damping_constant_y_s
alpha_sIs = tw_rad.rad_int_damping_constant_zeta_s

alpha_xXs = tw_rad.damping_constants_s[0]
alpha_yXs = tw_rad.damping_constants_s[1]
alpha_sXs = tw_rad.damping_constants_s[2]

RelErr_xs = alpha_xIs/alpha_xXs - 1
RelErr_ys = alpha_yIs/alpha_yXs - 1
RelErr_ss = alpha_sIs/alpha_sXs - 1

print("\n-----------------------------------------------------------------------------------")
print(f"Radiation Damping (Integral Method):\nalpha_x = {alpha_xIs} [s⁻¹]\nalpha_y = {alpha_yIs} [s⁻¹]\nalpha_s = {alpha_sIs} [s⁻¹]\n")
print(f"Radiation Damping (XSuite Method):\nalpha_x = {alpha_xXs} [s⁻¹]\nalpha_y = {alpha_yXs} [s⁻¹]\nalpha_s = {alpha_sXs} [s⁻¹]\n")
print(f"Relative Errors:\nalpha_xI/alpha_xX - 1 = {RelErr_xs:2e}\nalpha_xI/alpha_xX - 1 = {RelErr_ys:2e}\nalpha_sI/alpha_sX - 1 = {RelErr_ss:2e}\n")

# Radiation damping [turns⁻¹]
# NOTE: CORRECT
alpha_xIt = Integrals.radiation_damping_turns()[0]
alpha_yIt = Integrals.radiation_damping_turns()[1]
alpha_sIt = Integrals.radiation_damping_turns()[2]

alpha_xXt = tw_rad.damping_constants_turns[0]
alpha_yXt = tw_rad.damping_constants_turns[1]
alpha_sXt = tw_rad.damping_constants_turns[2]

RelErr_xt = alpha_xIt/alpha_xXt - 1
RelErr_yt = alpha_yIt/alpha_yXt - 1
RelErr_st = alpha_sIt/alpha_sXt - 1

print("\n-----------------------------------------------------------------------------------")
print(f"Radiation Damping (Integral Method):\nalpha_x = {alpha_xIt} [s⁻¹]\nalpha_y = {alpha_yIt} [s⁻¹]\nalpha_s = {alpha_sIt} [turns⁻¹]\n")
print(f"Radiation Damping (XSuite Method):\nalpha_x = {alpha_xXt} [s⁻¹]\nalpha_y = {alpha_yXt} [s⁻¹]\nalpha_s = {alpha_sXt} [turns⁻¹]\n")
print(f"Relative Errors:\nalpha_xI/alpha_xX - 1 = {RelErr_xt:2e}\nalpha_xI/alpha_xX - 1 = {RelErr_yt:2e}\nalpha_sI/alpha_sX - 1 = {RelErr_st:2e}\n")

# Quantum excitation in the longitudinal direction
eq_emittanceI = Integrals.equilibrium_emittance()
eq_emittanceX = tw_rad.eq_gemitt_zeta
print("\n-----------------------------------------------------------------------------------")
print(f"eq_emittanceI = {eq_emittanceI} [m]")
print(f"eq_emittanceX = {eq_emittanceX} [m]\n")
print(f"eq_emittanceI/eq_emittanceX - 1 = {eq_emittanceI/eq_emittanceX - 1:2e}")

# Quantum excitation in the radial direction
# NOTE: CORRECT
rms_betatronIx = tw_rad.rad_int_eq_gemitt_x
rms_betatronIy = tw_rad.rad_int_eq_gemitt_y
rms_betatronXx = tw_rad.eq_gemitt_x
rms_betatronXy = tw_rad.eq_gemitt_y
print("\n-----------------------------------------------------------------------------------")
print(f"rms_betatronIx = {rms_betatronIx} [-]")
print(f"rms_betatronIy = {rms_betatronIy} [-]")
print(f"rms_betatronXx = {rms_betatronXx} [-]")
print(f"rms_betatronXy = {rms_betatronXy} [-]\n")
print(f"rms_betatronIx/rms_betatronX - 1 = {rms_betatronIx/rms_betatronXx - 1:2e}")
print(f"rms_betatronIy/rms_betatronY - 1 = {rms_betatronIy/rms_betatronXy - 1:2e}")