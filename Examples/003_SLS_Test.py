# Basic SLS2.0 Simulation in XSuite.
# Date: 2024-09-25
#
# This simulation imports a .madx file and converts it into an XSuite line.
# Many lines of this code are copied from the XSuite User's Guide, so they may need adjustments specific to the SLS2.0's needs.

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# Packages ------------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

import xtrack as xt
import xobjects as xo

from scipy.constants import c as clight

# Import the .madx file and put it on a line --------------------------------------------------------------------------------------------------------
# Create a line object and import the SPS lattice from the sls2 sequence.
#
# NOTE: S.F. Buijsman, 2024-09-25
# It should be noted that, in order to make the lines below work, I needed to
# replace "hkicker" and "vkicker" with "kicker" in 'b075_2024.08.26_final_seq.madx'.
# The file is now called 'b075_2024.09.25.madx'.
# The original file (with "hkicker" and "vkicker") is stored as 'b075_2024.08.26_final_seq.madx'.
# Frankly, I do not quite know what the "hkicker" and "vkicker" where doing there, because they had no assigned "kick" value;
# Most of them were not assigned any value, some only had a length.
# I thought they might be important in beam correction later on, that is why I kept them.
#
# TODO: Ask Simona and Michael about the above.
line = xt.Line.from_json('../Example_Data/b075_2024.09.25.json')
line.particle_ref = xt.Particles(q0 = -1, mass0 = xt.ELECTRON_MASS_EV, kinetic_energy0=2.7e9)


# Frequency and voltage of the cavities from the technical design report.
cavfreq  = 499.6537e6
cavvolt  = 445e3


# Harmonic number, for later use.
harmonic = 480


# Create a table to inspect elements if needed.
tab = line.get_table()
print(tab.rows[190:210:'s'])


# This is the place where arc 04 ends and long straight 05 starts.
# The cavities should be located somewhere here.
# The drift region is about 4 meters long, so somewhere between s_start_long05 and s_start_long05 + 4.
s_start_long05 = tab['s', 'ars05_gmrk_0000']
print(f's_start_long05 = {s_start_long05} [m]')
s_start_maincavity = s_start_long05 + 1

# Insert main cavities.
line.insert_element('maincav01', at_s=s_start_maincavity, element=xt.Cavity(frequency=cavfreq, voltage=4*cavvolt, lag=180))

# Find the start of the secondary cavities
#s_start_long09 = tab['s', 'ars09_gmrk_0000']
#print(f's_start_long09 = {s_start_long09} [m]')
#s_start_3rdcavity = s_start_long09 + 1

# Insert secondary cavities.
#line.insert_element('3rdcav01', at_s=s_start_3rdcavity, element=xt.Cavity(frequency=1498.95e6, voltage=1e6))

slicing_strategies = [
    xt.Strategy(slicing=xt.Teapot(8)),  # Default
    xt.Strategy(slicing=xt.Teapot(2), element_type=xt.Bend),
    xt.Strategy(slicing=xt.Teapot(8), element_type=xt.Quadrupole),
    xt.Strategy(slicing=xt.Teapot(8), element_type=xt.Multipole),
]
line.slice_thick_elements(slicing_strategies)

# Build the tracker.
line.build_tracker(_context=xo.ContextCpu())
line.configure_radiation(model='mean')


# Tracking ------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the Twiss parameters.
tw = line.twiss(method='6d', strengths=True, eneloss_and_damping=True)


from synchrotron_integrals import SynchrotronIntegral as radint

Integrals = radint(line)


# Make a list of particles.
#n_part = 200
#particles = line.build_particles(
#                        x=np.random.uniform(-1e-3, 1e-3, n_part),
#                        px=np.random.uniform(-1e-5, 1e-5, n_part),
#                        y=np.random.uniform(-2e-3, 2e-3, n_part),
#                        py=np.random.uniform(-3e-5, 3e-5, n_part),
#                        zeta=np.random.uniform(-1e-2, 1e-2, n_part),
#                        delta=np.random.uniform(-1e-4, 1e-4, n_part))


# Track particles over 100 turns.
# Tracking time ~8.7 seconds for 200 particles and 100 turns on my desktop.
#num_turns = 100
#line.track(particles, num_turns=num_turns, time=True)
#print(f'Tracking time: {line.time_last_track}')



# Plotting ------------------------------------------------------------------------------------------------------------------------------------------
# Plot the Twiss parameters.
plt.close('all')


# Plot some Twiss parameters.
tw.plot('betx bety', 'dx dy')
plt.show()


# Survey the line and draw a FloorPlot.
# Commented this out, because I don't need it every time.
# Feel free to uncomment it if you want to see the floor plan.

sv = line.survey()
xplt.FloorPlot(sv, line)
plt.legend(fontsize='small', loc='upper left')
plt.show()


# Calculate associated physical quantities
# Momentum compaction factor
# For this case, the momentum compaction is negative. At first I was worried about this, but it is allowed to be negative.
# A negative momentum compaction factor means that particles of higher energy follow a longer path around the ring.
# This is because the particles follow a longer path than the design orbit.
# If gamma is great enough, the path length dominates the velocity.
# See Wiedemann page 248-250.
# Momentum Compaction Factor
# NOTE: CORRECT
alpha_cI = Integrals.momentum_compaction()
alpha_cX = tw.momentum_compaction_factor
alpha_cTDR = 1.05e-4
print("\n-----------------------------------------------------------------------------------")
print(f"Momentum Compaction Factors:\nalpha_cIx = {alpha_cI} [-]")
print(f"alpha_xs = {alpha_cX} [-]\n")
print(f"alpha_CDRp17 = {alpha_cTDR} [-]\n")
print(f"alpha_xs/alpha_cIx - 1 = {alpha_cI/alpha_cX - 1:2e}")

# Energy loss
# NOTE: CORRECT
U_0I = Integrals.energy_loss()
U_0X = tw.eneloss_turn
U_0TDR = 688e3
print("\n-----------------------------------------------------------------------------------")
print(f"Energy Loss:\nU_0I = {U_0I} [eV]")
print(f"U_0X = {U_0X} [eV]\n")
print(f"U_0TDRp17 = {U_0TDR} [eV]\n")
print(f"U_0I/U_0X - 1 = {U_0I/U_0X - 1:2e}")

# Radiation damping [s⁻¹]
# NOTE: CORRECT
alpha_xIs = Integrals.radiation_damping_s()[0]
alpha_yIs = Integrals.radiation_damping_s()[1]
alpha_sIs = Integrals.radiation_damping_s()[2]

alpha_xXs = tw.damping_constants_s[0]
alpha_yXs = tw.damping_constants_s[1]
alpha_sXs = tw.damping_constants_s[2]

alpha_xTDR = 1/4.2e-3
alpha_yTDR = 1/7.8e-3
alpha_sTDR = 1/6.8e-3

RelErr_xs = alpha_xIs/alpha_xXs - 1
RelErr_ys = alpha_yIs/alpha_yXs - 1
RelErr_ss = alpha_sIs/alpha_sXs - 1

print("\n-----------------------------------------------------------------------------------")
print(f"Radiation Damping (Integral Method):\nalpha_x = {alpha_xIs} [s⁻¹]\nalpha_y = {alpha_yIs} [s⁻¹]\nalpha_s = {alpha_sIs} [s⁻¹]\n")
print(f"Radiation Damping (XSuite Method):\nalpha_x = {alpha_xXs} [s⁻¹]\nalpha_y = {alpha_yXs} [s⁻¹]\nalpha_s = {alpha_sXs} [s⁻¹]\n")
print(f"Radiation Damping (TDRp23):\nalpha_x = {alpha_xTDR} [s⁻¹]\nalpha_y = {alpha_yTDR} [s⁻¹]\nalpha_s = {alpha_sTDR} [s⁻¹]\n")
print(f"Relative Errors:\nalpha_xI/alpha_xX - 1 = {RelErr_xs:2e}\nalpha_xI/alpha_xX - 1 = {RelErr_ys:2e}\nalpha_sI/alpha_sX - 1 = {RelErr_ss:2e}\n")

# Radiation damping [turns⁻¹]
# NOTE: CORRECT
alpha_xIt = Integrals.radiation_damping_turns()[0]
alpha_yIt = Integrals.radiation_damping_turns()[1]
alpha_sIt = Integrals.radiation_damping_turns()[2]

alpha_xXt = tw.damping_constants_turns[0]
alpha_yXt = tw.damping_constants_turns[1]
alpha_sXt = tw.damping_constants_turns[2]

RelErr_xt = alpha_xIt/alpha_xXt - 1
RelErr_yt = alpha_yIt/alpha_yXt - 1
RelErr_st = alpha_sIt/alpha_sXt - 1

print("\n-----------------------------------------------------------------------------------")
print(f"Radiation Damping (Integral Method):\nalpha_x = {alpha_xIt} [turns⁻¹]\nalpha_y = {alpha_yIt} [turns⁻¹]\nalpha_s = {alpha_sIt} [turns⁻¹]\n")
print(f"Radiation Damping (XSuite Method):\nalpha_x = {alpha_xXt} [turns⁻¹]\nalpha_y = {alpha_yXt} [turns⁻¹]\nalpha_s = {alpha_sXt} [turns⁻¹]\n")
print(f"Radiation Damping (TDRp23):\nalpha_x = {alpha_xTDR*tw['T_rev0']} [turns⁻¹]\nalpha_y = {alpha_yTDR*tw['T_rev0']} [turns⁻¹]\nalpha_s = {alpha_sTDR*tw['T_rev0']} [turns⁻¹]\n")
print(f"Relative Errors:\nalpha_xI/alpha_xX - 1 = {RelErr_xt:2e}\nalpha_xI/alpha_xX - 1 = {RelErr_yt:2e}\nalpha_sI/alpha_sX - 1 = {RelErr_st:2e}\n")


# Quantum excitation in the longitudinal direction
eq_emittanceI = Integrals.equilibrium_emittance()
eq_emittanceX = tw.eq_gemitt_zeta
eq_emittanceTDR = 1.103e-3
print("\n-----------------------------------------------------------------------------------")
print(f"eq_emittanceI = {eq_emittanceI} [m]")
print(f"eq_emittanceX = {eq_emittanceX} [m]\n")
print(f"eq_emittanceTDRp23 = {eq_emittanceTDR} [m]\n")
print(f"eq_emittanceI/eq_emittanceX - 1 = {eq_emittanceI/eq_emittanceX - 1:2e}")

# Quantum excitation in the radial direction
# NOTE: CORRECT
rms_betatronIx = Integrals.rms_betatron()[0]
rms_betatronIy = Integrals.rms_betatron()[1]
rms_betatronXx = tw.eq_gemitt_x
rms_betatronXy = tw.eq_gemitt_y
print("\n-----------------------------------------------------------------------------------")
print(f"rms_betatronIx = {rms_betatronIx} [-]")
print(f"rms_betatronIy = {rms_betatronIy} [-]")
print(f"rms_betatronX = {rms_betatronXx} [-]")
print(f"rms_betatronXy = {rms_betatronXy} [-]\n")
print(f"rms_betatronIx/rms_betatronX - 1 = {rms_betatronIx/rms_betatronXx - 1:2e}")
print(f"rms_betatronIy/rms_betatronY - 1 = {rms_betatronIy/rms_betatronXy - 1:2e}")

# This is here so I can place a debug break point at the end, so that I can still test some things.
print('Done')