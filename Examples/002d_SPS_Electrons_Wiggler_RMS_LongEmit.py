import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from matplotlib import pyplot as plt

import xtrack as xt

from synchrotron_integrals import SynchrotronIntegral as synint
from wiggler import wiggler


# Notes from figures:
# NOTE: There is a small discrepancy between the two methods.
# NOTE: There is a small discrepancy between the two angles as well.


line = xt.Line.from_json('Example_Data/001_sps.json')
line.particle_ref = xt.Particles(energy0=20e9, mass0=xt.ELECTRON_MASS_EV)

# Wiggler
# ----------------------------------------------------------------------------------------------------------------------------------------------
# Gets the position where we can place the wiggler.
tt = line.get_table()
s_start_wig = tt['s', 'actcsg.31780']

# Wiggler parameters
k0_wig = 0.000001
angle = 0#np.pi/2

k0_values = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

lenpole = 0.25
numpoles = 16
lenwig = lenpole * numpoles
numperiods = 4
lambdawig = lenwig / numperiods
rhowig = 1 / (k0_wig + 1e-9)
kwig = 2*np.pi / lambdawig

# Makes a chicane with specified parameters.
ChicaneObject = wiggler(Period=lambdawig, Amplitude=k0_wig, NumPeriods=numperiods, Angle_Rad=angle, Scheme='121a')
# The chicane is a dictionary containing the names, the elements and their positions.
Chicane = ChicaneObject.WigglerDict

# Loop over the elements of Chicane and insert them into the line.
for name, element in Chicane.items():
    line.insert_element(name=name, element=element['element'], at_s=s_start_wig + element['position'])

# Print the elements in te Chicane to check if everything went as desired.
tab = line.get_table()
print(tab.rows[list(Chicane.keys())[0] : list(Chicane.keys())[-1]])



# Slicing
# ----------------------------------------------------------------------------------------------------------------------------------------------
line.discard_tracker()
slicing_strategies = [
    xt.Strategy(slicing=xt.Teapot(1)),  # Default
    xt.Strategy(slicing=xt.Teapot(2), element_type=xt.Bend),
    xt.Strategy(slicing=xt.Teapot(8), element_type=xt.Quadrupole),
    xt.Strategy(slicing=xt.Teapot(20), name='mwp.*'),
]
line.slice_thick_elements(slicing_strategies)


# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = 0 rad.
iters = len(k0_values)
eq_emittanceEI0 = np.zeros(iters)
eq_emittanceEX0 = np.zeros(iters)
eq_emittanceERelEr0 = np.zeros(iters)

print(f'RMS Longitudinal Emittance for angle = {angle} [rad]:')
for iter in range(iters):
    line.discard_tracker()
    for number in range(len(Chicane)):
        sign = np.sign(line[list(Chicane.keys())[number]].k0)
        line[list(Chicane.keys())[number]].k0 = sign*k0_values[iter]

    print(f'k0 = {line[list(Chicane.keys())[0]].k0} [m⁻¹]')
    #print(f'angle = {line[list(Chicane.keys())[0]].rot_s_rad} [rad]')
    
    line.build_tracker()


    line.configure_radiation(model='mean')

    tw_rad = line.twiss(strengths=True, eneloss_and_damping=True)

    Integrals = synint(line)

    eq_emittanceEI0[iter] = Integrals.equilibrium_emittance()
    eq_emittanceEX0[iter] = tw_rad.eq_gemitt_zeta
    eq_emittanceERelEr0[iter] = np.abs(eq_emittanceEI0[iter]/eq_emittanceEX0[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'eq_emittanceEI = {eq_emittanceEI0} [s⁻¹]')
print(f'eq_emittanceEX = {eq_emittanceEX0} [s⁻¹]')
print(f'Relative Error = {eq_emittanceERelEr0}')

# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = pi/2 rad.
angle = np.pi/2

eq_emittanceEI90 = np.zeros(iters)
eq_emittanceEX90 = np.zeros(iters)
eq_emittanceERelEr90 = np.zeros(iters)

for number in range(len(Chicane)):
        line[list(Chicane.keys())[number]].rot_s_rad = angle


print(f'RMS Longitudinal Emittance for angle = {angle} [rad]:')
for iter in range(iters):
    line.discard_tracker()
    for number in range(len(Chicane)):
        sign = np.sign(line[list(Chicane.keys())[number]].k0)
        line[list(Chicane.keys())[number]].k0 = sign*k0_values[iter]
    
    print(f'k0 = {line[list(Chicane.keys())[0]].k0} [m⁻¹]')
    #print(f'angle = {line[list(Chicane.keys())[0]].rot_s_rad} [rad]')

    line.build_tracker()

    line.configure_radiation(model='mean')

    tw_rad = line.twiss(strengths=True, eneloss_and_damping=True)

    Integrals = synint(line)

    eq_emittanceEI90[iter] = Integrals.equilibrium_emittance()
    eq_emittanceEX90[iter] = tw_rad.eq_gemitt_zeta
    eq_emittanceERelEr90[iter] = np.abs(eq_emittanceEI90[iter]/eq_emittanceEX90[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'eq_emittanceEI = {eq_emittanceEI90} [s⁻¹]')
print(f'eq_emittanceEX = {eq_emittanceEX90} [s⁻¹]')
print(f'Relative Error = {eq_emittanceERelEr90}')


# -------------------------------------------------------------------------------------------------------------------------------
# Baseline: k_wig = 0
print(f'RMS Longitudinal Emittance for k_wig = 0 [m⁻¹]:')
line.discard_tracker()
for number in range(len(Chicane)):
    line[list(Chicane.keys())[number]].k0 = 0
    
print(f'k0 = {line[list(Chicane.keys())[0]].k0} [m⁻¹]')
    
line.build_tracker()

line.configure_radiation(model='mean')

tw_rad = line.twiss(strengths=True, eneloss_and_damping=True)

#tw_rad.plot('x y', 'dx dy')
#plt.show()

#print(tw_rad.keys())

# Compute radiation integrals
Integrals = synint(line)

eq_emittanceEIbase = Integrals.equilibrium_emittance()
eq_emittanceEXbase = tw_rad.eq_gemitt_zeta
eq_emittanceERelErbase = np.abs(eq_emittanceEIbase/eq_emittanceEXbase - 1)

print(f'eq_emittanceEI = {eq_emittanceEIbase} [s⁻¹]')
print(f'eq_emittanceEX = {eq_emittanceEXbase} [s⁻¹]')
print(f'Relative Error = {eq_emittanceERelErbase}')

line.discard_tracker()


# Plotting
# -------------------------------------------------------------------------------------------------------------------------------
# Plot momentum compaction
# Create the plot
RelErrs = plt.figure(figsize=(10, 6))

plt.plot(k0_values, eq_emittanceERelEr90, label='eq_emittanceERelEr90')
plt.plot(k0_values, eq_emittanceERelEr0, label='eq_emittanceERelEr0')

plt.axhline(y=eq_emittanceERelErbase, color='gray', linestyle='--', linewidth=1, label='No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('eq_emittanceE Relative Errors')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Relative Errors in the RMS Longitudinal Emittance')

plt.show()
RelErrs.savefig("Figures/RMS_Long_EmitRelativeErrors.pdf", bbox_inches='tight')

# Plot relative errors
# Create the plot
Values = plt.figure(figsize=(10, 6))

# Plot the first set of data
plt.plot(k0_values, eq_emittanceEI0, label='eq_emittanceEI0')
plt.plot(k0_values, eq_emittanceEX0, label='eq_emittanceEX0')

# Plot the second set of data
plt.plot(k0_values, eq_emittanceEI90, label='eq_emittanceEI90')
plt.plot(k0_values, eq_emittanceEX90, label='eq_emittanceEX90')

plt.axhline(y=eq_emittanceEIbase, color='gray', linestyle='--', linewidth=1, label='Integral, No Chicane')
plt.axhline(y=eq_emittanceEXbase, color='black', linestyle='--', linewidth=1, label='XSuite, No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('eq_emittanceE [-]')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('RMS Longitudinal Emittance for Two Angles')

# Show the plot
plt.show()
Values.savefig("Figures/RMS_Long_EmitValues.pdf", bbox_inches='tight')
