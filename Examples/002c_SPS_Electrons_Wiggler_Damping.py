#import sys
#import os
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from matplotlib import pyplot as plt

import xtrack as xt

from synchrotron_integrals import SynchrotronIntegral as synint
from wiggler import wiggler

line = xt.Line.from_json('001_sps.json')
line.particle_ref = xt.Particles(energy0=20e9, mass0=xt.ELECTRON_MASS_EV)


# Notes from figures:
# NOTE: The damping factors are almost independent of the angle of the wiggler.
# NOTE: The relative error is not simply quadratic.
# NOTE: The discrepancy between the two angles is smaller for the XSuite method than for the Integral one.


# Wiggler
# ----------------------------------------------------------------------------------------------------------------------------------------------
# Gets the position where we can place the wiggler.
tt = line.get_table()
s_start_wig = tt['s', 'actcsg.31780']

# Wiggler parameters
k0_wig = 0.000001
angle = 0#np.pi/2

k0_values = np.array([1e-5, 1e-4, 1e-3, 1e-2, 1e-1])

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
Alpha_EI0 = np.zeros(iters)
Alpha_EX0 = np.zeros(iters)
Alpha_ERelEr0 = np.zeros(iters)

print(f'Damping Factor for angle = {angle} [rad]:')
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

    Alpha_EI0[iter] = Integrals.radiation_damping_s()[2]
    Alpha_EX0[iter] = tw_rad.damping_constants_s[2]
    Alpha_ERelEr0[iter] = np.abs(Alpha_EI0[iter]/Alpha_EX0[iter] - 1)

#tw_rad.loglog('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'Alpha_EI = {Alpha_EI0} [s⁻¹]')
print(f'Alpha_EX = {Alpha_EX0} [s⁻¹]')
print(f'Relative Error = {Alpha_ERelEr0}')

# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = pi/2 rad.
angle = np.pi/2

Alpha_EI90 = np.zeros(iters)
Alpha_EX90 = np.zeros(iters)
Alpha_ERelEr90 = np.zeros(iters)

for number in range(len(Chicane)):
        line[list(Chicane.keys())[number]].rot_s_rad = angle


print(f'Damping Factor for angle = {angle} [rad]:')
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

    Alpha_EI90[iter] = Integrals.radiation_damping_s()[2]
    Alpha_EX90[iter] = tw_rad.damping_constants_s[2]
    Alpha_ERelEr90[iter] = np.abs(Alpha_EI90[iter]/Alpha_EX90[iter] - 1)

#tw_rad.loglog('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'Alpha_EI = {Alpha_EI90} [s⁻¹]')
print(f'Alpha_EX = {Alpha_EX90} [s⁻¹]')
print(f'Relative Error = {Alpha_ERelEr90}')


# -------------------------------------------------------------------------------------------------------------------------------
# Baseline: k_wig = 0
print(f'Damping Factor for k_wig = 0 [m⁻¹]:')
line.discard_tracker()
for number in range(len(Chicane)):
    line[list(Chicane.keys())[number]].k0 = 0
    
print(f'k0 = {line[list(Chicane.keys())[0]].k0} [m⁻¹]')
    
line.build_tracker()

line.configure_radiation(model='mean')

tw_rad = line.twiss(strengths=True, eneloss_and_damping=True)

#tw_rad.loglog('x y', 'dx dy')
#plt.show()

#print(tw_rad.keys())

# Compute radiation integrals
Integrals = synint(line)

Alpha_EIbase = Integrals.radiation_damping_s()[2]
Alpha_EXbase = tw_rad.damping_constants_s[2]
Alpha_ERelErbase = np.abs(Alpha_EIbase/Alpha_EXbase - 1)

print(f'Alpha_EI = {Alpha_EIbase} [s⁻¹]')
print(f'Alpha_EX = {Alpha_EXbase} [s⁻¹]')
print(f'Relative Error = {Alpha_ERelErbase}')

line.discard_tracker()


# Fitting
# -------------------------------------------------------------------------------------------------------------------------------
# Fit for horizontal wiggler
fitI0 = np.polyfit(k0_values, Alpha_EI0, 2)
fitX0 = np.polyfit(k0_values, Alpha_EX0, 2)
Alpha_EI0fit = fitI0[0]*k0_values**2 + fitI0[1]*k0_values + fitI0[2]
Alpha_EX0fit = fitX0[0]*k0_values**2 + fitX0[1]*k0_values + fitX0[2]

# Fit for vertical wiggler
fitI90 = np.polyfit(k0_values, Alpha_EI90, 2)
fitX90 = np.polyfit(k0_values, Alpha_EX90, 2)
Alpha_EI90fit = fitI90[0]*k0_values**2 + fitI90[1]*k0_values + fitI90[2]
Alpha_EX90fit = fitX90[0]*k0_values**2 + fitX90[1]*k0_values + fitX90[2]

# Relative errors
Alpha_ERelEr0fit = Alpha_EI0fit/Alpha_EX0fit - 1
Alpha_ERelEr90fit = Alpha_EI90fit/Alpha_EX90fit - 1


# Plotting
# -------------------------------------------------------------------------------------------------------------------------------
# Plot momentum compaction
# Create the plot
RelErrs = plt.figure(figsize=(10, 6))

plt.loglog(k0_values, Alpha_ERelEr90, label='Alpha_ERelEr90')
plt.loglog(k0_values, Alpha_ERelEr0, label='Alpha_ERelEr0')

# Plot the fits
plt.loglog(k0_values, Alpha_ERelEr0fit, label='Fit Alpha_ERelEr0')
plt.loglog(k0_values, Alpha_ERelEr90fit, label='Fit Alpha_ERelEr90')

plt.axhline(y=Alpha_ERelErbase, color='gray', linestyle='--', linewidth=1, label='No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('Alpha_E Relative Errors')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Relative Errors in the Damping Factors per Second')

plt.show()
RelErrs.savefig("Figures/Alpha_ERelativeErrors.pdf", bbox_inches='tight')

# loglog relative errors
# Create the loglog
Values = plt.figure(figsize=(10, 6))

# loglog the first set of data
plt.loglog(k0_values, Alpha_EI0, label='Alpha_EI0')
plt.loglog(k0_values, Alpha_EX0, label='Alpha_EX0')

# loglog the second set of data
plt.loglog(k0_values, Alpha_EI90, label='Alpha_EI90')
plt.loglog(k0_values, Alpha_EX90, label='Alpha_EX90')

# Plot the fits
plt.loglog(k0_values, Alpha_EI0fit, label='Fit Alpha_EI0')
plt.loglog(k0_values, Alpha_EX0fit, label='Fit Alpha_EX0')
plt.loglog(k0_values, Alpha_EI90fit, label='Fit Alpha_EI90')
plt.loglog(k0_values, Alpha_EX90fit, label='Fit Alpha_EX90')

plt.axhline(y=Alpha_EIbase, color='gray', linestyle='--', linewidth=1, label='Integral, No Chicane')
plt.axhline(y=Alpha_EXbase, color='black', linestyle='--', linewidth=1, label='XSuite, No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('Alpha_E [s⁻¹]')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Damping Factor for Two Angles')

# Show the loglog
plt.show()
Values.savefig("Figures/Alpha_EValues.pdf", bbox_inches='tight')
