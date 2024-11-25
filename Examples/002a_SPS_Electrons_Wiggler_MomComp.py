import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from matplotlib import pyplot as plt

import xtrack as xt

from synchrotron_integrals import SynchrotronIntegral as synint
from wiggler import wiggler

line = xt.Line.from_json('Example_Data/001_sps.json')
line.particle_ref = xt.Particles(energy0=20e9, mass0=xt.ELECTRON_MASS_EV)


# Notes from figures:
# NOTE: The momentum compaction factor as calculated using the integral method is independent on the strength of the wiggler.
# NOTE: There is a slight deviation in the momentum compaction factors, depending on the angle of the wiggler.
# NOTE: The discrepancy between the two methods is on the order of 1e-2.


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
alpha_cI0 = np.zeros(iters)
alpha_cX0 = np.zeros(iters)
alpha_cRelEr0 = np.zeros(iters)

print(f'Momentum Compaction Factors for angle = {angle} [rad]:')
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
    #print(tw_rad.keys())

    Integrals = synint(line)

    alpha_cI0[iter] = Integrals.momentum_compaction()
    alpha_cX0[iter] = tw_rad.momentum_compaction_factor
    alpha_cRelEr0[iter] = np.abs(alpha_cI0[iter]/alpha_cX0[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'alpha_cI = {alpha_cI0}')
print(f'alpha_cX = {alpha_cX0}')
print(f'Relative Error = {alpha_cRelEr0}')

# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = pi/2 rad.
angle = np.pi/2

alpha_cI90 = np.zeros(iters)
alpha_cX90 = np.zeros(iters)
alpha_cRelEr90 = np.zeros(iters)

for number in range(len(Chicane)):
        line[list(Chicane.keys())[number]].rot_s_rad = angle


print(f'Momentum Compaction Factors for angle = {angle} [rad]:')
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

    alpha_cI90[iter] = Integrals.momentum_compaction()
    alpha_cX90[iter] = tw_rad.momentum_compaction_factor
    alpha_cRelEr90[iter] = np.abs(alpha_cI90[iter]/alpha_cX90[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'alpha_cI = {alpha_cI90}')
print(f'alpha_cX = {alpha_cX90}')
print(f'Relative Error = {alpha_cRelEr90}')
# NOTE: Test with one active pole.

# -------------------------------------------------------------------------------------------------------------------------------
# Baseline: k_wig = 0
print(f'Momentum Compaction Factors for k_wig = 0 [m⁻¹]:')
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

alpha_cIbase = Integrals.momentum_compaction()
alpha_cXbase = tw_rad.momentum_compaction_factor
alpha_cRelErbase = np.abs(alpha_cIbase/alpha_cXbase - 1)

print(f'alpha_cI = {alpha_cIbase}')
print(f'alpha_cX = {alpha_cXbase}')
print(f'Relative Error = {alpha_cRelErbase}')

line.discard_tracker()


# Fitting
# -------------------------------------------------------------------------------------------------------------------------------
# Fit for horizontal wiggler
fitI0 = np.polyfit(k0_values, alpha_cI0, 2)
fitX0 = np.polyfit(k0_values, alpha_cX0, 2)
alpha_cI0fit = fitI0[0]*k0_values**2 + fitI0[1]*k0_values + fitI0[2]
alpha_cX0fit = fitX0[0]*k0_values**2 + fitX0[1]*k0_values + fitX0[2]

print(f'alpha_cI0fit = {fitI0[0]}*k0² + {fitI0[1]}*k0 + {fitI0[2]}')
print(f'alpha_cX0fit = {fitX0[0]}*k0² + {fitX0[1]}*k0 + {fitX0[2]}')

# Fit for vertical wiggler
fitI90 = np.polyfit(k0_values, alpha_cI90, 2)
fitX90 = np.polyfit(k0_values, alpha_cX90, 2)
alpha_cI90fit = fitI90[0]*k0_values**2 + fitI90[1]*k0_values + fitI90[2]
alpha_cX90fit = fitX90[0]*k0_values**2 + fitX90[1]*k0_values + fitX90[2]

print(f'alpha_cI90fit = {fitI90[0]}*k0² + {fitI90[1]}*k0 + {fitI90[2]}')
print(f'alpha_cX90fit = {fitX90[0]}*k0² + {fitX90[1]}*k0 + {fitX90[2]}')

# Relative errors
alpha_cRelEr0fit = alpha_cI0fit/alpha_cX0fit - 1
alpha_cRelEr90fit = alpha_cI90fit/alpha_cX90fit - 1

# Plotting
# -------------------------------------------------------------------------------------------------------------------------------
# Plot momentum compaction
# Create the plot
RelErrs = plt.figure(figsize=(10, 6))

plt.plot(k0_values, alpha_cRelEr0, label='alpha_cRelEr0')
plt.plot(k0_values, alpha_cRelEr90, label='alpha_cRelEr90')
plt.plot(k0_values, alpha_cRelEr0fit, label='Fit alpha_cRelEr0')
plt.plot(k0_values, alpha_cRelEr90fit, label='Fit alpha_cRelEr90')

plt.axhline(y=alpha_cRelErbase, color='gray', linestyle='--', linewidth=1, label='No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('alpha_c Relative Errors')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Relative Errors in the Momentum Compaction Factor')

plt.show()
RelErrs.savefig("Figures/alpha_cRelativeErrors.pdf", bbox_inches='tight')

# Plot relative errors
# Create the plot
Values = plt.figure(figsize=(10, 6))

# Plot the first set of data
plt.plot(k0_values, alpha_cI0, label='alpha_cI0')
plt.plot(k0_values, alpha_cX0, label='alpha_cX0')
plt.plot(k0_values, alpha_cI0fit, label='Fit alpha_cI0')
plt.plot(k0_values, alpha_cX0fit, label='Fit alpha_cX0')

# Plot the second set of data
plt.plot(k0_values, alpha_cI90, label='alpha_cI90')
plt.plot(k0_values, alpha_cX90, label='alpha_cX90')
plt.plot(k0_values, alpha_cI90fit, label='Fit alpha_cI90')
plt.plot(k0_values, alpha_cX90fit, label='Fit alpha_cX90')

plt.axhline(y=alpha_cIbase, color='gray', linestyle='--', linewidth=1, label='Integral, No Chicane')
plt.axhline(y=alpha_cXbase, color='black', linestyle='--', linewidth=1, label='XSuite, No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('alpha_c')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Momentum Compaction Factors for Two Angles')

# Show the plot
plt.show()
Values.savefig("Figures/alpha_cValues.pdf", bbox_inches='tight')
