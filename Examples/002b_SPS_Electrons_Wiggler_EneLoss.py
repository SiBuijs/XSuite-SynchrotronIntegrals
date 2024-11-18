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

# Wiggler
# ----------------------------------------------------------------------------------------------------------------------------------------------
# Gets the position where we can place the wiggler.
tt = line.get_table()
s_start_wig = tt['s', 'actcsg.31780']

# Wiggler parameters
k0_wig = 0.000001
angle = 0#np.pi/2

k0_values =np.array([1e-5, 1e-4, 1e-3, 1e-2, 1e-1]) #np.array([1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1])

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
U_0I0 = np.zeros(iters)
U_0X0 = np.zeros(iters)
U_0RelEr0 = np.zeros(iters)

print(f'Energy Loss per Turn for angle = {angle} [rad]:')
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

    U_0I0[iter] = Integrals.energy_loss()
    U_0X0[iter] = tw_rad.eneloss_turn
    U_0RelEr0[iter] = np.abs(U_0I0[iter]/U_0X0[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'U_0I = {U_0I0} [eV]')
print(f'U_0X = {U_0X0} [eV]')
print(f'Relative Error = {U_0RelEr0}')

# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = pi/2 rad.
angle = np.pi/2

U_0I90 = np.zeros(iters)
U_0X90 = np.zeros(iters)
U_0RelEr90 = np.zeros(iters)

for number in range(len(Chicane)):
        line[list(Chicane.keys())[number]].rot_s_rad = angle


print(f'Energy Loss per Turn for angle = {angle} [rad]:')
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

    U_0I90[iter] = Integrals.energy_loss()
    U_0X90[iter] = tw_rad.eneloss_turn
    U_0RelEr90[iter] = np.abs(U_0I90[iter]/U_0X90[iter] - 1)

#tw_rad.plot('x y', 'dx dy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'U_0I = {U_0I90} [eV]')
print(f'U_0X = {U_0X90} [eV]')
print(f'Relative Error = {U_0RelEr90}')


# -------------------------------------------------------------------------------------------------------------------------------
# Baseline: k_wig = 0
print(f'Energy Loss per Turn for k_wig = 0 [m⁻¹]:')
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

U_0Ibase = Integrals.energy_loss()
U_0Xbase = tw_rad.eneloss_turn
U_0RelErbase = np.abs(U_0Ibase/U_0Xbase - 1)

print(f'U_0I = {U_0Ibase} [eV]')
print(f'U_0X = {U_0Xbase} [eV]')
print(f'Relative Error = {U_0RelErbase}')

line.discard_tracker()


# Rescale
# -------------------------------------------------------------------------------------------------------------------------------
#U_0I0 = U_0I0 - U_0Ibase
#U_0X0 = U_0X0 - U_0Xbase
#U_0I90 = U_0I90 - U_0Ibase
#U_0X90 = U_0X90 - U_0Xbase

# Fitting
# -------------------------------------------------------------------------------------------------------------------------------
# Fit for horizontal wiggler
fitI0 = np.polyfit(k0_values, U_0I0, 2)
fitX0 = np.polyfit(k0_values, U_0X0, 2)
U_0I0fit = fitI0[0]*k0_values**2 + fitI0[1]*k0_values + fitI0[2]
U_0X0fit = fitX0[0]*k0_values**2 + fitX0[1]*k0_values + fitX0[2]

# Print the formulas
print(f'U_0I0fit = {fitI0[0]}*k0_values**2 + {fitI0[1]}*k0_values + {fitI0[2]}')
print(f'U_0X0fit = {fitX0[0]}*k0_values**2 + {fitX0[1]}*k0_values + {fitX0[2]}')

# Fit for vertical wiggler
fitI90 = np.polyfit(k0_values, U_0I90, 2)
fitX90 = np.polyfit(k0_values, U_0X90, 2)
U_0I90fit = fitI90[0]*k0_values**2 + fitI90[1]*k0_values + fitI90[2]
U_0X90fit = fitX90[0]*k0_values**2 + fitX90[1]*k0_values + fitX90[2]

# Print the formulas
print(f'U_0I90fit = {fitI90[0]}*k0_values**2 + {fitI90[1]}*k0_values + {fitI90[2]}')
print(f'U_0X90fit = {fitX90[0]}*k0_values**2 + {fitX90[1]}*k0_values + {fitX90[2]}')

# Relative errors
U_0RelEr0fit = U_0I0fit/U_0X0fit - 1
U_0RelEr90fit = U_0I90fit/U_0X90fit - 1

# Plotting
# -------------------------------------------------------------------------------------------------------------------------------
# Create the plot
RelErrs = plt.figure(figsize=(10, 6))

plt.loglog(k0_values, U_0RelEr90, label='U_0RelEr90')
plt.loglog(k0_values, U_0RelEr0, label='U_0RelEr0')

# Plot the fits
plt.loglog(k0_values, U_0RelEr0fit, label='Fit U_0RelEr0')
plt.loglog(k0_values, U_0RelEr90fit, label='Fit U_0RelEr90')

plt.axhline(y=U_0RelErbase, color='gray', linestyle='--', linewidth=1, label='No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('U_0 Relative Errors')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Relative Errors in the Energy Loss per Turn')

plt.show()
RelErrs.savefig("Figures/U_0RelativeErrors.pdf", bbox_inches='tight')

# Plot relative errors
# Create the plot
Values = plt.figure(figsize=(10, 6))

# Plot the first set of data
plt.loglog(k0_values, U_0I0, label='U_0I0')
plt.loglog(k0_values, U_0X0, label='U_0X0')

# Plot the second set of data
plt.loglog(k0_values, U_0I90, label='U_0I90')
plt.loglog(k0_values, U_0X90, label='U_0X90')

# Plot the fits
plt.loglog(k0_values, U_0I0fit, label='Fit U_0I0')
plt.loglog(k0_values, U_0X0fit, label='Fit U_0X0')
plt.loglog(k0_values, U_0I90fit, label='Fit U_0I90')
plt.loglog(k0_values, U_0X90fit, label='Fit U_0X90')

plt.axhline(y=U_0Ibase, color='gray', linestyle='--', linewidth=1, label='Integral, No Chicane')
plt.axhline(y=U_0Xbase, color='black', linestyle='--', linewidth=1, label='XSuite, No Chicane')

# Set x-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
plt.xlabel('k_wig [m⁻¹]')
plt.ylabel('U_0 [eV]')
plt.legend(loc='best')  # Automatically place the legend in an optimal location
plt.title('Energy Loss per Turn for Two Angles')

# Show the plot
plt.show()
Values.savefig("Figures/U_0Values.pdf", bbox_inches='tight')
