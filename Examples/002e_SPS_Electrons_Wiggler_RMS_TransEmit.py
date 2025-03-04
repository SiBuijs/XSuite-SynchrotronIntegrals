import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from matplotlib import pyplot as plt

import xtrack as xt

from synchrotron_integrals import SynchrotronIntegral as synint
from wiggler import wiggler


# NOTE: x-direction: There is a strong discrepancy between the integral and XSuite method for an angle of 90 degrees.
# NOTE: x-direction: This implies some error in the calculation of I5.
# NOTE: y-direction: XSuite and Integral surprisingly agree very well!
# NOTE: y-direction: The vertical emittance for an angle of 0 degrees is approximately zero, as it should be.
# NOTE: y-direction: The vertical emittance from both methods for an angle of 90 degrees agree very well for larger k_wig values.
# NOTE: y-direction: For low values of k_wig, the discrepancy between XSuite and Integrals becomes larger.


line = xt.Line.from_json('Example_Data/001_sps.json')
line.particle_ref = xt.Particles(energy0=20e9, mass0=xt.ELECTRON_MASS_EV)

# Wiggler
# ----------------------------------------------------------------------------------------------------------------------------------------------
# Gets the position where we can place the wiggler.
tt = line.get_table()
s_start_wig = tt['s', 'actcsg.31780']
print(f's_start_wig = {s_start_wig}')

# Wiggler parameters
k0_wig = 0.000001
angle = 0#np.pi/2

k0_values = [5e-4]

lenpole = 0.5
numpoles = 8
lenwig = lenpole * numpoles
numperiods = 2
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
RMS_beta_yI0 = np.zeros(iters)
RMS_beta_yX0 = np.zeros(iters)
RMS_beta_yRelEr0 = np.zeros(iters)

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

    RMS_beta_yI0[iter] = Integrals.rms_betatron()[1]
    RMS_beta_yX0[iter] = tw_rad.eq_gemitt_y
    RMS_beta_yRelEr0[iter] = RMS_beta_yI0[iter]/RMS_beta_yX0[iter] - 1

#tw_rad.plot('dx dy', 'dpx dpy')
#plt.show()
#print(line.get_table(attr=True).cols['name', 'element_type', 'k0l', 'rot_s_rad'].rows['mwp.*'])

print(f'RMS_beta_yI = {RMS_beta_yI0} [s⁻¹]')
print(f'RMS_beta_yX = {RMS_beta_yX0} [s⁻¹]')
print(f'Relative Error = {RMS_beta_yRelEr0}')

# -------------------------------------------------------------------------------------------------------------------------------
# Momentum Compaction Factor for Angle = pi/2 rad.
angle = np.pi/2

RMS_beta_yI90 = np.zeros(iters)
RMS_beta_yX90 = np.zeros(iters)
RMS_beta_yRelEr90 = np.zeros(iters)

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

    RMS_beta_yI90[iter] = Integrals.rms_betatron()[1]
    RMS_beta_yX90[iter] = tw_rad.eq_gemitt_y
    RMS_beta_yRelEr90[iter] = RMS_beta_yI90[iter]/RMS_beta_yX90[iter] - 1


print(f'RMS_beta_yI = {RMS_beta_yI90} [s⁻¹]')
print(f'RMS_beta_yX = {RMS_beta_yX90} [s⁻¹]')
print(f'Relative Error = {RMS_beta_yRelEr90}')


# Plotting
# -------------------------------------------------------------------------------------------------------------------------------
# Plot momentum compaction
# Create the plot
#RelErrs = plt.figure(figsize=(10, 6))
#plt.plot(k0_values, RMS_beta_yRelEr90, label='RMS_beta_yRelEr90')
#plt.plot(k0_values, RMS_beta_yRelEr0, label='RMS_beta_yRelEr0')

#plt.axhline(y=RMS_beta_yRelErbase, color='gray', linestyle='--', linewidth=1, label='No Chicane')

# Set x-axis to logarithmic scale
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
#plt.xlabel('k_wig [m⁻¹]')
#plt.ylabel('RMS_beta_y Relative Errors')
#plt.legend(loc='best')  # Automatically place the legend in an optimal location
#plt.title('Relative Errors in the RMS Betatron Oscillations')

#plt.show()
#RelErrs.savefig("Figures/RMS_beta_yRelErrs.pdf", bbox_inches='tight')

# Plot relative errors
# Create the plot
#Values = plt.figure(figsize=(10, 6))

# Plot the first set of data
#plt.plot(k0_values, RMS_beta_yI0, label='RMS_beta_yI0')
#plt.plot(k0_values, RMS_beta_yX0, label='RMS_beta_yX0')

# Plot the second set of data
#plt.plot(k0_values, RMS_beta_yI90, label='RMS_beta_yI90')
#plt.plot(k0_values, RMS_beta_yX90, label='RMS_beta_yX90')

#plt.axhline(y=RMS_beta_yIbase, color='gray', linestyle='--', linewidth=1, label='Integral, No Chicane')
#plt.axhline(y=RMS_beta_yXbase, color='black', linestyle='--', linewidth=1, label='XSuite, No Chicane')

# Set x-axis to logarithmic scale
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(1e-5, 1e-1)  # Set the limits for the x-axis

# Add labels and legend
#plt.xlabel('k_wig [m⁻¹]')
#plt.ylabel('RMS_beta_y [-]')
#plt.legend(loc='best')  # Automatically place the legend in an optimal location
#plt.title('RMS Betatron Oscillations for Two Angles')

# Show the plot
#plt.show()
#Values.savefig("Figures/RMS_beta_yValues.pdf", bbox_inches='tight')