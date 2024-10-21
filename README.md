**Synchrotron Integrals**

This class computes the synchrotron integrals for a given line object in XSuite and calculates associated physical quantities from these integrals.

**Initialization:**
- Pass a line object.
- The class calculates the Twiss-Parameters of the line
- The class extracts the required parameters from the line and the Twiss table

**Public methods:**

*momentum_compaction()*

Calculates the momentum compaction factor.
Returns a (2, ) array of the momentum compaction in the x- and y-direction for the first and second element respectively.

*energy_loss()*

Calculates the energy loss per turn due to synchrotron radiation.
Returns a scalar.

*radiation_damping_s()*

Calculates the exponential damping factors per second in the x-, y- and s-direction.
Returns a (3, ) array of the damping factors in the x-, y- and s-direction for the first, second and third element respectively.

*radiation_damping_turns()*

Calculates the exponential damping factors per turn in the x-, y- and s-direction.
Calls radiation_damping_s() and multiplies with the revolution period.
Returns a (3, ) array of the damping factors in the x-, y- and s-direction for the first, second and third element respectively.

*rms_energies()*

Calculates the RMS of the energy.
Returns a scalar.

*rms_betatron()*

Calculates the RMS betatron amplitude.
Returns a (2, ) array of the RMS betatron amplitudes in the x- and y-direction for the first and second element respectively.

*equilibrium_emittance()*

Calculates the equilibrium emittance in the x- and y-direction.
Calls rms_betatron() and the voltage and frequency of the RF cavities.
Returns a (2, ) array of the equilibrium emittances in the x- and y-direction for the first and second element respectively.
