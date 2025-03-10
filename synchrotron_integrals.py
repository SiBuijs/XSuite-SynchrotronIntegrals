# RadiationIntegral Class
# Works for both thick and thin elements
# By Simon Fanica Buijsman
# Date: 2024-11-28
# simon.fanica.buijsman@cern.ch


# TODO: Check if the curvature is calculated correctly for off-momentum particles.
# The current challenge is to see if this simulation is valid for particles that do not follow the design orbit.
# So far, I've only tested the simulation for the design trajectory. In this case, the curvature is trivial.


#### ---------------------------------------------------------------------------------------------------------------------------------
#### ---------------------------------------------------------------------------------------------------------------------------------
import numpy as np

# Physical constants required to calculate the physical quantities associated with the radiation integrals.
from scipy.constants import c as clight, hbar, electron_volt


#dict_keys(['name', 's', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'ptau', 'W_matrix', 'kin_px', 'kin_py', 'kin_ps', 'kin_xprime', 'kin_yprime',
#           'betx', 'bety', 'alfx', 'alfy', 'gamx', 'gamy', 'dx', 'dpx', 'dy', 'dpy', 'dx_zeta', 'dpx_zeta', 'dy_zeta', 'dpy_zeta', 'betx1', 
#           'bety1', 'betx2', 'bety2', 'mux', 'muy', 'muzeta', 'nux', 'nuy', 'nuzeta', 'dzeta', 'only_markers', 'particle_on_co', 'circumference', 
#           'orientation', 'R_matrix', 'steps_r_matrix', 'R_matrix_ebe', 'slip_factor', 'momentum_compaction_factor', 'bets0', 'T_rev0', 'gamma0', 
#           'beta0', 'p0c', 'qx', 'qy', 'qs', 'c_minus', 'c_r1_avg', 'c_r2_avg', 'eigenvalues', 'rotation_matrix', 'dmux', 'dmuy', 'bx_chrom', 'by_chrom', 
#           'ax_chrom', 'ay_chrom', 'wx_chrom', 'wy_chrom', 'ddx', 'ddpx', 'ddy', 'ddpy', 'dqx', 'dqy', 'ddqx', 'ddqy', 'eneloss_turn', 'damping_constants_turns', 
#           'damping_constants_s', 'partition_numbers', 'eq_gemitt_x', 'eq_gemitt_y', 'eq_gemitt_zeta', 'eq_nemitt_x', 'eq_nemitt_y', 'eq_nemitt_zeta', 
#           'dl_radiation', 'n_dot_delta_kick_sq_ave', 'values_at', 'k0l', 'k1l', 'k2l', 'k3l', 'k4l', 'k5l', 'k0sl', 'k1sl', 'k2sl', 'k3sl', 'k4sl', 'k5sl', 
#           'angle_rad', 'rot_s_rad', 'hkick', 'vkick', 'ks', 'length', 'element_type', 'isthick', 'parent_name', 'method', 'radiation_method', 'reference_frame', 
#           'line_config', '_action'])


class SynchrotronIntegral:

    def __init__(self, line):

        # The class only needs a line object to be instantiated.
        # The line object contains all the information about the elements, including the Twiss.
        self.line = line
        self.tw = line.twiss(method='4d', strengths=True)#, eneloss_and_damping=True)


        # Properties of the particle -------------------------------------------------------------------------------------------------
        # The values for gamma0 and energy0 are indexed with [0] because they are arrays of shape (1, ).
        # To prevent getting arrays as the output of functions where we just wamt to get scalars, the indexing is performed here.
        self.mass0 = self.line.particle_ref.mass0
        self.gamma0 = self.line.particle_ref.gamma0[0]
        self.energy0 = self.line.particle_ref.energy0[0]
        self.r_e = self.line.particle_ref.get_classical_particle_radius0()

        # Properties of the beam line ------------------------------------------------------------------------------------------------
        self.length = self.tw['length']                 # The length of each element
        self.circum = self.tw['circumference']          # The circumference of the ring
        self.T_0 = self.tw['T_rev0']                    # The revolution period of the beam.

        # Particle coordinates
        self.x = self.tw['x']                           # The x-coordinate of the particle
        self.y = self.tw['y']                           # The y-coordinate of the particle
        self.s = self.tw['s']                           # The s-coordinate of the particle

        self.px = self.tw['px']                         # The x-momentum of the particle
        self.py = self.tw['py']                         # The y-momentum of the particle

        # Curvatures for each element
        self.k = self._get_curvature_()                 # The curvature of each element
        self.kmag = np.linalg.norm(self.k, axis=1)      # The magnitude of the curvature of each element
        #self.kx = self._get_curvature_()[0]             # The integrated normal curvature of each element
        #self.ky = self._get_curvature_()[1]             # The integrated skew curvature of each element

        # Dispersions ----------------------------------------------------------------------------------------------------------------
        # x-direction
        self.dx   = self.tw['dx']                       # Dispersion in the x-direction
                                                        # Derivative of the dispersion in the x-direction, w.r.t. s
        self.dxprime  = self.tw['dpx'] * (1 - self.tw['delta']) - self.tw['px']

        # y-direction
        self.dy   = self.tw['dy']                       # Dispersion in the y-direction
                                                        # Derivative of the dispersion in the y-direction, w.r.t. s
        self.dyprime  = self.tw['dpy'] * (1 - self.tw['delta']) - self.tw['py']


    # Private methods ----------------------------------------------------------------------------------------------------------------
    # Below are functions that calculate physical quantities that are later used in the integrals.


    # This function calculates the horizontal and vertical acceleration of the particle.
    def _get_derivative_(self, yaxis):
            derivative_values = np.zeros(shape=(len(self.length)))
            dx = np.zeros(shape=(len(self.length)))
            dx[:-1] = self.length[1:] + self.length[0:-1]
            dx[-1] = self.length[0] + self.length[-1]

            derivative_values[1:-1] = yaxis[2:] - yaxis[0:-2]
            derivative_values[0] = yaxis[1] - yaxis[-1]
            derivative_values[-1] = yaxis[0] - yaxis[-2]


            #derivative_values = np.gradient(yaxis)
            mask = dx != 0
            derivative_values[mask] = derivative_values[mask] / dx[mask]
            
            return derivative_values


    # Calculate the curvature of the design orbit.
    def _orbit_curvature_(self):
        kappa_xy = np.zeros(shape=(2, len(self.length)))
        
        # TODO: I think the self.tw['k0l'] should be replaced with self.tw['angle_rad'].
        angle_rad = self.tw['angle_rad']
        rot_s_rad = self.tw['rot_s_rad']
        mask = self.length != 0

        kappa_xy[0, :][mask] = angle_rad[mask] * np.cos(rot_s_rad[mask]) / self.length[mask]
        kappa_xy[1, :][mask] = angle_rad[mask] * np.sin(rot_s_rad[mask]) / self.length[mask]

        return kappa_xy


   # Calculate the curvature of the design orbit.
    def _get_curvature_(self):
        kxy = np.zeros(shape=(len(self.length), 3))
        
        #self._get_derivative_(self.x, self.s)
        #self._get_derivative_(self.y, self.s)
        
        pxprime = self._get_derivative_(self.px)
        pyprime = self._get_derivative_(self.py)

        kappa_0xy = self._orbit_curvature_()
        h = 1 + kappa_0xy[0, :] * self.x + kappa_0xy[1, :] * self.y
        hprime = kappa_0xy[0, :] * self.px + kappa_0xy[1, :] * self.py

        #kxy[:, 0] = - (h * (pxprime - h * kappa_0xy[0, :]) - 2 * self.px * (kappa_0xy[0, :] * self.px + kappa_0xy[1, :] * self.py)) / (self.px**2 + h**2)**(3/2)
        #kxy[:, 1] =   (h * (pyprime - h * kappa_0xy[1, :]) - 2 * self.py * (kappa_0xy[0, :] * self.px + kappa_0xy[1, :] * self.py)) / (self.py**2 + h**2)**(3/2)

        rprime = np.stack([self.px, self.py, h], axis=1)  # Shape (N, 3)
        rdoubleprime = np.stack([pxprime - h*kappa_0xy[0, :], pyprime - h * kappa_0xy[1, :], 2*hprime], axis=1)  # Shape (N, 3)

        #cross_products = np.cross(rprime, rdoubleprime)

        magnitudesvelocities = np.linalg.norm(rprime, axis=1).reshape(-1, 1)

        mask = (magnitudesvelocities != 0).reshape(-1, )

        kxy[mask] = rdoubleprime[mask] / magnitudesvelocities[mask]**2

        return kxy


    # For the field index
    # NOTE: Finally works. The radius in the field index is the radius of curvature of the design orbit, not of the particle orbit.
    def _get_fieldindex_(self):
        # The integrated normal and skew quadrupole strengths of each element
        quadkn = self.tw['k1l']            # The integrated normal quadrupole strength of each element
        #quadks = self.tw['k1sl']          # The integrated skew quadrupole strength of each element

        reforbitkappax = self._orbit_curvature_()[0]
        reforbitkappay = self._orbit_curvature_()[1]


        fieldindex = np.zeros(shape=(2, len(self.length)))
        mask = self.length * self.kmag * reforbitkappax != 0
        fieldindex[0, :][mask] = - quadkn[mask] / (self.length[mask] * self.kmag[mask] * reforbitkappax[mask])
        mask = self.length * self.kmag * reforbitkappay != 0
        fieldindex[1, :][mask] = - quadkn[mask] / (self.length[mask] * self.kmag[mask] * reforbitkappay[mask])
            
        return fieldindex


    # Calculates the H function
    def _H_function_(self):
        betx = self.tw['betx']             # Twiss beta function x
        alfx = self.tw['alfx']             # Twiss alpha x
        gamx = self.tw['gamx']             # Twiss gamma x
        bety = self.tw['bety']             # Twiss beta function y
        alfy = self.tw['alfy']             # Twiss alpha y
        gamy = self.tw['gamy']             # Twiss gamma y

        H_values_xy = np.zeros(shape=(2, len(self.length)))
            
        H_values_xy[0, :] = gamx * self.dx**2 + 2*alfx * self.dx * self.dxprime + betx * self.dxprime**2
        H_values_xy[1, :] = gamy * self.dy**2 + 2*alfy * self.dy * self.dyprime + bety * self.dyprime**2

        return H_values_xy


    # The integrands -----------------------------------------------------------------------------------------------------------------
    # Below, we calculate the integrands. The integrals are arrays of length equal to the number of elements.
    # To get the integrals, we need to sum over al the elements in te integrals.
    # This is done in the methods that calculate the physical quantities.


    # Integrand 1
    # Units: [m]
    # Returns an array of shape (2, number of elements)
    # The two rows correspond to the x- and y-values respectively.
    # NOTE: Don't know why the minus sign is there.
    def _integrand1_(self):
        I1xy_values = np.zeros(shape=len(self.length))

        disp = -np.stack([self.dx, self.dy, np.zeros(shape=len(self.dx))], axis=1)

        dot_product = np.sum(disp * self.k, axis=1)

        I1xy_values = dot_product * self.length

        return I1xy_values
        
    
    # Integrand 2
    # Units: [m⁻¹]
    def _integrand2_(self):
        I2_values = np.zeros(shape=len(self.length))

        I2_values = self.length * self.kmag**2
        
        return I2_values


    # Integrand 3
    # Units: [m⁻²]
    def _integrand3_(self):
        I3xy_values = np.zeros(shape=len(self.length))

        I3xy_values = self.length * self.kmag**3

        return I3xy_values
            
    
    # Integrand 4
    # Units: [m⁻¹]
    # Returns an array of shape (2, number of elements)
    # The two rows correspond to the x- and y-values respectively.
    # NOTE: I am going to mess up the field index here, so if it doesn't work, change that back first!
    def _integrand4_(self):
        I4xy_values = np.zeros(shape=(2, len(self.length)))


        fieldindex = self._get_fieldindex_()
        disp = -np.stack([self.dx, self.dy, np.zeros(shape=len(self.dx))], axis=1)
        dot_product = np.sum(disp * self.k, axis=1)

        I1xy_values = dot_product * self.length

        I4xy_values[0, :] = self.length * self.kmag**3 * self.dx * (1 - 2 * fieldindex[0])
        I4xy_values[1, :] = self.length * self.kmag**3 * self.dy * (1 - 2 * fieldindex[1])

        return I4xy_values
        
    
    # Integrand 5
    # Units: [m⁻¹]
    def _integrand5_(self):
        I5xy_values = np.zeros(shape=(2, len(self.length)))
        H = self._H_function_()
            
        I5xy_values[0, :] = self.length * self.kmag**3 * H[0, :]
        I5xy_values[1, :] = self.length * self.kmag**3 * H[1, :]
            
        return I5xy_values


    # The physical quantities --------------------------------------------------------------------------------------------------------
    # Below are the public methods that calculate the physical quantities using the integrands.

    # Momentum Compaction Factor
    # This function gets the (2, number of elements) array of the _Integrand1_
    # and returns a (2, 1) array with the x- and y-values for the momentum compaction factor respectively.
    # Call momentum_compaction()[0] for the x-value
    # and momentum_compaction()[1] for the y-value.
    def momentum_compaction(self):
        I1 = np.sum(self._integrand1_())

        return I1 / self.circum


    # Energy Loss
    # Calculates the energy loss using _Integrand2_().
    # Integrand 2 only depends on the length and bending radii of the elements.
    # This means that for the total energy loss, we can simply sum the x- and y-integrals.
    # This amounts to summing all elements of _Integrand2_() together.
    def energy_loss(self):
        I2 = np.sum(self._integrand2_())

        return 2/3 * self.r_e * self.gamma0**3 * self.energy0 * I2


    # Radiation Damping per Second
    # Calculates the radiation damping using _Integrand2_() and _Integrand4_().
    # I2xy = I2x + I2y, which is just the sum over all elements of _Integrand2_().
    # I4x and I4y are the x- and y-values of the sum over all elements of _Integrand4_().
    # I4xy = I4x + I4y.Equation
    #
    # The radiation damping coefficients are calculated from symmetry considerations.
    # For example, if only horizontal bending etc. is present, alpha_x is proportional to I_2 - I_4.
    # Also, alpha_y is proportional to just I_2.
    # So far, I_2 and I_4 have only been in the x-direction.
    # If we want to include the effects of bending in the y-direction, we can switch the roles of alpha_x and alpha_y.
    # This would result in alpha_x being proportional to I_2 and alpha_y being proportional to I_2 - I_4, which are now calculated in the y-direction.
    # The total effect is assumed to be the sum of the effects in the x- and y-directions.
    # This results in adding an I_2y proportionality to alpha_x and an I_2x proportionality to alpha_y.
    def radiation_damping_s(self):
        rad_damp_coefs = np.zeros(shape=(3,))
        
        I2xy = np.sum(self._integrand2_())
        I4x  = np.sum(self._integrand4_()[0])
        I4y  = np.sum(self._integrand4_()[1])
        I4xy = np.sum(self._integrand4_())
            
        rad_damp_coefs[0] = self.r_e/3 * self.gamma0**3 * clight/self.circum * (I2xy - I4x)
        rad_damp_coefs[1] = self.r_e/3 * self.gamma0**3 * clight/self.circum * (I2xy - I4y)
        rad_damp_coefs[2] = self.r_e/3 * self.gamma0**3 * clight/self.circum * (2 * I2xy + I4xy)

        return rad_damp_coefs
        

    # Radiation Damping per Turn
    # The radiation damping per turn is calculated by multiplying the radiation damping per second with the revolution period.
    def radiation_damping_turns(self):
        return self.radiation_damping_s() * self.T_0


    # RMS Energies
    # It is assumed that the RMS of the energies depends both on the x- and y-integrals.
    def rms_energies(self):
        I2 = np.sum(self._integrand2_())
        I3 = np.sum(self._integrand3_())
        I4 = np.sum(self._integrand4_())

        if 2 * I2 + I4 != 0:
            return 55/(32 * 3**(1/2)) * hbar / electron_volt * clight / self.mass0 * self.gamma0**2 * I3 / (2 * I2 + I4)
        
        else:
            raise ValueError("The denominator of the RMS energy, 2I2 - I4 = 0. The RMS energy is not defined.")
        

    # Radial Quantum Emission
    # Since the denominator depends on the damping contants, it is expected that it contains the same integrals as their respective damping constants.
    # That is, the x-component has I2xy - I4x in the denominator and the y-component has I2xy - I4y.
    def rms_betatron(self):
        rms_betatronxy = np.zeros(shape=(2,))
        
        I2xy = np.sum(self._integrand2_())
        I4x  = np.sum(self._integrand4_(), axis=1)[0]
        I4y  = np.sum(self._integrand4_(), axis=1)[1]
        I5x  = np.sum(self._integrand5_(), axis=1)[0]
        I5y  = np.sum(self._integrand5_(), axis=1)[1]

        if I2xy - I4x != 0:
            rms_betatronxy[0] = 55/(32 * 3**(1/2)) * hbar / electron_volt * clight / self.mass0 * self.gamma0**2 * I5x / (I2xy - I4x)

        if I2xy - I4y != 0:
            rms_betatronxy[1] = 55/(32 * 3**(1/2)) * hbar / electron_volt * clight / self.mass0 * self.gamma0**2 * I5y / (I2xy - I4y)

        return rms_betatronxy
    

    # Equilibrium emittance
    # This quantity is derived from the RMS Energy.
    # See Hoffmann Section 14.3.2 for more information about the equation used here.
    # TODO: Find out how we can correctly implement the effect of phi_s.
        # Above transition, we must have phi_s = 180°, but this gives a negative cosine and thus an imaginary number.
        # I put it to zero for now, but there is something to be done here.
    # TODO: Is creating this large table not a bit too much for such a small function?
    def equilibrium_emittance(self):
        tab = self.line.get_table(attr = True)
        V_cav = tab.rows[tab.element_type == 'Cavity']['voltage']
        f_cav = tab.rows[tab.element_type == 'Cavity']['frequency']
        lag_cav = 0#tab.rows[tab.element_type == 'Cavity']['lag']

        # self.line.attr['voltage'][mask]
        # self.line.attr['frequency'][mask]

        V_cav_tot = np.sum(V_cav)                           # The voltage of all cavities. See the TODO.
        f_cav = np.average(f_cav[V_cav != 0])               # The frequency of the cavity.
        p_0 = self.line.particle_ref.p0c[0] / clight        # The [0] index is to convert this to a scalar instead of a (1, ) array.
        harmonic = f_cav * self.T_0                         # The harmonic number is the ratio of the cavity frequency to the revolution frequency.
        phi_s = lag_cav                                     # The synchronous phase is assumed to be 0. See the TODO.
        alpha_cI = self.momentum_compaction()               # The momentum compaction factor.

        rms_E_conv_factor = self.T_0 * self.energy0 / p_0 * (alpha_cI * self.energy0 / (2*np.pi * harmonic * V_cav_tot * np.cos(phi_s)))**(1/2)

        return rms_E_conv_factor * self.rms_energies()