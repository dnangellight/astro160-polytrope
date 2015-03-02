#!/usr/bin/env python2.7

################################################################################
## This script is for finding the density profile of a polytrope.
## Copyright (C) 2013  Rachel Domagalski: idomagalski@berkeley.edu
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os    as _os
import sys   as _sys
import numpy as _np

class Polytrope:
    """
    Polytrope class.
    """
    def __init__(self, nsteps, index, rho_0, pressure_0, radius):
        """
        Initialize the polytrope and assign values to the initial
        conditions, as well as setting the polytropic index.
        """
        self.nsteps     = int(nsteps)
        self.index      = float(index)
        self.rho_0      = float(rho_0)
        self.pressure_0 = float(pressure_0)
        self.radius     = float(radius)

        alpha, beta = self.get_ab(self.index, self.rho_0, self.pressure_0)
        self.alpha  = alpha
        self.beta   = beta

    def get_ab(self, index, rho_0, pressure_0):
        """
        Compute the alpha and beta constants used in the solutions.
        """
        # Physical constants
        G = 6.67259e-8 # Require CGS units

        gamma = (self.index + 1.0) / self.index
        K = self.pressure_0 / self.rho_0 ** gamma
        beta_sq = K * gamma / (4 * _np.pi * G)
        beta = _np.sqrt(beta_sq)

        # Ratio of my beta constant to the book's alpha constant
        index_power = (1.0 - self.index) / self.index
        ratio_ab = self.index * self.rho_0 ** index_power
        ratio_ab = _np.sqrt(ratio_ab)
        alpha = ratio_ab * beta

        return (alpha, beta)

    def update_index(self, index):
        """
        Update the polytropic index.
        """
        self.index  = float(index)
        alpha, beta = self.get_ab(self.index, self.rho_0, self.pressure_0)
        self.alpha  = alpha
        self.beta   = beta

        # Remove old solutions
        if hasattr(self, 'radpts'):
            del self.radpts
        if hasattr(self, 'density'):
            del self.density
        if hasattr(self, 'gradient'):
            del self.gradient

    def save_data(self, radpts, density, gradient):
        """
        Save a solution to a numpy array.
        """
        self.radpts   = _np.array(radpts)
        self.density  = _np.array(density)
        self.gradient = _np.array(gradient)

    def derivative(self, radius, rho, u):
        """
        Second derivative of density defined in terms of the helper
        function u.
        """
        beta_sq     = self.beta ** 2
        mult_index  = (self.index - 1.0) / self.index
        exp_index   = (2*self.index - 1.0) / self.index
        derivative  = mult_index * u**2 / rho
        derivative -= 2 * u / radius
        derivative -= rho ** exp_index / beta_sq
        return derivative

    def rk4_terms(self, radius, rho, u, step_size):
        """
        Calculates the terms for an iteration in RK4.
        """
        coefs = [1, 2, 2, 1]
        exp_index = (2*self.index - 1.0) / self.index

        # Set the first term
        k1 = [u]
        k2 = [self.derivative(radius, rho, u)]

        # Set the second and third terms.
        for i in range(1,3):
            rad_step = radius + step_size / 2.0
            rho_step = rho    + step_size * k1[-1] / 2.0
            u_step   = u      + step_size * k2[-1] / 2.0

            # Python will object to fractional exponents of negative numbers
            if rho_step < 0 and exp_index != int(exp_index):
                return (None, None)

            k1.append(u_step)
            k2.append(self.derivative(rad_step, rho_step, u_step))

        # Calculate the last term
        rad_step = radius + step_size
        rho_step = rho    + step_size * k1[-1]
        u_step   = u      + step_size * k2[-1]

        if rho_step < 0 and exp_index != int(exp_index):
            return (None, None)

        k1.append(u_step)
        k2.append(self.derivative(rad_step, rho_step, u_step))

        # Scale each term by the proper amount so I can use them in a sum.
        k1 = [coefs[i] * k1[i] for i in range(4)]
        k2 = [coefs[i] * k2[i] for i in range(4)]
        return (k1, k2)

    def rk4_outward(self):
        """
        This function runs a runge-kutta 4 solver on the system outward
        from the center of the star.
        """
        # Initialize the solver
        self.radpts   = [0.007 * self.radius]
        self.density  = [self.rho_0]
        self.gradient = [0.0] # The center derivative
        step_size     = (self.radius - self.radpts[0]) / self.nsteps

        # Run the solver and save the results.
        for i in range(self.nsteps):
            last_rad = self.radpts[-1]
            last_rho = self.density[-1]
            last_u   = self.gradient[-1]

            # Get the coefficients and check that they are valid. If the
            # coefficients cannot be computed, save the current data and stop
            # the solver.
            k1, k2   = self.rk4_terms(last_rad, last_rho, last_u, step_size)
            if k1 == None or k2 == None:
                self.save_data(self.radpts, self.density, self.gradient)
                break

            # Save the step in the iteration
            next_rad = last_rad + step_size
            next_rho = last_rho + step_size * sum(k1) / 6.0
            next_u   = last_u   + step_size * sum(k2) / 6.0

            # Negative density makes no physical sense.
            if next_rho <= 0:
                self.save_data(self.radpts, self.density, self.gradient)
                break
            else:
                self.radpts.append(next_rad)
                self.density.append(next_rho)
                self.gradient.append(next_u)

        # Save the finished data as numpy arrays.
        self.save_data(self.radpts, self.density, self.gradient)
