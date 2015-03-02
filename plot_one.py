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

import os
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt

from Polytrope import *

def usage(code):
    print 'Usage:',
    print 'python2 plot_one.py -i index -n nsteps [-d data.txt] [-x] [-y]'
    print 'Flags:'
    print '    -h: Print this message and exit.'
    print '    -d: Data file for comparing the model.'
    print '    -i: Polytropic index.'
    print '    -n: Number of steps in the radius.'
    print '    -x: Use a log scale on the x-axis.'
    print '    -y: Use a log scale on the y-ayis.'
    sys.exit(code)

def exact1(radius, rho_0, alpha):
    """
    Computes the exact solution for a polytrope of n = 1. This function
    should also work on numpy arrays.
    """
    xi = radius / alpha
    theta = np.sin(xi) / xi
    return rho_0 * theta

def exact5(radius, rho_0, alpha):
    """
    Computes the exact solution for a polytrope of n = 5. This function
    should also work on numpy arrays.
    """
    xi = radius / alpha
    theta = 1 / np.sqrt(1.0 + xi**2 / 3.0)
    return rho_0 * theta ** 5

## RUN EVERYTHING ##
if __name__ == '__main__':
    # Parse arguments and get options
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:hi:n:xy')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Check that the options inputted to the script are valid.
    data_file = None
    index     = None
    nsteps    = None
    log_x     = False
    log_y     = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = arg
        # Enforce that the polytropic index is a real number.
        elif opt == '-i':
            try:
                index = float(arg)
            except ValueError:
                print 'ERROR: Polytopic index must be a number.'
                sys.exit(1)
        # Enforce that the number of steps is an integer.
        elif opt == '-n':
            try:
                nsteps = int(arg)
            except ValueError:
                print 'ERROR: Number of steps must be an integer.'
                sys.exit(1)
        elif opt == '-x':
            log_x = True
        elif opt == '-y':
            log_y = True

    # Check for mandatory arguments.
    if index == None:
        print 'Error: No polytropic index supplied!'
        sys.exit(1)
    if nsteps == None:
        print 'Error: No step number supplied!'
        sys.exit(1)

    # Initialize the polytrope.
    rho_0      = 150.0
    pressure_0 = 10**17.369
    radius     = 6.95508e10
    star = Polytrope(nsteps, index, rho_0, pressure_0, radius)

    # Run the solver and plot the results.
    star.rk4_outward()
    radpts  = star.radpts / star.radius
    density = star.density / star.rho_0

    # Plot the results.
    axes = plt.axes()
    plt.plot(radpts, density, label='Numerical solution')
    plt.xlabel(r'$R / R_*$')
    plt.ylabel(r'$\rho / \rho_c$')
    if star.index == int(star.index):
        plt.title(r'$n = ' + str(int(star.index)) + '$')
    else:
        plt.title(r'$n = ' + str(star.index) + '$')

    # Set the scales and configure the axes
    if log_x:
        plt.xscale('log')
    if log_y:
        plt.yscale('log')

    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.xaxis.set_ticks_position('bottom')
    axes.yaxis.set_ticks_position('left')
    plt.axis([min(radpts), max(radpts), min(density), 1])

    # If the program has been supplied with a data file, plot it too.
    if data_file != None:
        solar_data = np.genfromtxt(data_file)
        solar_radius  = np.array([data[0] for data in solar_data])
        solar_density = np.array([data[1] for data in solar_data])
        solar_density /= solar_density[0]
        plt.plot(solar_radius, solar_density, label='Solar data')
        if log_x:
            plt.legend(loc=(0.02, 0.84))
        else:
            plt.legend()
        if log_y:
            plt.axis([min(radpts), max(radpts), min(solar_density), 1])

    # For two polytropic indexes, there are exact solutions
    has_exact = False
    if index == 1.0:
        exact = exact1(star.radpts, star.rho_0, star.alpha)
        has_exact = True
    elif index == 5.0:
        exact = exact5(star.radpts, star.rho_0, star.alpha)
        has_exact = True

    # Make plots for exact solutions.
    if has_exact and data_file == None:
        # Add the exact solution to the plot.
        plt.plot(radpts, exact / star.rho_0, label='Exact solution')
        plt.legend()

        # Create error figure
        plt.figure()
        axes = plt.axes()
        plt.plot(radpts, (star.density - exact) / exact)
        plt.xlabel(r'$R / R_*$')
        plt.ylabel('Relative difference')
        plt.title(r'$n = ' + str(int(star.index)) + '$')

        # Set the axes
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.xaxis.set_ticks_position('bottom')
        axes.yaxis.set_ticks_position('left')

    # If the solution must stop early, report where it stops.
    if len(radpts) != nsteps + 1:
        print 'Solution stops at R/R_* =', str(radpts[-1]) + '.'

    plt.show()
