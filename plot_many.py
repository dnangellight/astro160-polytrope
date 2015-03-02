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
    print 'Usage: python2 plot_many.py -i max_index -n nsteps [-x] [-y]'
    print 'Flags:'
    print '    -h: Print this message and exit.'
    print '    -i: Maximum polytopic index.'
    print '    -n: Number of steps in the radius.'
    print '    -x: Use a log scale on the x-axis.'
    print '    -y: Use a log scale on the y-ayis.'
    sys.exit(code)

## RUN EVERYTHING ##
if __name__ == '__main__':
    # Parse arguments and get options
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:n:xy')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Check that the options inputted to the script are valid.
    index     = None
    nsteps    = None
    log_x     = False
    log_y     = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-i':
            try:
                index = int(arg)
                # No valid solutions to Poisson's equation exist for n < 1.
                if index < 1:
                    print 'ERROR: Maximum index must be at least 1.'
                    sys.exit(1)
            except ValueError:
                print 'ERROR: Maximum index must be an integer.'
                sys.exit(1)
        # Enforce that the number of steps is an integer.
        elif opt == '-n':
            try:
                nsteps = int(arg)
            except ValueError:
                print 'ERROR: Number of steps must be an integer.'
                sys.exit(1)
        # Log scale
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
    indices = range(1, index+1) + [1.5]
    indices.sort()
    for index in indices:
        star.update_index(index)
        star.rk4_outward()
        radpts  = star.radpts / star.radius
        density = star.density / star.rho_0
        plot_label = r'$n = ' + str(index) + '$'
        plt.plot(radpts, density, label=plot_label)

    axes = plt.axes()
    plt.xlabel(r'$R / R_*$')
    plt.ylabel(r'$\rho / \rho_c$')
    plt.title('Solutions to Poisson\'s equation')

    # Set the scales and configure the axes
    if log_x:
        plt.xscale('log')
    if log_y:
        plt.yscale('log')

    # Configure the axes
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.xaxis.set_ticks_position('bottom')
    axes.yaxis.set_ticks_position('left')
    plt.axis([min(radpts), max(radpts), min(density), 1])

    plt.legend()

    plt.show()
