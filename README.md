Densities of polytropic stars
=============================

This code is a solution to a problem in my stellar astrophysics class where I am
computing numerical solutions to Poisson's eqation for a polytopic star. The
scripts that I am using to create plots are using data about the sun to generate
the solutions.

Running the scripts:
-----------------------------------

Running the Python scripts is fairly simple and straightforward. I am using
ipython for these, as it is good with matplotlib.

To create a plot of a solution for Poisson's equation with one polytropic index,
run the following:  
$ ipython -- plot\_one.py -i index -n nsteps  
It's also possible add a curve generated from a data file to the plot as well as
setting the x and y axes to a log scale. To see all of the available options
that the script allows, enter the following:  
$ python plot\_one.py -h (or with ipython: ipython -- plot\_one.py -h)  

To create a plot with many solutions on a single plot, enter the following
command:  
$ ipython -- plot\_many.py -i maxindex -n nsteps  
To see the all of the available options, pass the -h flag to the plot\_many.py
script. 

The plot\_one.py and plot\_many.py scripts use solar data to set the initial
conditions for the numerical solver. Actual solar data, which I took from my
textbook, can be found in solar\_interior.txt, which consists of radii in units
of solar radius and density in cgs units.

Have fun!
