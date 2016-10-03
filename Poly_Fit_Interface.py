# Python script for running the Slab-Waveguide Solver executable
# and plotting the results of the solutions generated
# R. Sheehan 7 - 9 - 2016

"""
Import Libraries
"""

import os # operating system interface
import sys # system specific parameters and functions
import glob # use this to check whether certain files are in a directory
import subprocess # use this to run executables

import matplotlib.pyplot as plt #  this is needed for plotting solutions

import PyUseful

"""
Function Definitions
"""

def run_poly_fit(x_name, y_name, degree, storage):
    # run the executable for the least squares polynomial fit code

    # where is the executable located?
    exe_dir = "C:/Users/Robert/Programming/C++/Demo_Projects/Poly_Fit/Release/"

    # what is the name of the executable?
    # prog_name needs a space after it to distinguish it from arg_vals
    prog_name = "Poly_Fit.exe" + " "

    # convert arguments to a string
    # need a space between arguments and "\\" added to storage
    arg_vals = x_name + " " + y_name + " " + "{:d}".format(degree) + " " + storage + "\\"

    print arg_vals

    # args is the value that is passed to subprocess
    args = exe_dir + prog_name + arg_vals

    # subprocess.call to run the program without printing to the python shell
    # shell=False is to be used as standard practice unless you
    # really know what you're doing!
    # output = subprocess.call(args, stdin=None, stdout=None, stderr=None, shell=False)

    # use subprocess.check if you want to print the output from the
    # program to the python shell
    # shell=False is to be used as standard practice unless you
    # really know what you're doing!
    print subprocess.check_output(args, shell=False)

def plot_input_data(x_data, y_data, fit_coeffs, loud = False):
    # plot the TE and TM dispersion equations

    if x_data is not None and y_data is not None and fit_coeffs is not None:

        col_label = 'r*'
        lin_label = 'g-'
        plt_label = "Data"
        fig_label = "Data.png"

        fit_x = [x_data[0], x_data[-1]]
        fit_y = [fit_coeffs[0] + x_data[0]*fit_coeffs[1], fit_coeffs[0] + x_data[-1]*fit_coeffs[1]]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.plot(x_data, y_data, col_label)
        ax.plot(fit_x, fit_y, lin_label)
        

        plt.xlabel('x data', fontsize = 16)
        plt.ylabel('y data', fontsize = 16)
        plt.title(plt_label)
        spacer = 1.0
        plt.axis( [x_data[0] - spacer, x_data[-1] + spacer, min(y_data) - spacer, max(y_data) + spacer] )

        plt.savefig(fig_label)
        if loud: plt.show()
        plt.clf()
        plt.cla()
        plt.close()
        
"""
Call and run script from inside main
"""

def main():
    pass

if __name__ == '__main__':
    main()
    
    # what arguments are needed by prog_name?
    # Slab_WG_Slv needs
    x_name = "X_data_Ex_4.txt"
    y_name = "Y_data_Ex_4.txt"
    degree = 1
    
    storage = os.getcwd() # store the solutions in the current directory

    run_poly_fit(x_name, y_name, degree, storage)

    x_data = PyUseful.read_data(x_name)
    y_data = PyUseful.read_data(y_name)

    fit_coeffs = PyUseful.read_data("Fit_Coefficients.txt")

    plot_input_data(x_data, y_data, fit_coeffs, loud = True)


