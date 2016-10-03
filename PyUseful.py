# This is a python module that contains several useful python methods
# R. Sheehan 26 - 9 - 2016


"""
Import Libraries
"""

import os # operating system interface
import sys # system specific parameters and functions
import glob # use this to check whether certain files are in a directory
import subprocess # use this to run executables

import matplotlib.pyplot as plt #  this is needed for plotting solutions


"""
Function Definitions
"""

def poly_eval(coeffs, pos):
    # evaluate the polynomial defined by the set coeffs at position pos

    if coeffs is not None and len(coeffs) > 0:
        

        return poly_val
    else:
        return 0.0

"""
List Manipulation Methods
"""

def list_1D_array(size):
	# return a python list of length size whose elements are all set to zero
    # indexing starts at zero and ends at size-1
	# alternative is to use [None]*size
	# this method uses list comprehensions

    return [0 for i in range(size)] if size > 0 else None

def list_2D_array(rows, columns):
    # return a 2D python list of size rows*columns whose elements are all set to zero
    # indexing starts at zero and ends at rows-1 / columns-1
    # array accessed as Z[row_index, col_index]
	
    if rows > 0 and columns > 0:
        
        ret_lst = list_1D_array(rows)

        for i in range(rows):
            ret_lst[i] = list_1D_array(columns)

        return ret_lst
    else:
        return None

def transpose_multi_col(data):
    # generate the transpose of a multi-column list
    # this will return the transpose of an array in a manner
    # that is convenient for plotting
    # this is arguably more convenient, and pythonic, than using the get_col or get_row method
    # this method assumes that each column has the same length

    return list( map( list, zip(*data) ) )

def get_matrix_dims(matrix):
    # retrieve the dimensions of a 2D array
    # R. Sheehan 18 - 5 - 2016

    if matrix is not None:
        row_size = len( matrix )
        col_size = len( matrix[0] )

        return [row_size, col_size]
    else:
        return None

"""
File IO
"""
def count_lines(thedata, thepath, loud = False):
    # count the number of lines in a file that has been opened and read into memory
    # thedata is the stream that contains the data from an open file
    # how do you know if thedata contains data?
    # assume that you only call count_lines inside another function for reading data
    # thepath is the name of the file containing the data
    # R. Sheehan 26 - 4 - 2014

    nlines=0
    for lines in thedata:
        nlines = nlines + 1

    if loud:
        print "There are %(nlines)d lines in %(path)s"%{"nlines":nlines,"path":thepath}
    
    return nlines

def read_data(thepath, loud = False):
    # read a single column of data from an open file
    # R. Sheehan 26 - 4 - 2014

    if glob.glob(thepath):

        thefile = file(thepath,"r") # open file for reading

        # check that the files are available
        if thefile.closed:
            print "%(path)s could not be opened"%{"path":thepath}
            return None # return an empty array
        else:
            if loud: print "%(path)s is open"%{"path":thepath}

            thedata = thefile.readlines() # read the data from the file

            nlines = count_lines(thedata, thepath) # count the number of data points

            if nlines > 0:

                datapts = list_1D_array(nlines) # create an array of zeros of length nlines

                i=0
                for lines in thedata:
                    datapts[i] = float(lines)
                    i = i + 1

                del thedata # clear thedata stream

                del thefile

                return datapts
            else:
                return None

def read_matrix(thepath, ignore_first_line = False, loud = False):
    # read an array of data from a file
    # if ignore_first_line == True, this means the first line of the file
    # contains text and should not be counted when reading in data
    # R. Sheehan 8 - 8 - 2014

    thefile = file(thepath,"r") # open file for reading

    # check that the files are available
    if thefile.closed:
        print "%(path)s could not be opened"%{"path":thepath}
        return None
    else:
        if loud: print "%(path)s is open"%{"path":thepath}

        thedata = thefile.readlines() # read the data from the file

        nrows = count_lines(thedata, thepath) # count the number of rows

        print "Nrows = ",nrows
        
        if ignore_first_line:
            ncols = len(thedata[1].split(',')) # count the number of columns
            nrows -= 1 # substract one from the number of rows
        else:
            ncols = len(thedata[0].split(',')) # count the number of columns

        if loud: print "rows = ",nrows,"cols = ",ncols

        datapts = list_2D_array(nrows, ncols) # use the native list object instead of numpy
            
        for i in range(0, nrows, 1 ):
            for j in range(0, ncols,1):
                if ignore_first_line:
                    datapts[i][j] = float( thedata[i+1].split(',')[j] )
                else:
                    datapts[i][j] = float( thedata[i].split(',')[j] )
        
    del thefile

    return datapts
