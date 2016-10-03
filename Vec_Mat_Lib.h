#ifndef VEC_MAT_LIB_H
#define VEC_MAT_LIB_H

// This library contains functions that dynamically create vectors and matrices of type double, of arbitrary size
// Indexing of arrays is one-based, i.e. indices start at 1 and end at N
// This is different from the C/C++ default of zero-based array indexing, i.e. indices start at zero and end at N-1
// R. Sheehan 31 - 1 - 2013

namespace vec_mat_funcs{

	double *vector(int size); // this function will return a vector that can hold size elements of type double

	double *zero_vector(int size); // dynamically allocate an array filled with zeroes

	double *matrix(int rows, int columns); // this function will return a matrix that can hold row*cols elements of type double

	double *zero_matrix(int rows, int columns); // dynamically allocate an array filled with zeroes

	double *identity_matrix(int size); // dynamically allocate an array to return the order N identity matrix

	void print_vector(double *vec, int size); // print a vector to the screen

	void print_matrix(double *mat, int rows, int columns); // print a matrix to the screen

	double *read_vector_from_file(std::string filename, int &n_data); 

	void read_vector_from_file(std::string filename, int n_data, double *arr); 
}

#endif