#ifndef ATTACH_H
#include "Attach.h"
#endif

double *vec_mat_funcs::vector(int size) 
{
	// this function will return a vector that can hold size elements of type double
	// the function dynamically creates an array
	// array indexing starts at i = 0 and ends at i = N - 1
	// R. Sheehan 31 - 1 - 2013

	try{
		
		double *vec_ptr = new (double [size]); 

		return vec_ptr; 

		// bad_alloc thrown by default if there is a problem with the dimensions, 
		// or not enough memory is available

	}
	catch(std::bad_alloc &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double *vec_mat_funcs::zero_vector(int size)
{
	try{

		if(size > 0){

			double *arr = new (double [size]); 

			// assign all elements in arr to zero using fill_n from <algorithm>
			// http://www.cplusplus.com/reference/algorithm/fill_n/
			// I would have liked to have used this but the compiler is giving me a really annoying error message
			// that I can't get rid of, so I'll write it myself. 

			if(arr != nullptr){
			
				int i=0; 
				while(i<size){
					*(arr+i) = 0.0; 
					i++; 
				}

			}
	
			return arr;
		
		}
		else{
			std::string reason = "Error: double *vec_mat_funcs::zero_vector(int size)\n"; 
			reason += "size = " + template_funcs::toString(size) + " is not positive\n"; 
			throw std::invalid_argument(reason); 
			return nullptr; 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double *vec_mat_funcs::matrix(int rows, int columns)
{	
	// this function will return a matrix that can hold row*cols elements of type double
	// use dynamic memory allocation to create a matrix
	// row indexing starts at i = 0 and ends at i = rows-1
	// column indexing starts at j = 0 and ends at j = columns-1
	// R. Sheehan 31 - 1 - 2013

	try{
		
		double *mat_ptr = new (double [rows*columns]); 

		return mat_ptr; 

		// bad_alloc thrown by default if there is a problem with the dimensions, 
		// or not enough memory is available

	}
	catch(std::bad_alloc &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

double *vec_mat_funcs::zero_matrix(int rows, int columns)
{	
	// this function will return a matrix that can hold row*cols elements of type double
	// use dynamic memory allocation to create a matrix
	// row indexing starts at i = 0 and ends at i = rows-1
	// column indexing starts at j = 0 and ends at j = columns-1
	// R. Sheehan 31 - 1 - 2013

	try{
		
		double *mat = zero_vector( (rows*columns) ); 

		return mat; 

		// bad_alloc thrown by default if there is a problem with the dimensions, 
		// or not enough memory is available

	}
	catch(std::bad_alloc &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

double *vec_mat_funcs::identity_matrix(int size)
{
	try{

		if(size > 0){

			double *arr = zero_matrix(size, size); 

			int i=0; 
			while(i<size){
				*(arr+i*size+i) = 1.0; 
				i++; 
			}

			return arr; 		
		}
		else{
			std::string reason = "Error: double *vec_mat_funcs::identity_matrix(int size)\n"; 
			reason += "size = " + template_funcs::toString(size) + " is not positive\n"; 
			throw std::invalid_argument(reason); 
			return nullptr; 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void vec_mat_funcs::print_vector(double *vec, int size)
{
	// print a vector to the screen
	// R. Sheehan 31 - 1 - 2013

	if(vec != nullptr){

		//std::cout<<"\nYour vector is\n"; 
		for(int i=0; i<size; i++){
			std::cout<<*(vec+i)<<"\n";
		}
		std::cout<<"\n";

	}

}

void vec_mat_funcs::print_matrix(double *mat, int rows, int columns)
{
	// print a matrix to the screen
	// R. Sheehan 31 - 1 - 2013

	if(mat != nullptr){

		//std::cout<<"\nYour matrix is\n"; 
		for(int i=0; i<rows; i++){
			for(int j=0; j<columns; j++)
				std::cout<<*(mat + i*columns + j)<<" ";
			std::cout<<"\n";
		}
		std::cout<<"\n";

	}
}

double *vec_mat_funcs::read_vector_from_file(std::string filename, int &n_data)
{
	// read a column of data from the file with name filename
	// store the data in the array arr
	// this function returns the data as a correctly sized array

	try{

		std::ifstream thefile;
		thefile.open(filename, std::ios_base::in);
	
		if(thefile.is_open()){
			// Since you are assuming a single column of numerical data you can use the stream extraction operator
			
			n_data=-1; // for some reason an extra \n is being added to the end single column files

			// http://www.cplusplus.com/reference/istream/istream/ignore/
			// istream& ignore (streamsize n = 1, int delim = EOF);
			// Extracts characters from the input sequence and discards them, until either n characters have been extracted, 
			// or one compares equal to delim.

			char endline = '\n'; 

			// Count the number of lines in the file
			while(thefile.ignore(1280, endline)){
				n_data++;
			}

			thefile.clear(); // empty the buffer
			thefile.seekg(0,std::ios::beg); // move to the start of the file

			// allocate the memory needed to store the data
			double *arr = new(double [n_data]); 

			// Loop over the lines and read the data into memory
			for(int i=0; i < n_data; i++){
				thefile>>arr[i];
			}
		
			thefile.close();

			return arr; 
		}
		else{
			std::string reason = "Error: void vec_mat_funcs::read_vector_from_file(std::string filename, int &n_data, double *arr)\n"; 
			reason += "File " + filename + " could not be opened\n"; 
			throw std::invalid_argument(reason); 
			return nullptr; 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void vec_mat_funcs::read_vector_from_file(std::string filename, int n_data, double *arr)
{
	// read a column of data from the file with name filename
	// store the data in the array arr
	// memory for arr must be allocated, this is not a very effective solution since you must know how much data is in the file
	// being read, which is not always the case
	// change function so that it returns an appropriately sized array

	try{

		std::ifstream thefile;
		thefile.open(filename, std::ios_base::in);
	
		if(thefile.is_open()){
			// Since you are assuming a single column of numerical data you can use the stream extraction operator
			// The number of data points must be known a-priori

			// Loop over the lines and read the data into memory
			for(int i=0; i < n_data; i++){
				thefile>>arr[i];
			}
		
			thefile.close();
		}
		else{
			std::string reason = "Error: void vec_mat_funcs::read_vector_from_file(std::string filename, int &n_data, double *arr)\n"; 
			reason += "File " + filename + " could not be opened\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}