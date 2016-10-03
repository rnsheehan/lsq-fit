#ifndef ATTACH_H
#include "Attach.h"
#endif

void fitting_object_example(); 

int main(int argc, char *argv[])
{
	try{
	
		if(argc >= 4){
			
			// List off the input parameters
			// Program needs 5 or more parameters to run, remember that the name of the program is also considered a parameter
			// argv[0] = program name
			// argv[1] = name of file containing horizontal data points (ordinates)
			// argv[2] = name of file containing vertical data points (abcissae)
			// argv[3] = degree of fit to be computed m = 1 => linear, m = 2 => quadratic
			// argv[4] = location at which data is stored
			
			std::cout<<argc-1<<" parameters were input into the program\n"; 
			for(int count = 1; count < argc; count++){
				std::cout<<"argv["<<count<<"] = "<<argv[count]<<"\n"; 
			}
			std::cout<<"\n";

			// Assign the values of the arguments to the input parameters for the fitting object
			std::string x_file(argv[1]); 
			std::string y_file(argv[2]);
			std::string storage(argv[4]); 

			int fit_degree = atoi(argv[3]); 

			bool dir_exists = false; 

			useful_funcs::set_directory(storage, dir_exists); 

			if(dir_exists){

				lsfit the_fit(x_file, y_file, fit_degree); // instantiate the fitting class
	
				// to perform the fit to the data set we use the fit_data method
				// which is defined in the LSFit.cpp source file
				the_fit.fit_data(); 

				// With the lsfit object we can also determine a measure of the quality of the fit
				// This is measured by the R^{2} coefficient, if R^{2} ~ 1 then the fit is good
				// R^{2} is computed using the total error, often called the chi^{2} sum
				// and another value called the SST, R^{2} = 1 - chi^{2} / SST
				the_fit.chi_sqr_sum(); // this is a measure of the total error

				the_fit.R_sqr_sum(); // this computes the R^{2} coefficient once chi^{2} has been computed

				the_fit.chi_sqr_probability(); // Calculates the chi^{2} probability associated with the fit

				the_fit.residuals(); // Obtain the residuals, differences between measured and predicted values

				the_fit.write_coeffs(); // print the fitted coefficients to a file

				the_fit.write_report(); // Output the results to a file
			}
			else{
				std::string reason = "Error: Poly_Fit\n";
				reason += "Named storage directory: " + storage + " does not exist\n"; 

				throw std::invalid_argument(reason); 
			}
		}
		else{
			std::string reason = "Error: Poly_Fit\n";
			reason += "Insufficient number of input arguments\n"; 

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	
	return 0; 
}

void fitting_object_example()
{
	// This function demonstrates the features available to the fitting object lsfit
	// This class is used to compute a fitting polynomial for some data set
	// where the error is minimised in the least squares sense. 
	// For Theory see "Data Reduction and Error Analysis" by P. R. Bevington and D. K. Robinson, Ch. 6 and 7.
	// or any standard statistics text. 
	// R. Sheehan 15 - 7 - 2014

	// Use the stl::string to name the files
	std::string x_file = "X_Data.txt"; 
	std::string y_file = "Y_Data.txt"; 

	// Instantiate an lsfit object
	// The constructor for the class requires that we tell the class
	// what files we're reading the data from and also the degree of the fitting polynomial
	
	int fit_degree = 1; // start with a linear fit
	
	// these variables could be modified to be input from keyboard

	lsfit the_fit(x_file, y_file, fit_degree); // instantiate the fitting class
	
	// to perform the fit to the data set we use the fit_data method
	// which is defined in the LSFit.cpp source file
	the_fit.fit_data(); 

	// With the lsfit object we can also determine a measure of the quality of the fit
	// This is measured by the R^{2} coefficient, if R^{2} ~ 1 then the fit is good
	// R^{2} is computed using the total error, often called the chi^{2} sum
	// and another value called the SST, R^{2} = 1 - chi^{2} / SST
	the_fit.chi_sqr_sum(); // this is a measure of the total error

	the_fit.R_sqr_sum(); // this computes the R^{2} coefficient once chi^{2} has been computed

	the_fit.chi_sqr_probability(); // Calculates the chi^{2} probability associated with the fit

	the_fit.residuals(); // Obtain the residuals, differences between measured and predicted values

	the_fit.write_report(); // Output the results to a file

	// It would be a small matter to adapt this program for use in a command line setting
	// For more details please contact me via email at r.sheehan@gmx.com
}