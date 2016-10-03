#ifndef LSFIT_H
#define LSFIT_H

// R. Sheehan 11 - 5 - 2009
// Class defining the functions required to perform a least squares fit to a set of data
// For Theory see "Data Reduction and Error Analysis" by P. R. Bevington and D. K. Robinson, Ch. 6 and 7.

class lsfit{
	// Methods
public:
	lsfit(); // Default constructor
	lsfit(std::string ord_file, std::string abs_file, int deg, bool assume_errors = false); // Use this constructor when reading the data from a file
	lsfit(int deg, int ndat, double *xdat, double *ydat, double *yerdat); // Use this constructor when the data is defined in code
	~lsfit(); // destructor

	void perform_fit(); // perform all the steps necessary to generate the fit
	void chi_sqr_sum(); // Calculate the chi-squared value from the fit
	void chi_sqr_probability(); // Probablity that the observed chi-square will exceed the chi-squared value for the correct fit
	void R_sqr_sum(); // Compute the R^{2} coefficient for the fit
	void set_degree(int deg); // set the degree of the fitting polynomial
	void fit_data(); // Generate the fitting matrix and rhs vector, solve the system for the fitting parameters	
	void residuals(); // Compute the residuals associated with the fit
	void write_report(); // print the results to a file
	void write_coeffs(); // print the fitting coefficients to a file

private:
	void make_alpha();//Creates the matrix to be used in determining the fit
	void make_beta();//Creates the rhs vector to be used in determining the fit

	double poly(double x, double *b, int m); // returns the value of a polynomial of degree m whose coefficients are stored in b

	double beta_elem(int k); // return element k of the rhs fit vector
	double alpha_elem(int l, int k); // return element l, k of the fit matrix	
	
	// Class Members
private:
	int N; // Number of data points
	int m; // Degree of the interpolating polynomial
	int r; // r = m + 1, deg(3) poly => 4 parameters to be found
	int nu; // degrees of freedom of the fit

	bool include_errors; // perform fit taking statistical errors into account
	
	double chi_sqr; // chi-squared value computed from the fit
	double chi_nu; // Ratio of the chi-squared value and the degrees of freedom of the fit
	double R_sqr; // R^{2} coefficient, measures "goodness" of fit

	double *x; // vector holding the x-data
	double *y; // vector holding the y-data
	double *yerr; // vector holding the errors in the y-data
	double *beta; // RHS Fitting vector

	double *alpha; // Fitting Matrix
	double *covar; // Co-Variance Matrix

	std::ofstream report; // object for writing results to a file
	std::string report_name; // name of file to hold the results of the calculation
};

#endif