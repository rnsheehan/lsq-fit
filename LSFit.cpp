#ifndef ATTACH_H
#include "Attach.h"
#endif

// Method Definitions for the lsfit class
// R. Sheehan 11 - 5 - 2009

lsfit::lsfit()
{
	// default constructor
	
	N = m = r = nu = 0;

	chi_sqr = chi_nu = R_sqr = 0.0;

	include_errors = false; 

	x = nullptr; y = nullptr; yerr = nullptr; 
	beta = nullptr; alpha = nullptr; covar = nullptr;   
}

lsfit::lsfit(std::string ord_file, std::string abs_file, int deg, bool assume_errors)
{
	// Use this constructor when reading the data from a file
	// the ordinate data, x-values, are held in a file named ord_file
	// the abscissa data, y-values, are held in a file named abs_file
	// the data from each file can be read using the read_vector_from_file function

	try{

		bool c1 = ord_file != "" ? true : false; 
		bool c2 = abs_file != "" ? true : false; 
		bool c3 = deg > 0 ? true : false; 

		if(c1 && c2 && c3){

			// This approach is wasterful as it requires you to open, read, close each file twice. 
			//N = std::min( useful_funcs::count_lines(ord_file), useful_funcs::count_lines(abs_file) ); 
			//x = vec_mat_funcs::zero_vector(N); 
			//y = vec_mat_funcs::zero_vector(N); 
			//vec_mat_funcs::read_vector_from_file(ord_file, N, x); // read the horizontal data and the number of data points
			//vec_mat_funcs::read_vector_from_file(abs_file, N, y); // read the vertical data and the number of data points

			int n1, n2; // use temporary variables to determine the size of the arrays, 
			// the user may attempt to input files which contain unequal numbers of data points

			x = vec_mat_funcs::read_vector_from_file(ord_file, n1); 
			y = vec_mat_funcs::read_vector_from_file(abs_file, n2); 

			N = std::min(n1, n2); // the size that is referenced throughout is the smaller, thereby avoiding potential beyond array access errors

			std::cout<<"N = "<<N<<"\n"; 

			m = deg; // specify the degree of the fitting polynomial
			r=m+1; 
			nu=N-r; 

			// Allocate the required memory
			alpha = vec_mat_funcs::zero_matrix(r, r); 
			beta = vec_mat_funcs::zero_vector(r); 
			covar = vec_mat_funcs::identity_matrix(r);

			include_errors = assume_errors; 

			if(include_errors){

				// include errors in the measured data in fit calculation
				yerr = vec_mat_funcs::zero_vector(N); 
				for(int i=0; i<N; i++){
					// Assume that error is sqrt of fitting value
					yerr[i] = sqrt(y[i]); 
				}

			}

			report_name = "Fitting_Report.txt";
		}
		else{
			std::string reason = "Error: lsfit::lsfit(std::string ord_file, std::string abs_file, int deg)\n"; 
			if(!c1) reason += ord_file + " is not a valid file name\n"; 
			if(!c2) reason += abs_file + " is not a valid file name\n"; 
			if(!c3) reason += "Degree: " + template_funcs::toString(deg) + " is not a valid value\n"; 

			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

lsfit::lsfit(int deg, int ndat, double *xdat, double *ydat, double *yerdat)
{
	// Use this constructor when data is defined in the code
	// xdat, ydat are pointers to dynamically allocated that hold the data in question

	try{

		bool c1 = deg > 0 && deg < ndat - 1 ? true : false;
		bool c2 = ndat > 2 ? true : false;
		bool c3 = xdat != nullptr ? true: false; 
		bool c4 = ydat != nullptr ? true: false; 
		bool c5 = yerdat != nullptr ? true: false; 

		if(c1 && c2 && c3 && c4 && c5){
			
			m=deg; // specify the degree of the fitting polynomial

			r=m+1;

			N=ndat;

			nu=N-r; // degrees of freedom

			x=xdat;

			y=ydat;

			include_errors = true; 

			yerr=yerdat;

			// Allocate the required memory
			alpha = vec_mat_funcs::zero_matrix(r, r); 
			beta = vec_mat_funcs::zero_vector(r); 
			covar = vec_mat_funcs::identity_matrix(r);

			report_name = "Fitting_Report.txt";
		}
		else{
			std::string reason = "Error: lsfit::lsfit(int deg, int ndat, double *xdat, double *ydat, double *yerdat)\n"; 
			if(!c1) reason += "Degree: " + template_funcs::toString(deg) + " is not a valid value\n"; 
			if(!c2) reason += "Number of data points: " + template_funcs::toString(ndat) + " is not a valid value\n"; 
			if(!c3) reason += "xdat contains no data\n";  
			if(!c4) reason += "ydat contains no data\n"; 
			if(!c5) reason += "yerdat contains no data\n"; 

			throw std::invalid_argument(reason);
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

lsfit::~lsfit()
{
	// destructor

	if(covar != nullptr) delete[] covar;

	if(alpha != nullptr) delete[] alpha; 

	if(beta != nullptr) delete[] beta; 

	if(x != nullptr) delete[] x;

	if(y != nullptr) delete[] y;

	if(include_errors && yerr != nullptr) delete[] yerr; 
}

void lsfit::make_alpha()
{
	// This function calculates the matrix to be used in solving for the parameters
	// The matrix has size r*r

	//alpha = vec_mat_funcs::zero_matrix(r, r); 

	try{

		if(alpha != nullptr){

			for(int l = 0; l < r; l++){
				for(int k = 0; k < r; k++){

					if(l<=k){
						*(alpha + l*r + k) = alpha_elem(l,k);
					}
					else{
						*(alpha + l*r + k) = *(alpha + k*r + l); // The calculation is symmetric
					}
				}
			}			
		}
		else{
			std::string reason = "Error: void lsfit::make_alpha()\n"; 
			reason += "alpha has not been allocated\n"; 

			throw std::invalid_argument(reason);
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void lsfit::make_beta()
{
	//This function calculates the rhs vector to be used in solving for the parameters
	//The vector has size r

	//beta = vec_mat_funcs::zero_vector(r); 

	try{

		if(beta != nullptr){
			
			for(int k = 0; k < r; k++){
				*(beta + k) = beta_elem(k);
			}

		}
		else{
			std::string reason = "Error: void lsfit::make_beta()\n"; 
			reason += "beta has not been allocated\n"; 

			throw std::invalid_argument(reason);
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void lsfit::perform_fit()
{
	//Do everything required for a fit

	fit_data();

	chi_sqr_sum();

	R_sqr_sum();

	chi_sqr_probability();

	residuals();
}

void lsfit::set_degree(int deg)
{
	//This specifies the degree of the fitting polynomial and some other parameters
	//Can be used to change the degree etc
	//deg = 1 => line
	//deg = 2 => quadratic
	//deg = 3 => cubic.....

	m=deg;
	r=m+1;
	nu=N-r;
}

void lsfit::fit_data()
{
	//This function performs the calculations necessary for determining the fit parameters
	// The main calculation is the solution of the system of linear equations alpha.a = beta
	// since alpha symmetric and positive definite this is done using Cholesky decomposition
	// which is implemented in LAPACK, you could use dposv, but since you also want to determine the
	// co-variance matrix you need to store the Cholesky decomposition locally in alpha
	// Hence use dpotrf followed by two calls to dpotrs so that you only have to compute the Cholesky decomposition
	// of alpha once. 
	// This will not save you a whole lot of time since the systems of equations will most often only be of size 2*2, 
	// for linear fits, or 3*3 for quadratic fits. 

	try{

		bool c1 = r > 1 ? true : false; 
		bool c2 = alpha != nullptr ? true : false; 
		bool c3 = beta != nullptr ? true : false; 

		if(c1 && c2 && c3){

			// 1. Calculate the fit matrices alpha and beta	
			make_alpha();
	
			make_beta();

			std::cout<<"The fitting matrix is \n";
			vec_mat_funcs::print_matrix(alpha,r,r);
	
			std::cout<<"The fitting vector is \n";
			vec_mat_funcs::print_vector(beta,r);

			char UPLO = 'L'; // access lower triangular portion	
			int INFO = 0; // exit code
			int nrhs = 1; 

			// 1. Compute the Cholesky Decomposition of alpha, this is used in 
			dpotrf_(&UPLO, &r, alpha, &r, &INFO); 

			if(INFO < 0){
				std::cerr<<"Error: void dpotrf_(const char *uplo, const int *N, double *A, const int *lda, int *info)\n";
				std::cerr<<"INFO = "<<INFO<<": i^{th} argument had an illegal value\n"; 
			}
			else if(INFO > 0){
				std::cerr<<"Error: void dpotrf_(const char *uplo, const int *N, double *A, const int *lda, int *info)\n";
				std::cerr<<"INFO = "<<INFO<<": the leading minor of order i is not positive definite, and the factorization could not be completed\n"; 
			}
			else{
				// INFO  == 0
				// Successful calculation of the Cholesky Decomposition
			
				// 2. Solve the system alpha.a=beta using Cholesky Decomposition of alpha
				// alpha now contains its own Cholesky Decomposition
				dpotrs_(&UPLO, &r, &nrhs, alpha, &r, beta, &r, &INFO);

				if(INFO != 0){
					std::cerr<<"INFO = "<<INFO<<": the i-th argument had an illegal value when calculating beta\n";
				}

				// 3. Calculate the covariance matrix by computing the inverse of alpha using its Cholesky Decomposition
				// covar is assigned to be the order r identity matrix
				// alpha now contains its own Cholesky Decomposition
				dpotrs_(&UPLO, &r, &r, alpha, &r, covar, &r, &INFO);

				if(INFO != 0){
					std::cerr<<"INFO = "<<INFO<<": the i-th argument had an illegal value when calculating covar\n";
				}

				std::cout<<"The error matrix is \n";
				vec_mat_funcs::print_matrix(covar,r,r);

				std::cout<<"Upon solution of this system of equations\n"; 
				std::cout<<"The fitting parameters are \n";
				for(int i = 0; i < r; i++){
					std::cout<<"a_"<< i <<" = "<<beta[i]<<" +/- "<<sqrt(*(covar + i*r + i))<<"\n";
				}
				std::cout<<"\n";
			}		
		}
		else{
			std::string reason = "Error: void lsfit::fit_data()\n"; 
			if(!c1) reason += "r: " + template_funcs::toString(r) + " is not a valid value\n"; 
			if(!c2) reason += "alpha has not been allocated\n"; 
			if(!c3) reason += "beta has not been allocated\n"; 

			throw std::invalid_argument(reason);
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void lsfit::chi_sqr_probability()
{
	// Calculates the chi^{2} probability associated with the fit

	try{

		bool c1 = nu > m ? true : false; // degrees of freedom must be greater than degree of fitting polynomial
		bool c2 = chi_sqr > 0 ? true : false; // chi^{2} parameter for the fit has not been computed

		if(c1 && c2){

			double prob;

			prob = special::gammq(0.5*nu, 0.5*chi_sqr );
	
			std::cout<<"\n";
			std::cout<<"The probability of exceeding the chi-squared value \n";
			std::cout<<"for the fit in some measurement is "<<100.0*prob<<" %\n";
			std::cout<<"\n";
		
		}
		else{
			std::string reason = "Error: void lsfit::chi_sqr_probability()\n"; 
			if(!c1) "Degrees of freedom: " + template_funcs::toString(nu) + " is not valid\n"; 
			if(!c2) "chi^{2}: " + template_funcs::toString(chi_sqr, 2) + " is not valid\n"; 
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
	}
}

void lsfit::chi_sqr_sum()
{
	// Computes the value of chi^{2} for the fit

	try{

		bool c1 = r > 1 ? true : false;
		bool c2 = nu > m ? true : false; // degrees of freedom must be greater than degree of fitting polynomial
		bool c3 = x != nullptr ? true : false; 
		bool c4 = y != nullptr ? true : false; 
		bool c5 = beta != nullptr ? true : false; 

		if(c1 && c2 && c3 && c4 && c5){

			double tt,pval;

			tt=chi_sqr=0.0;

			for(int i = 0; i < N; i++){

				pval = useful_funcs::poly(x[i], beta, r); // P( x_{i} ) value of the data according to the fitting model

				//tt = ( y[i] - pval )/(yerr[i]); // y_{i} - P( x_{i} )
				tt = (include_errors ? ( y[i] - pval )/(yerr[i]) : ( y[i] - pval ) ); // y_{i} - P( x_{i} )

				chi_sqr += template_funcs::DSQR(tt); // compute the chi^{2} parameter for the fit
			}

			chi_nu = (chi_sqr/nu);

			std::cout<<"The chi-squared value for the fit = "<<chi_sqr<<"\n";
			std::cout<<"The number of degrees of freedom for the fit = "<<nu<<"\n";
			std::cout<<"chi-squared / nu = "<<chi_nu<<"\n";
			std::cout<<"\n";
		
		}
		else{
			std::string reason = "Error: void lsfit::chi_sqr_sum()\n"; 
			if(!c1) "r: " + template_funcs::toString(r) + " is not valid\n"; 
			if(!c2) "Degrees of freedom: " + template_funcs::toString(nu) + " is not valid\n"; 
			if(!c3) "x has not been assigned\n"; 
			if(!c4) "y has not been assigned\n"; 
			if(!c5) "beta has not been assigned\n"; 
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
	}
}

void lsfit::R_sqr_sum()
{
	// Compute the R^{2} coefficient for the fit

	try{

		bool c1 = N > 2 ? true : false; 
		bool c2 = chi_sqr > 0 ? true : false;  
		bool c3 = y != nullptr ? true : false;

		if(c1 && c2 && c3){

			double t,s1,s2,sst;

			t=sst=s1=s2=0.0;

			for(int i = 0; i < N; i++){

				//t = y[i]/(yerr[i]);
				t = (include_errors ? y[i]/(yerr[i]) : y[i]);

				s1 += template_funcs::DSQR(t);

				s2 += t;
			}

			sst = s1 - ( template_funcs::DSQR(s2) / N );

			R_sqr = 1.0 - (chi_sqr/sst);

			std::cout<<"The R-squared value for the fit is "<<R_sqr<<"\n";
		
		}
		else{
			std::string reason = "Error: void lsfit::R_sqr_sum()\n"; 
			if(!c1) "N: " + template_funcs::toString(N) + " is not valid\n"; 
			if(!c2) "chi^{2}: " + template_funcs::toString(chi_sqr, 3) + " is not valid\n"; 
			if(!c3) "y has not been assigned\n"; 
			
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
	}
}

void lsfit::residuals()
{
	try{
			
		bool c1 = N > 2 ? true : false; 
		bool c2 = x != nullptr ? true : false;
		bool c3 = y != nullptr ? true : false;

		if(c1 && c2 && c3){

			std::cout<<"The residuals for the fitted data r = y_{i} - P(x_{i})\n"; 
			for(int i = 0; i < N; i++){
				std::cout<<x[i]<<" , "<<y[i] - useful_funcs::poly(x[i],beta,r)<<"\n";
			}
			std::cout<<"\n";
		
		}
		else{
			std::string reason = "Error: void lsfit::residuals()\n"; 
			if(!c1) "N: " + template_funcs::toString(N) + " is not valid\n"; 
			if(!c2) "x has not been assigned\n"; 
			if(!c3) "y has not been assigned\n"; 
			
			throw std::invalid_argument(reason); 		
		}	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
	}
}

void lsfit::write_report()
{
	// Create a report for the computed data
	// R. Sheehan 15 - 7 - 2014

	try{
		// open the file for writing
		report.open(report_name,std::ios_base::out|std::ios_base::trunc);

		if(report.is_open()){
			std::string time = useful_funcs::TheTime();

			report<<"Fit completed "<<time<<"\n";
			report<<"The coefficients of a fitting\npolynomial of degree m = "<<m<<" were computed\n"; 
			report<<"( m = 1 => linear fit, m = 2 => quadratic fit, etc. )\n\n";
			report<<"Fit polynomial is of the form \\sum_{i=0}^{m} a_{i} x^{i}\n\n";

			report<<"The fitting parameters in this case are \n";
			for(int i = 0;i < r; i++){
				report<<"a_"<<(i)<<" = "<<std::setprecision(15)<<beta[i]<<" +/- "<<sqrt( *(covar +i*r + i) )<<"\n";
			}
			report<<"\n";

			report<<"The co-variance matrix for the fit is\n";
			for(int i = 0; i< r; i++){
				for(int j = 0; j< r; j++)
					report<<*(covar +i*r + j)<<" ";
				report<<"\n";
			}
			report<<"\n";

			report<<"The chi-squared value for the fit = "<<chi_sqr<<"\n";
			report<<"The number of degrees of freedom for the fit = "<<nu<<"\n";
			report<<"chi-squared / nu = "<<chi_nu<<"\n";
			report<<"\n";

			report<<"The R-squared value for the fit is "<<R_sqr<<"\n";
			report<<"\n";

			report<<"The probability of exceeding the chi-squared value \n";
			report<<"for the fit in some measurement is "<<100.0*(special::gammq(0.5*nu, 0.5*chi_sqr ))<<" %\n";
			report<<"\n";

			report<<"The residuals for the fitted data \nassuming r = y_{i} - P(x_{i})\n"; 
			for(int i = 0; i < N; i++){
				report<<x[i]<<"\t"<<y[i] - useful_funcs::poly(x[i],beta,r)<<"\n";
			}
			report<<"\n";

			report<<"This code has been brought to you by Robert Sheehan, BSc (Hons), PhD\n";
			report<<"For more information please contact robernsheehan@gmail.com\n";
			report.close(); 

		}
		else{
			std::string reason = "Error: void lsfit::write_report()\n"; 
			reason += "Could not open file " + report_name + "\n";

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lsfit::write_coeffs()
{
	try{
		// open the file for writing
		report.open("Fit_Coefficients.txt",std::ios_base::out|std::ios_base::trunc);

		if( report.is_open() ){

			for(int i = 0; i < r; i++){
				report<<std::setprecision(15)<<beta[i]<<"\n";
			}

			report.close(); 
		}
		else{
			std::string reason = "Error: void lsfit::write_report()\n"; 
			reason += "Could not open file " + report_name + "\n";

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// private methods

double lsfit::beta_elem(int k)
{
	// This returns the sum that belongs in element k of the vector beta

	try{

		bool c1 = N > 2 ? true : false; 
		bool c2 = x != nullptr ? true : false;
		bool c3 = y != nullptr ? true : false;

		if(c1 && c2 && c3){
		
			double t,s;
			t=0.0; // term in a sum
			s=0.0; // sum variable

			for(int j = 0; j < N; j++){
				if(k==0){
					//t = ( y[j] / template_funcs::DSQR( yerr[j] ) );
					t = ( include_errors ? ( y[j] / template_funcs::DSQR( yerr[j] ) ) : y[j] ); 
				}
				else if(k==1){
					//t = ( y[j] / template_funcs::DSQR( yerr[j] ) )*x[j];
					t = (include_errors ? ( y[j] / template_funcs::DSQR( yerr[j] ) )*x[j] : y[j]*x[j] ); 
				}
				else{
					//t = ( y[j] / template_funcs::DSQR( yerr[j] ) )*pow(x[j],k);
					t = ( include_errors ? ( y[j] / template_funcs::DSQR( yerr[j] ) )*pow(x[j],k) : y[j]*pow(x[j],k));
				}

				s+=t; 
			}

			return s;
		}
		else{
			std::string reason = "Error: double lsfit::beta_elem(int k)\n"; 
			if(!c1) "N: " + template_funcs::toString(N) + " is not valid\n"; 
			if(!c2) "x has not been assigned\n"; 
			if(!c3) "y has not been assigned\n"; 
			
			throw std::invalid_argument(reason); 	

			return -1; 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double lsfit::alpha_elem(int l, int k)
{
	// This returns the sum that belongs to element (l,k) of the matrix alpha

	try{

		bool c1 = N > 2 ? true : false; 
		bool c2 = x != nullptr ? true : false;
		bool c3 = y != nullptr ? true : false;

		if(c1 && c2 && c3){
		
			double t,s;

			t=s=0.0;

			for(int i = 0;i< N;i++){
				if(l==0 && k==0){
					//t = 1.0 / ( yerr[i]*yerr[i] );
					t = ( include_errors ? 1.0 / ( template_funcs::DSQR(yerr[i]) ) : 1.0 ); 
				}
				else{
					//t = ( pow(x[i],l + k) ) / ( template_funcs::DSQR(yerr[i]) );
					t = ( include_errors ? ( pow(x[i],l + k) ) / ( template_funcs::DSQR(yerr[i]) ) : pow( x[i], l + k) ); 
				}

				s += t; 
			}

			return s;

		}
		else{
			std::string reason = "Error: double lsfit::alpha_elem(int k)\n"; 
			if(!c1) "N: " + template_funcs::toString(N) + " is not valid\n"; 
			if(!c2) "x has not been assigned\n"; 
			if(!c3) "y has not been assigned\n"; 
			
			throw std::invalid_argument(reason); 	

			return -1; 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}