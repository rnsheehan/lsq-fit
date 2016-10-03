#ifndef ATTACH_H
#include "Attach.h"
#endif

std::string useful_funcs::TheTime()
{
	// Implementation of a function returning the current time as a string
	// This is just a neater way of ensuring that the time can be correctly and easily accessed
	// without being hassled about whether or not you've remembered to use non-deprecated versions 
	// of certain functions
	// R. Sheehan 4 - 7 - 2011
	
	const int N=30;	
	char time_str[N];	
	size_t bytes=( N*sizeof(char) );
	
	time_t rawtime;
	
	struct tm timeinfo;
	struct tm *timeinfo_ptr;
	
	timeinfo_ptr=&timeinfo;
	
	// Get current time information
	time(&rawtime);
	
	localtime_s(timeinfo_ptr,&rawtime);
	
	asctime_s(time_str,bytes,timeinfo_ptr);
	
	// Deprecated calls
	//timeinfo=localtime(&rawtime);
	//asctime(timeinfo);
	
	std::string the_time;
	the_time.append(time_str);
	
	return the_time;
}

int useful_funcs::count_lines(std::string filename)
{
	// Count the number of lines in a file containing data

	try{

		std::ifstream thefile;
		thefile.open(filename, std::ios_base::in);
	
		if(thefile.is_open()){				
			
			int n_data=-1; // for some reason an extra \n is being added to the end single column files

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
		
			thefile.close();

			return n_data; 
		}
		else{
			std::string reason = "Error: int useful_funcs::count_lines(std::string filename)\n"; 
			reason += "File " + filename + " could not be opened\n"; 
			throw std::invalid_argument(reason);
			return -1; 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void useful_funcs::exit_failure_output(std::string reason)
{
	// Code that creates a file and writes a reason in it why the program crashed
	// If it is called of course
	// Call before using the exit(EXIT_FAILURE) command

	// This function outputs to a file an explanation of why the program exited with an EXIT_FAILURE
	// R. Sheehan 17 - 5 - 2011
	
	// Get current time information
	std::string time = TheTime();

	std::ofstream write; // open file for writing
	
	write.open("Exit_Failure_Explanation.txt",std::ios_base::out|std::ios_base::trunc);
	
	//if(!write){
	//	std::cout<<"You're not going to see this statement\n";
	//	std::cout<<"\n";
	//}
	//else{
	//	//printf ( "Current local time and date: %s", asctime (timeinfo) );
	//	write<<"Program Exit Explanation\n\n";
	//	write<<"Error occurred "<<time<<"\n";
	//	write<<reason<<"\n";
	//	write.close();
	//}

	if( write.is_open() ){
		
		write<<"Program Exit Explanation\n\n";
		write<<"Error occurred: "<<time<<"\n";
		write<<reason<<"\n";

		write.close();
	}
}

void useful_funcs::set_directory(std::string &dir_name, bool &dir_set)
{
	// Set the current working directory
	// _chdir return a value of 0 if successful. 
	// A return value of –1 indicates failure. If the specified path could not be found, errno is set to ENOENT. 
	// If dirname is NULL, the invalid parameter handler is invoked
	// R. Sheehan 6 - 8 - 2012
	
	if(_chdir( dir_name.c_str() ) ){
		switch (errno){
			case ENOENT:
				//printf( "Unable to locate the directory: %s\n", dir_name );
				std::cout<<"Unable to locate the directory: "<<dir_name<<"\n";
				break;
			case EINVAL:
				std::cout<<"Invalid buffer.\n";
				break;
			default:
				std::cout<<"Unknown error.\n";
		}
		
		dir_set = false; 
	}
	else{
		
		//cout<<"Directory has been changed\n"; 

		dir_set = true; 
	}
}

double useful_funcs::poly(double x, double *b, int m)
{
	// This evaluates the polynomial of degree m at the point x
	// The coefficients of the polynomial are stored in the vector b
	
	double p;

	p = b[m]; // initialise the value of p

	// loop over the remaining coefficients
	for(int j=m-1;j>=0;j--){// For zero-offset arrays replace j>0 by j>=0

		p = p*x + b[j];
	}

	return p;
}