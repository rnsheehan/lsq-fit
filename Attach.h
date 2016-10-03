#ifndef ATTACH_H
#define ATTACH_H

#include <cstdlib>
#include <iostream>
#include <iomanip>

// need these for directory manipulation
#include <direct.h>
#include <errno.h>

#include <string>
#include <sstream>
#include <fstream>

#include <complex>
#include <cmath>

// Constants
static const double EPS=(3.0e-12);

static const double p=(atan(1.0));
static const double Two_PI=(8.0*p);
static const double PI=(4.0*p);
static const double PI_2=(2.0*p);
static const double PI_3=((4.0/3.0)*p);
static const double PI_4=(p);
static const double PI_5=((4.0/5.0)*p);
static const double PI_6=((2.0/3.0)*p);

//static const std::complex<double> one(1.0, 0.0); 

//static const std::string dottxt=".txt";

#include "Useful.h"
#include "Templates.h"
#include "Special_Functions.h"

#include "Vec_Mat_Lib.h"
#include "LAPACK_Header.h"

#include "LSFit.h"

#include "Examples.h"

#endif