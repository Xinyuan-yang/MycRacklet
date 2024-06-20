#ifndef FIND_STABLE_ALPHA_HH
#define FIND_STABLE_ALPHA_HH

// Include necessary libraries
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <cmath> // For sqrt, pow, etc.

struct Parameters {
    double T;
};

struct integration_params {
    double lambda;
    Parameters params;
};
// Function prototypes
double Phi(double lambda, double x);
double integrand(double x, void *p);
double func(double lambda);
double derivative_func(double lambda);
double newton_raphson(double initial_guess, double tolerance, int max_iterations, Parameters &params);

#endif // MYSOLVE_HH