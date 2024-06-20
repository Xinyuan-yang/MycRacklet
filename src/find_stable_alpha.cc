#include "find_stable_alpha.hh"
#include <iostream> // For std::cout

double Phi(double lambda, double x) {
    gsl_set_error_handler_off();
    return gsl_sf_expint_E1(pow(lambda, 2) * pow(x, 2)) / sqrt(1 - pow(x, 2)) * x;
}

double integrand(double x, void *p) {
    integration_params *params = (integration_params *)p;
    return Phi(params->lambda, x);
}

double func(double lambda, Parameters &params) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    integration_params int_params = {lambda, params};

    gsl_function F;
    F.function = &integrand;
    F.params = &int_params;

    gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result - params.T;
}

double derivative_func(double lambda, Parameters &params) {
    const double delta = 1e-5;
    double f1 = func(lambda + delta, params);
    double f2 = func(lambda - delta, params);
    return (f1 - f2) / (2 * delta);
}

double newton_raphson(double initial_guess, double tolerance, int max_iterations, Parameters &params) {
    double lambda = initial_guess;
    for (int i = 0; i < max_iterations; i++) {
        double f = func(lambda, params);
        double df = derivative_func(lambda, params);
        double delta = f / df;
        lambda -= delta;
        if (std::abs(delta) < tolerance) {
            std::cout << "Root found: " << lambda << std::endl;
            break;
        }
    }
    return lambda;
}