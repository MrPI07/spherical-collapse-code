#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>

using namespace std;

namespace integration {

// Function prototype
double integrateSmp(double a, double b, double (*f)(double), int n);
double integrateTz(double a, double b, double y, double (*f)(double, double), int n);
void integrateTz2(double a, double b, double (*f)(double), int n, double& integral);

// Función de integración de Simpson
double integrateSmp(double a, double b, double (*f)(double), int n) {
    double h = (b-a)/n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        sum += 2*f(x + h/2.0);
        sum += 4*f(x);
    }

    sum *= h/6.0;
    return h*sum;
}

// Función de integración del trapecio
double integrateTz0(double a, double b, double (*f)(double), int n) {
    double h = (b-a)/n;
    double sum = 0.5*(f(a) + f(b));
    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        sum += f(x);
    }
    return h*sum;
}

double integrateTz(double a, double b, double y, double (*f)(double, double), int n) {
    double h = (b-a)/n;
    double sum = 0.5*(f(a, y) + f(b, y));
    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        sum += f(x, y);
    }
    return h*sum;
}

void integrateTz2(double a, double b, double (*f)(double), int n, double& integral) {
    double h = (b-a)/n;
    double sum = 0.5*(f(a) + f(b));
    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        sum += f(x);
    }
    integral = h*sum;
}


} //integration