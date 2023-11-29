#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include<fstream>

#include "rkf45.hpp"
#include "ode.hpp"
#include "integral.cpp"
#include "spline.hpp"

using namespace std;
using namespace integration;

namespace EdS{ //Modelo EdS

    int modo = 0;
    double z0 = 1e7;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    void Model_initialization();
    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga);

    double dEdaE(double a);
    double E(double a);

    double Om_m, Om_w, Om_r, Om_k, Om_b, w0;

    void modelo(){
        cout<<"\n \t\t Modelo EdS \n"<<endl;
        return;
    }

    void Model_initialization(){

        cout<<"\n \t\t Modelo EdS \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w);
        eos = -1;
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a));
        ga = pow(exp(a),-3.0);
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 1.0;
        r0_rho = 0;
        de0_rho = 0.0;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = -1;

        return;
    }

    double E2(double a){ //a = lna
        double E2;
        E2 = (Om_m / pow(exp(a),3)) + Om_w;
        return E2;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
        //double Om_m = 0.3;
        return Om_m / (pow(exp(a),3)*E2(a));
    }

    double Omega_w(double a){
        //double Om_w = 0.7;
        return Om_w / E2(a);
    }

    double w(double a){
        return -1;
    }

    double g(double a){
        return pow(exp(a),-3.0);
    }

}

namespace LCDM{ //Modelo LCDM

    int modo = 0;
    double z0 = 1e7;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    double dEdaE(double a);
    double E(double a);
    void Model_initialization();
    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga);
    
    void modelo(){
        cout<<"\n \t\t Modelo LCDM \n"<<endl;
        return;
    }
    
    double Om_m, Om_w, Om_r, Om_k, Om_b, w0;

    void Model_initialization(){

        cout<<"\n \t\t Modelo LCDM \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        //cambiar Exp_func por cosmologyFunc
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w);
        eos = -1;
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a));
        ga = 1.0;
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 0.3103;
        r0_rho = 0;
        de0_rho = 0.6897;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = -1;

        return;
    }

    double E2(double a){ //a = lna
        double E;
        E = (Om_m / pow(exp(a),3)) + Om_w*g(a);
        return E;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
        //double Om_m = 0.3;
        return (Om_m)/(pow(exp(a),3)*E2(a));
    }

    double Omega_w(double a){
        //double Om_w = 0.7;
        return (Om_w / E2(a))*g(a);
    }

    double w(double a){
        return -1;
    }

    double g(double a){
        return 1;
    }

}

namespace CPL{ //Modelo CPL

    int modo = 0;
    double z0 = 1e7;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    double dEdaE(double a);
    double E(double a);
    
    void modelo(){
        cout<<"\n \t\t Modelo CPL "<<endl;
        cout<<" \t w(a) = w0 + wa(1-a) \n"<<endl;
        return;
    }
    
    double Om_m, Om_b, Om_w, Om_r, Om_k, w0;

    void Model_initialization(){

        cout<<"\n \t\t Modelo CPL "<<endl;
        cout<<" \t w(a) = w0 + wa(1-a), w0=-0.75, wa=0.40. \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        //cambiar Exp_func por cosmologyFunc
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w*g(a));
        eos = w(a);
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a))*g(a);
        ga = g(a);
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 0.3103;
        r0_rho = 0;
        de0_rho = 0.6897;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = w(0);

        return;
    }

    double E2(double a){ //a = lna
        double E;
        E = (Om_m / pow(a,3)) + Om_w*g(a);
        return E;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
       
        return Om_m/(pow(a,3)*E2(a));
    }

    double Omega_w(double a){
      
        return (Om_w / E2(a))*g(a);
    }

    double w(double a){
        double w0=-0.75, wa=0.40, w;
        w = w0 + wa*(1.0-exp(a));
        return w;
    }

    double g(double a){
        double w0=-0.75, wa=0.40, gp;
        gp = (pow(exp(a),-3*(1.0+w0+wa)))*exp(-3.0*wa*(1.0-exp(a)));
        return gp;
    }

}

namespace CPLF{ //Modelo CPL Phanton

    int modo = 0;
    double z0 = 1e7;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    double dEdaE(double a);
    double E(double a);
    
    void modelo(){
        cout<<"\n \t\t Modelo CPL Phanton "<<endl;
        cout<<" \t w(a) = w0 + wa(1-a) \n"<<endl;
        return;
    }

    double Om_m, Om_b, Om_w, Om_r, Om_k, w0;

    void Model_initialization(){

        cout<<"\n \t\t Modelo CPL Phanton "<<endl;
        cout<<" \t w(a) = w0 + wa(1-a), w0=-1.1, wa=-1. \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        //cambiar Exp_func por cosmologyFunc
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w*g(a));
        eos = w(a);
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a))*g(a);
        ga = g(a);
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 0.3103;
        r0_rho = 0;
        de0_rho = 0.6897;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = w(0);

        return;
    }

    double E2(double a){ //a = lna
        double E;
        E = (Om_m / pow(a,3)) + Om_w*g(a);
        return E;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
        //double Om_m = 0.3;
        return Om_m/(pow(a,3)*E2(a));
    }

    double Omega_w(double a){
        //double Om_w = 0.7;
        return (Om_w / E2(a))*g(a);
    }

    double w(double a){
        double w0=-1.1, wa=-1.0, w;
        w = w0 + wa*(1-exp(a));
        return w;
    }

    double g(double a){
        double w0=-1.1, wa=-1.0, gp;
        gp = pow(exp(a),-3*(1+w0+wa))*exp(3*wa*(exp(a)-1));
        return gp;
    }

}

namespace neDE{ //Modelo de energía oscura temprana

    int modo = 0;
    double z0 = 1e7;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    double dEdaE(double a);
    double E(double a);
    
    void modelo(){
        cout<<"\n \t\t Modelo neDE \n"<<endl;
        return;
    }
    
    double Om_m, Om_b, Om_w, Om_r, Om_k, w0;

    void Model_initialization(){

        cout<<"\n \t\t Modelo neDE \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        //cambiar Exp_func por cosmologyFunc
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w*g(a));
        eos = w(a);
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a))*g(a);
        ga = g(a);
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 0.3103;
        r0_rho = 0;
        de0_rho = 0.6897;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = w(0);

        return;
    }

    double E2(double a){ //a = lna
        double E;
        E = (Om_m / pow(a,3)) + Om_w*g(a);
        return E;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
        //double Om_m = 0.3;
        return Om_m/(pow(a,3)*E2(a));
    }

    double Omega_w(double a){
        //double Om_w = 0.7;
        return (Om_w / E2(a))*g(a);
    }

    double w(double a){
        double ade=5.91331086e-4, s=3.20, w;
        w = ((4.0/3.0)/((pow(exp(a)/ade, s)+1.0)))-1.0;
        return w;
    }

    double g(double a){
        double ade=5.91331086e-4, s=3.20, gp;
        gp = (pow(exp(a),-4))*(pow((1.0+pow(exp(a)/ade,s))/(1.0+pow(1.0/ade,s)), 4.0/s));
        return gp;
    }

}

namespace EXP2{ //Modelo de energía oscura 2EXP doble exponencial

    int modo = 0;
    double z0 = 1e3;

    double g(double a);
    double w(double a);
    double Omega_w(double a);
    double Omega_m(double a);
    double E2(double a);
    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0);
    double dEdaE(double a);
    double E(double a);
    
    void modelo(){
        cout<<"\n \t\t Modelo 2EXP \n"<<endl;
        return;
    }
    
    double Om_m, Om_b, Om_w, Om_r, Om_k, w0;

    void Model_initialization(){

        cout<<"\n \t\t Modelo 2EXP \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << Om_m << endl;
        cout << "Barions: " << Om_b << endl;
        cout << "Radiation: " << Om_r << endl;
        cout << "Dark Energy: " << Om_w << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << Om_m+Om_r+Om_w << endl;        

        cout << endl;

        return;
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        //cambiar Exp_func por cosmologyFunc
        present_data(Om_m,Om_r,Om_w,Om_k, Om_b, w0);

        E_a = sqrt((Om_m / pow(exp(a),3)) + Om_w*g(a));
        eos = w(a);
        OmegaM = (Om_m / (pow(exp(a),3)*E_a*E_a));
        OmegaDE = (Om_w / (E_a*E_a))*g(a);
        ga = g(a);
        
        
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {

        m0_rho = 0.3103;
        r0_rho = 0;
        de0_rho = 0.6897;
        k0_rho = 0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        w0 = w(0);

        return;
    }

    double E2(double a){ //a = lna
        double E;
        E = (Om_m / pow(a,3)) + Om_w*g(a);
        return E;
    }

    double dEdaE(double a){
        double dE = -(3/2)*(Om_m*pow(exp(a),-3)+(1.0+w(a))*Om_w*g(a))/E2(a);
        return dE;
    }
    double E(double a){
        double E = sqrt(E2(a));
        return E;
    }

    double Omega_m(double a){ //a = lna
        //double Om_m = 0.3;
        return Om_m/(pow(a,3)*E2(a));
    }

    double Omega_w(double a){
        //double Om_w = 0.7;
        return (Om_w / E2(a))*g(a);
    }

    double w(double a){
        double w0=-1.0, wm=0.01, am=0.19,dm=0.043, w;
        double p1, p2, p3, p4;
        p1=1.0+exp(am/dm), p2=1.0-exp(-(exp(a)-1)/dm);
        p3=1.0+exp(-(exp(a)-am)/dm), p4=1.0-exp(1.0/dm); 
        w = w0 + (wm - w0)*((p1*p2)/(p3*p4));
        return w;
    }

    double g(double a){
        double gp, tf=log(1.0), w_integral;
        integrateTz2(a, tf, w, 1000, w_integral);
        gp = pow(exp(a),-3)*exp(3*w_integral);
        return gp;
    }

}

namespace devgf { //mode de melissa: Dark Energy from vector gauge fields

    int modo = 1;

    // Function prototype
    double w(double lna);
    void dynamic_eq(double lna, double z[], double dz_out[]);
    void present_data(double &m0_rho, double &r0_rho, double &de0_rho, double &k0_rho, double &b0_rho, double& w0);
    void evol_var(double lna, double &m_density, double &r_density, double &de_density, double &w, double &wp, double &b_density);
    double E2(double lna);
    void evol_Omega_m(double lna, double& m_density, double& m_density_p);
    double Omega_m(double lna);
    //double Omega_r(double lna);
    void evol_Omega_w(double lna, double& de_density, double& de_density_p);
    double Omega_w(double lna);
    double g(double lna);
    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]);
    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]);
    void g_function(double a, double& ga);
    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga);

    #define Number_of_splines 500;

    const int EOS_nsplines=Number_of_splines;

    const int intEOS_nsplines=Number_of_splines;

    void modelo(){
        cout<<"\n \t\t Modelo W dinámica-numérica \n"<<endl;
        return;
    }

    // Declaración de constantes
    const double lambda = 1.0;
    const double mu = 2.7;
    const double n = 0.001;
    const double xi = 1.0e-19;
    const double yi = 6.14e-15;
    const double zi = 1.52e-10;
    const double omega_bi = 0.7e-5;
    const double omega_ri = 0.99994;
    const double z0 = exp(18.84) - 1;

    // Equations of motion for the cosmological model
    void dynamic_eq(double lna, double z[], double dz_out[]) {
        double q, Om_m;

        // conventions:
        // z[0] --> x, z[1] --> y, z[2] --> z
        // z[3] --> rho_b, z[4] --> rho_r, Om --> rho_m

        Om_m = (1 - z[4] - z[3] - (z[0] * z[0]) - (z[1] * z[1])) / (1 + 2 * n);

        q = (1 + (z[0] * z[0]) - (3 * (z[1] * z[1])) * (1 + (4.0 / 3.0) * lambda * (z[2] * z[2])) + 4 * mu * Om_m * (z[2] * z[2]) - 2 * n * Om_m + z[4]) / 2.0;

        dz_out[0] = z[0] * (-n * (2.0 + q - 4 * mu * (z[2] * z[2])) * Om_m +
            2 * z[0] * z[2] * (-lambda * (z[1] * z[1]) + (1 + 2 * n) * mu * Om_m) +
            (z[0] * z[0]) * (1 - q + 8 * n * (1 + q) * Om_m)) / (-n * Om_m + (z[0] * z[0]) * (-1 + 8 * n * Om_m));

        dz_out[1] = z[1] * (q + 1) - 2 * lambda * z[1] * z[2] * (z[0] - z[2]);

        dz_out[2] = (z[0] - z[2]);

        dz_out[3] = 2 * z[3] * (q - 1. / 2.); //bariones

        dz_out[4] = 2 * (-1.0 + q) * z[4]; //radiación

        return;
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {
        double wp0;
        evol_var(0.0, m0_rho, r0_rho, de0_rho, w0, wp0, b0_rho);
        k0_rho = 0.0;
        return;
    }

    void evol_var(double lna, double& m_density, double& r_density, double& de_density, double& w, double& wp, double& b_density) {

        int flag = 1;
        int neqn = 5;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[5];
        dz = new double[5];
        //work = new double[100+21*neqn];

        double zp[3];
        double q, Om_m;

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        double t0 = log(1.0/(1.0+z0)); 
        double tf = lna;

        z[0] = xi;
        z[1] = yi;
        z[2] = zi;
        z[3] = omega_bi;
        z[4] = omega_ri;

        dynamic_eq(t0, z, dz);

        //ode ( dynamic_eq, neqn, z, t0, tf, r_error, abs_error, flag, work, iwork );
        flag = r8_rkf45(dynamic_eq, neqn, z, dz, &t0, tf, &r_error, abs_error, flag);

        m_density = ((1 - z[4] - z[3] - pow(z[0],2) - pow(z[1],2)) /(1 + 2*n)) + z[3];
        r_density = z[4];
        de_density = pow(z[0],2) + pow(z[1],2);
        b_density = z[3];
        w = (1.0/3.0) * (pow(z[0],2) - 4 * lambda * pow(z[1],2) * pow(z[2],2) - 3 * pow(z[1],2))/ (pow(z[0],2) + pow(z[1],2));
        //w = (2*q -1 - r_density)/(3*de_density);
        //=======================================Derivada de w===================================================
        Om_m = (1 - z[4] - z[3] - (z[0] * z[0]) - (z[1] * z[1])) / (1 + 2 * n);

        q = (1 + (z[0] * z[0]) - (3 * (z[1] * z[1])) * (1 + (4.0 / 3.0) * lambda * (z[2] * z[2])) + 4 * mu * Om_m * (z[2] * z[2]) - 2 * n * Om_m + z[4]) / 2.0;

        zp[0] = z[0] * (-n * (2.0 + q - 4 * mu * (z[2] * z[2])) * Om_m +
            2 * z[0] * z[2] * (-lambda * (z[1] * z[1]) + (1 + 2 * n) * mu * Om_m) +
            (z[0] * z[0]) * (1 - q + 8 * n * (1 + q) * Om_m)) / (-n * Om_m + (z[0] * z[0]) * (-1 + 8 * n * Om_m));
        
        zp[1] = z[1] * (q + 1) - 2 * lambda * z[1] * z[2] * (z[0] - z[2]);

        zp[2] = (z[0] - z[2]);

        wp = (8.0/3.0)*( z[1]*( -z[0]*z[1]*(1+lambda*pow(z[2],2))*zp[0]+lambda*pow(z[1],3)*z[2]*zp[2]+pow(z[0],2)*( (1+lambda*pow(z[1],2))*zp[1]+z[1]*z[2]*zp[2] ) ) )/pow((pow(z[0],2)+pow(z[1],2)),2);
        //=======================================================================================================
        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }

    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]){
    // # define NUMBER 500

        double w1, wp0, wpf, dx,dy,dw, db;
        double t0 = log(1.0/(1.0+z0));
        //double t0_0 = log(1.0/(1.0+(exp(18.8400001) - 1)));
        double tf=0.0;
        double dt=(tf-t0)/(EOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=1;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < EOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            //evol_var(t0 + (i)*dt,dx,dy,dw,w1,wp,db);
            vec_w[i] =  w(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        evol_var(t0,dx,dy,dw,w1,wp0,db);
        ybcbeg=wp0;
        evol_var(tf,dx,dy,dw,w1,wpf,db);
        ybcend=wpf;
        //================================================!

        vec_wpp = spline_cubic_set ( EOS_nsplines, vec_t, vec_w, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w


        delete [] vec_wpp;
    }

    // Equation of state for the model
    double w(double lna) {
        double matter, radiation, dark_e, w, wp, barions;
        evol_var(lna, matter, radiation, dark_e, w, wp, barions);
        return w;
    }

    double w_int(double lna){
        double w_int, tf=0.0;
        integrateTz2(lna, tf, w, 1000, w_int);
        return w_int;
    }

    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]){
    // # define NUMBER 500

        //double wp, w0, wf;
        double t0 = log(1.0/(1.0+z0));
        double tf=0.0;
        double dt=(tf-t0)/(intEOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=10;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < intEOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            vec_wInt[i] =  w_int(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        ybcbeg=w(t0);
        ybcend=w(tf);
        //================================================!

        vec_wIntpp = spline_cubic_set ( intEOS_nsplines, vec_t, vec_wInt, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w_int

        delete [] vec_wIntpp;
    }

    double E2(double lna){ //a = lna
        double m_rho, r_rho, de_rho, b_rho, w, wp;
        evol_var(0.0, m_rho, r_rho, de_rho, w, wp, b_rho);
        double E;
        E = (m_rho/(pow(exp(lna),3))) + de_rho*g(lna);
        return E;
    }

    void evol_Omega_m(double lna, double& m_density, double& m_density_p) {

        int flag = 1;
        int neqn = 5;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[5];
        dz = new double[5];
        //work = new double[100+21*neqn];

        double zp[5];
        double q, Om_m;

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        //int flag = 1;
        double t0 = log(1.0/(1.0+z0)); 
        double tf = lna;

        z[0] = xi;
        z[1] = yi;
        z[2] = zi;
        z[3] = omega_bi;
        z[4] = omega_ri;

        dynamic_eq(t0, z, dz);
        //ode ( dynamic_eq, neqn, z, t0, tf, r_error, abs_error, flag, work, iwork );
        flag = r8_rkf45(dynamic_eq, neqn, z, dz, &t0, tf, &r_error, abs_error, flag);

        m_density = ((1 - z[4] - z[3] - pow(z[0],2) - pow(z[1],2)) /(1 + 2*n)) + z[3];

        //=======================================Derivada de Omega_m===================================================
        Om_m = (1 - z[4] - z[3] - (z[0] * z[0]) - (z[1] * z[1])) / (1 + 2 * n);

        q = (1 + (z[0] * z[0]) - (3 * (z[1] * z[1])) * (1 + (4.0 / 3.0) * lambda * (z[2] * z[2])) + 4 * mu * Om_m * (z[2] * z[2]) - 2 * n * Om_m + z[4]) / 2.0;

        zp[0] = z[0] * (-n * (2.0 + q - 4 * mu * (z[2] * z[2])) * Om_m +
            2 * z[0] * z[2] * (-lambda * (z[1] * z[1]) + (1 + 2 * n) * mu * Om_m) +
            (z[0] * z[0]) * (1 - q + 8 * n * (1 + q) * Om_m)) / (-n * Om_m + (z[0] * z[0]) * (-1 + 8 * n * Om_m));
        
        zp[1] = z[1] * (q + 1) - 2 * lambda * z[1] * z[2] * (z[0] - z[2]);

        zp[2] = (z[0] - z[2]);

        zp[3] = 2 * z[3] * (q - 1. / 2.); //bariones

        zp[4] = 2 * (-1.0 + q) * z[4]; //radiación

        m_density_p = (2*z[0]*zp[0]+2*z[1]*zp[1]+zp[3]+zp[4])/(1+2*n);
        //=======================================================================================================
        
        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }

    double Omega_m(double lna){ //a = lna
        double m_rho, m_rho_p;
        evol_Omega_m(lna, m_rho, m_rho_p);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return m_rho;
    }

    void Om_spline_set(double vec_t[], double vec_Om[], double vec_Ompp[]){
    // # define NUMBER 500

        double om, omp0, ompf;
        double t0 = log(1.0/(1.0+z0));
        //double t0_0 = log(1.0/(1.0+(exp(18.85) - 1)));
        double tf=0.0;
        double dt=(tf-t0)/(EOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=1;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < EOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            //evol_var(vec_time[i],dx,dy,dw,w1,wp,db);
            vec_Om[i] =  Omega_m(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        evol_Omega_m(t0, om, omp0);
        ybcbeg=omp0;
        evol_Omega_m(t0, om, ompf);
        ybcend=ompf;
        //================================================!

        vec_Ompp = spline_cubic_set ( EOS_nsplines, vec_t, vec_Om, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w


        delete [] vec_Ompp;
    }

    double Omega_w(double lna){ //a = lna
        double de_rho, de_rho_p;
        evol_Omega_w(lna, de_rho, de_rho_p);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return de_rho;
    }

    void evol_Omega_w(double lna, double& de_density, double& de_density_p){

        int flag = 1;
        int neqn = 5;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[5];
        dz = new double[5];
        //work = new double[100+21*neqn];

        double zp[5];
        double q, Om_m;

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        //int flag = 1;
        double t0 = log(1.0/(1.0+z0)); 
        double tf = lna;

        z[0] = xi;
        z[1] = yi;
        z[2] = zi;
        z[3] = omega_bi;
        z[4] = omega_ri;

        dynamic_eq(t0, z, dz);
        //ode ( dynamic_eq, neqn, z, t0, tf, r_error, abs_error, flag, work, iwork );
        flag = r8_rkf45(dynamic_eq, neqn, z, dz, &t0, tf, &r_error, abs_error, flag);

        de_density = z[0]*z[0] + z[1]*z[1];

        //=======================================Derivada de Omega_de===================================================
        Om_m = (1 - z[4] - z[3] - (z[0] * z[0]) - (z[1] * z[1])) / (1 + 2 * n);

        q = (1 + (z[0] * z[0]) - (3 * (z[1] * z[1])) * (1 + (4.0 / 3.0) * lambda * (z[2] * z[2])) + 4 * mu * Om_m * (z[2] * z[2]) - 2 * n * Om_m + z[4]) / 2.0;

        zp[0] = z[0] * (-n * (2.0 + q - 4 * mu * (z[2] * z[2])) * Om_m +
            2 * z[0] * z[2] * (-lambda * (z[1] * z[1]) + (1 + 2 * n) * mu * Om_m) +
            (z[0] * z[0]) * (1 - q + 8 * n * (1 + q) * Om_m)) / (-n * Om_m + (z[0] * z[0]) * (-1 + 8 * n * Om_m));
        
        zp[1] = z[1] * (q + 1) - 2 * lambda * z[1] * z[2] * (z[0] - z[2]);

        zp[2] = (z[0] - z[2]);

        zp[3] = 2 * z[3] * (q - 1. / 2.); //bariones

        zp[4] = 2 * (-1.0 + q) * z[4]; //radiación

        de_density_p = 2*(z[0]*zp[0]+z[1]*zp[1]);
        //=======================================================================================================
        
        //delete [] work;
        delete [] z;
        delete [] dz;

    }

    void DE_spline_set(double vec_t[], double vec_Ode[], double vec_Odepp[]){
    // # define NUMBER 500

        double ode, odep0, odepf;
        double t0 = log(1.0/(1.0+z0));
        //double t0_0 = log(1.0/(1.0+(exp(18.85) - 1)));
        double tf=0.0;
        double dt=(tf-t0)/(EOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=1;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < EOS_nsplines; i++ ){
            vec_t[i] = exp(t0 + (i)*dt);
            //evol_var(vec_time[i],dx,dy,dw,w1,wp,db);
            vec_Ode[i] =  Omega_w(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        evol_Omega_w(t0, ode, odep0);
        ybcbeg=odep0;
        evol_Omega_w(t0, ode, odepf);
        ybcend=odepf;
        //================================================!

        vec_Odepp = spline_cubic_set ( EOS_nsplines, vec_t, vec_Ode, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w


        delete [] vec_Odepp;
    }

    double g(double lna){
        double gp, tf=0.0, eos_integral;
        integrateTz2(lna, tf, w, 1000, eos_integral);
        gp = exp(-3*lna)*exp(3*eos_integral);
        return gp;
    }

    //------------------------------Interpolación spline----------------------------
    //Asignar el número de splines
    int G_spline_numbe = Number_of_splines;

    double G_m0_density, G_r0_density, G_b0_density, G_k0_density, G_de0_density, w0;

    double* vec_w= new double[G_spline_numbe];
    double* vec_t= new double[G_spline_numbe];
    double* vec_wpp= new double[G_spline_numbe];
    double* vec_w_int= new double[G_spline_numbe];
    double* vec_w_int_pp= new double[G_spline_numbe];
    double* vec_om= new double[G_spline_numbe];
    double* vec_ompp= new double[G_spline_numbe];
    double* vec_ode= new double[G_spline_numbe];
    double* vec_odepp= new double[G_spline_numbe];

    //Función para inicializar el modelo
    void Model_initialization(){

        cout<<"\n \t\t Modelo W dinámica-numérica DEVGF \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density, G_b0_density, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << G_m0_density << endl;
        cout << "Barions: " << G_b0_density << endl;
        cout << "Radiation: " << G_r0_density << endl;
        cout << "Dark Energy: " << G_de0_density << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << G_m0_density+G_r0_density+G_de0_density << endl;        

        cout << endl;
        cout << "==========SPLINE INITIALIZATIONS================" << endl;
            
        //Inicializar el spline de w
        w_spline_set(vec_t, vec_w, vec_wpp);
            
        //Inicializar el spline de w_int_spline
        w_int_spline_set(vec_t, vec_w_int, vec_w_int_pp);

        //Om_spline_set(vec_t, vec_om, vec_ompp);

        //DE_spline_set(vec_t, vec_ode, vec_odepp);

        for (int i=0; i<=20; i++){
            cout<<vec_t[i]<< "  "<<vec_w[i]<< "  "<<vec_w_int[i]<<endl; 
        }

        cout << endl;

        cout << "================================================" << endl;
        

        //Liberar memoria dinámica
        // delete[] vec_w;
        // delete[] vec_w;
        // delete[] vec_wpp;
        // delete[] vec_w_int;
        // delete[] vec_w_int_pp;
        return;
    }

    void g_function(double a, double& ga){
        double y, yp, ypp;

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        //ga=exp(-3.0*a)*exp(-3.0*(y));
        ga=pow(exp(a),-3.0)*exp(3.0*(y));
        
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density, G_b0_density, w0);
        double Omega_M0 = G_m0_density;

        double Omega_DE0 = G_de0_density;
        double wp, wpp; 
        double y, yp, ypp;

        double m_rho, r_rho, de_rho, b_rho, w, wp1;
        evol_var(a, m_rho, r_rho, de_rho, w, wp1, b_rho);

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        ga=pow(exp(a),-3.0)*exp(3.0*(y));

        E_a = sqrt((Omega_M0/(pow(exp(a),3.0))) + (Omega_DE0 * ga));
        OmegaM = Omega_M0/(E_a*E_a*pow(exp(a),3));
        OmegaDE = (Omega_DE0/(E_a*E_a)) * ga;

        eos = spline_cubic_val(G_spline_numbe, vec_t, vec_w, vec_wpp, a, &wp, &wpp);        
        
    }

} //cierra namespace w-de-m

namespace nonabelian { //arXiv:2007.12964v2

    int modo = 1;

    // Function prototype
    double w(double lna);
    void dynamic_eq(double lna, double z[], double dz_out[]);
    void present_data(double &m0_rho, double &r0_rho, double &de0_rho, double &k0_rho, double& b0_rho, double& w0);
    void evol_var(double lna, double &m_density, double &r_density, double &de_density, double &w, double &wp);
    double E2(double lna);
    double Omega_m(double lna);
    double Omega_r(double lna);
    double Omega_w(double lna);
    double g(double lna);
    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]);
    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]);
    void g_function(double a, double& ga);
    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga);

    #define Number_of_splines 500;

    const int EOS_nsplines=Number_of_splines;

    const int intEOS_nsplines=Number_of_splines;

    void modelo(){
        cout<<"\n \t\t Modelo W dinámica-numérica non Abelian \n"<<endl;
        return;
    }

    // Declaración de constantes
    const double alpha = 1.0e9;
    const double xi = 3.0e-30;
    const double yi = 0.001;
    const double omega_i = 0.0; //6.0e-28 //3.0e-27 //5.0e-26;
    const double zi = 1.1e5;
    const double omega_ri = 0.99998;
    const double z0 = 2.18e8;

    // Equations of motion for the cosmological model
    void dynamic_eq(double lna, double z[], double dz_out[]) {
        
        double q;

        // conventions:
        // z[0] --> x, z[1] --> y, z[2] --> w
        // z[3] --> z, z[4] --> rho_r

        q = (1.0/2.0)*( 1+pow(z[0],2)*(1.0-12*alpha*pow(z[3],4))+pow(z[1],2)-pow(z[2],2)+z[4] );

        dz_out[0] = z[0]*q-( 2*(pow(z[1],2)/z[3])+8*alpha*pow(z[0],2)*pow(z[3],3)+z[0]*(1-12*alpha*pow(z[3],4))+(pow(z[2],2)/z[3]) )/(1+4*alpha*pow(z[3],4));

        dz_out[1] = z[1]*( 2*(z[0]/z[3])+q-1 );

        dz_out[2] = z[2]*( (z[0]/z[3])+q );

        dz_out[3] = z[0]-z[3]; 

        dz_out[4] = 2*z[4]*( q-1 );

        return;
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {
        double wp0;
        evol_var(0.0, m0_rho, r0_rho, de0_rho, w0, wp0);
        k0_rho = 0.0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        return;
    }

    void evol_var(double lna, double& m_density, double& r_density, double& de_density, double& w, double& wp) {

        int flag = 1;
        int neqn = 5;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[5];
        dz = new double[5];
        //work = new double[100+21*neqn];

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        double num, nump, den, denp; //para wp=dw/dlna

        double t0 = log(1.0/(1.0+z0)); 
        double tf = lna;

        z[0] = xi;
        z[1] = yi;
        z[2] = omega_i;
        z[3] = zi;
        z[4] = omega_ri;

        dynamic_eq(t0, z, dz);

        //ode ( dynamic_eq, neqn, z, t0, tf, r_error, abs_error, flag, work, iwork );
        flag = r8_rkf45(dynamic_eq, neqn, z, dz, &t0, tf, &r_error, abs_error, flag);

        m_density = 1-pow(z[0],2)*(1+4*alpha*pow(z[3],4))-pow(z[1],2)-pow(z[2],2)-z[4];
        r_density = z[4];
        de_density = pow(z[0],2)*(1+4*alpha*pow(z[3],4))+pow(z[1],2)+pow(z[2],2);
        
        w = (1.0/3.0)*( pow(z[0],2)*(1-12*alpha*pow(z[3],4))+pow(z[1],2)-pow(z[2],2) )/( pow(z[0],2)*(1+4*alpha*pow(z[3],4))+pow(z[1],2)+pow(z[2],2) );

        //=======================================Derivada de w===================================================

        num = pow(z[0],2)*(1-12*alpha*pow(z[3],4))+pow(z[1],2)-pow(z[2],2);

        nump = 2*z[0]*dz[0]*(1-12*alpha*pow(z[3],4))-48*pow(z[0],2)*alpha*pow(z[3],3)*dz[3]+2*z[1]*dz[1]-2*z[2]*dz[2];

        den = pow(z[0],2)*(1+4*alpha*pow(z[3],4))+pow(z[1],2)+pow(z[2],2);

        denp = 2*z[0]*dz[0]*(1+4*alpha*pow(z[3],4))+16*pow(z[0],2)*alpha*pow(z[3],3)*dz[3]+2*z[1]*dz[1]+2*z[2]*dz[2];

        wp = (1/3)*(nump*den - num*denp)/(pow(den,2));
        //=======================================================================================================
        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }

    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]){
    // # define NUMBER 500

        double w1, wp0, wpf, dx,dy,dw;
        double t0 = log(1.0/(1.0+z0));
        //double t0_0 = log(1.0/(1.0+(exp(18.8400001) - 1)));
        double tf=0.0;
        double dt=(tf-t0)/(EOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=1;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < EOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            //evol_var(t0 + (i)*dt,dx,dy,dw,w1,wp,db);
            vec_w[i] =  w(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        evol_var(t0,dx,dy,dw,w1,wp0);
        ybcbeg=wp0;
        evol_var(tf,dx,dy,dw,w1,wpf);
        ybcend=wpf;
        //================================================!

        vec_wpp = spline_cubic_set ( EOS_nsplines, vec_t, vec_w, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w


        delete [] vec_wpp;
    }

    // Equation of state for the model
    double w(double lna) {
        double matter, radiation, dark_e, w, wp;
        evol_var(lna, matter, radiation, dark_e, w, wp);
        return w;
    }

    double w_int(double lna){
        double w_int, tf=0.0;
        integrateTz2(lna, tf, w, 1000, w_int);
        return w_int;
    }

    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]){
    // # define NUMBER 500

        //double wp, w0, wf;
        double t0 = log(1.0/(1.0+z0));
        double tf=0.0;
        double dt=(tf-t0)/(intEOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=10;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < intEOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            vec_wInt[i] =  w_int(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        ybcbeg=w(t0);
        ybcend=w(tf);
        //================================================!

        vec_wIntpp = spline_cubic_set ( intEOS_nsplines, vec_t, vec_wInt, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w_int

        delete [] vec_wIntpp;
    }

    double E2(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(0.0, m_rho, r_rho, de_rho, w, wp);
        double E;
        E = (m_rho/(pow(exp(lna),3))) + de_rho*g(lna);
        return E;
    }

    double Omega_m(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return m_rho;
    }

    double Omega_r(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return r_rho;
    }

    double Omega_w(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return de_rho;
    }

    double g(double lna){
        double gp, tf=0.0, eos_integral;
        integrateTz2(lna, tf, w, 1000, eos_integral);
        gp = exp(-3*lna)*exp(3*eos_integral);
        return gp;
    }

    //------------------------------Interpolación spline----------------------------
    //Asignar el número de splines
    int G_spline_numbe = Number_of_splines;

    double G_m0_density, G_r0_density, G_k0_density, G_de0_density, G_b0_density, w0;

    double* vec_w= new double[G_spline_numbe];
    double* vec_t= new double[G_spline_numbe];
    double* vec_wpp= new double[G_spline_numbe];
    double* vec_w_int= new double[G_spline_numbe];
    double* vec_w_int_pp= new double[G_spline_numbe];
    double* vec_om= new double[G_spline_numbe];
    double* vec_ompp= new double[G_spline_numbe];
    double* vec_ode= new double[G_spline_numbe];
    double* vec_odepp= new double[G_spline_numbe];

    //Función para inicializar el modelo
    void Model_initialization(){

        cout<<"\n \t\t Modelo W dinámica-numérica non Abelian \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density, G_b0_density, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << G_m0_density << endl;
        cout << "Barions: " << G_b0_density << endl;
        cout << "Radiation: " << G_r0_density << endl;
        cout << "Dark Energy: " << G_de0_density << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << G_m0_density+G_r0_density+G_de0_density << endl;        

        cout << endl;
        cout << "==========SPLINE INITIALIZATIONS================" << endl;
            
        //Inicializar el spline de w
        w_spline_set(vec_t, vec_w, vec_wpp);
            
        //Inicializar el spline de w_int_spline
        w_int_spline_set(vec_t, vec_w_int, vec_w_int_pp);

        for (int i=0; i<=20; i++){
            cout<<vec_t[i]<< "  "<<vec_w[i]<< "  "<<vec_w_int[i]<<endl; 
        }

        cout << endl;

        cout << "================================================" << endl;
        

        //Liberar memoria dinámica
        // delete[] vec_w;
        // delete[] vec_w;
        // delete[] vec_wpp;
        // delete[] vec_w_int;
        // delete[] vec_w_int_pp;
        return;
    }

    void g_function(double a, double& ga){
        double y, yp, ypp;

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        //ga=exp(-3.0*a)*exp(-3.0*(y));
        ga=pow(exp(a),-3.0)*exp(3.0*(y));
        
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density, G_b0_density, w0);
        double Omega_M0 = G_m0_density;

        double Omega_DE0 = G_de0_density;
        double wp, wpp;
        double y, yp, ypp;

        double m_rho, r_rho, de_rho, w, wp1;
        evol_var(a, m_rho, r_rho, de_rho, w, wp1);

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        ga=pow(exp(a),-3.0)*exp(3.0*(y));

        E_a = sqrt((Omega_M0/(pow(exp(a),3.0))) + (Omega_DE0 * ga));
        OmegaM = Omega_M0/(E_a*E_a*pow(exp(a),3));
        OmegaDE = (Omega_DE0/(E_a*E_a)) * ga;

        eos = spline_cubic_val(G_spline_numbe, vec_t, vec_w, vec_wpp, a, &wp, &wpp);        
        
    }

} //cierra non_abelian

namespace twoform { //arXiv:1902.05846v1

    int modo = 1;

    // Function prototype
    double w(double lna);
    void dynamic_eq(double lna, double z[], double dz_out[]);
    void present_data(double &m0_rho, double &r0_rho, double &de0_rho, double &k0_rho, double& b0_rho, double& w0);
    void evol_var(double lna, double &m_density, double &r_density, double &de_density, double &w, double &wp);
    double E2(double lna);
    double Omega_m(double lna);
    double Omega_r(double lna);
    double Omega_w(double lna);
    double g(double lna);
    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]);
    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]);
    void g_function(double a, double& ga);
    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga);

    #define Number_of_splines 500;

    const int EOS_nsplines=Number_of_splines;

    const int intEOS_nsplines=Number_of_splines;

    void modelo(){
        cout<<"\n \t\t Modelo W dinámica-numérica 2-form \n"<<endl;
        return;
    }

    // Declaración de constantes
    const double lambda = 2.0;
    const double mu = 50.0;
    const double x1i = 1.0e-13;
    const double x2i = 1.0e-14;
    const double sigma = 0.0; 
    const double omega_bi = 1.0e-10;
    const double omega_ri = 0.999961;
    const double z0 = 8.3e7;

    // Equations of motion for the cosmological model
    void dynamic_eq(double lna, double z[], double dz_out[]) {

        // conventions:
        // z[0] --> x1, z[1] --> x2, z[2] --> sigma
        // z[3] --> rho_b, z[4] --> rho_r

        dz_out[0] = (3.0/2.0)*z[0]*(pow(z[0],2)-pow(z[1],2)+pow(z[2],2)-1.0-(1.0/3.0)*z[3]+(1.0/3.0)*z[4])+(sqrt(6)/2.0)*(lambda*pow(z[1],2)-mu*z[3]);

        dz_out[1] = (1.0/2.0)*z[1]*(3.0*pow(z[0],2)-3.0*pow(z[1],2)+3.0*pow(z[2],2)+3.0-sqrt(6.0)*lambda*z[0]-z[3]+z[4]);

        dz_out[2] = (1.0/2.0)*z[2]*(3.0*pow(z[0],2)-3.0*pow(z[1],2)+3.0*pow(z[2],2)-3.0-z[3]+z[4])-2.0*z[3];

        dz_out[3] = z[3]*(3.0*pow(z[0],2)-3.0*pow(z[1],2)+3.0*pow(z[2],2)+4*z[2]+1.0+sqrt(6.0)*mu*z[0]-z[3]+z[4]); 

        dz_out[4] = z[4]*(3.0*pow(z[0],2)-3.0*pow(z[1],2)+3.0*pow(z[2],2)-1.0-z[3]+z[4]);

        return;
    }

    void present_data(double& m0_rho, double& r0_rho, double& de0_rho, double& k0_rho, double& b0_rho, double& w0) {
        double wp0;
        evol_var(0.0, m0_rho, r0_rho, de0_rho, w0, wp0);
        k0_rho = 0.0;
        b0_rho = (0.02234/0.14205)*m0_rho;
        return;
    }

    void evol_var(double lna, double& m_density, double& r_density, double& de_density, double& w, double& wp) {

        int flag = 1;
        int neqn = 5;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[5];
        dz = new double[5];
        //work = new double[100+21*neqn];

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        double num, nump, den, denp; //para wp=dw/dlna

        double t0 = log(1.0/(1.0+z0)); 
        double tf = lna;

        z[0] = x1i;
        z[1] = x2i;
        z[2] = sigma;
        z[3] = omega_bi;
        z[4] = omega_ri;

        dynamic_eq(t0, z, dz);

        //ode ( dynamic_eq, neqn, z, t0, tf, r_error, abs_error, flag, work, iwork );
        flag = r8_rkf45(dynamic_eq, neqn, z, dz, &t0, tf, &r_error, abs_error, flag);

        m_density = 1-pow(z[0],2)-pow(z[1],2)-pow(z[2],2)-z[3]-z[4];
        r_density = z[4];
        de_density = 1-m_density-r_density;
        
        w = (3.0*(pow(z[0],2)-pow(z[1],2)+pow(z[2],2))-z[3])/(3.0*(pow(z[0],2)+pow(z[1],2)+pow(z[2],2)+z[3]));

        //=======================================Derivada de w===================================================

        num = 3.0*(pow(z[0],2)-pow(z[1],2)+pow(z[2],2))-z[3];

        nump = 6.0*(z[0]*dz[0]-z[1]*dz[1]+z[2]*dz[2])-dz[3];

        den = 3.0*(pow(z[0],2)+pow(z[1],2)+pow(z[2],2)+z[3]);

        denp = 6.0*(z[0]*dz[0]+z[1]*dz[1]+z[2]*dz[2]+dz[3]);

        wp = (nump*den - num*denp)/(pow(den,2));
        //=======================================================================================================
        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }

    void w_spline_set(double vec_t[], double vec_w[], double vec_wpp[]){
    // # define NUMBER 500

        double w1, wp0, wpf, dx,dy,dw;
        double t0 = log(1.0/(1.0+z0));
        //double t0_0 = log(1.0/(1.0+(exp(18.8400001) - 1)));
        double tf=0.0;
        double dt=(tf-t0)/(EOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=1;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < EOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            //evol_var(t0 + (i)*dt,dx,dy,dw,w1,wp,db);
            vec_w[i] =  w(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        evol_var(t0,dx,dy,dw,w1,wp0);
        ybcbeg=wp0;
        evol_var(tf,dx,dy,dw,w1,wpf);
        ybcend=wpf;
        //================================================!

        vec_wpp = spline_cubic_set ( EOS_nsplines, vec_t, vec_w, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w


        delete [] vec_wpp;
    }

    // Equation of state for the model
    double w(double lna) {
        double matter, radiation, dark_e, w, wp;
        evol_var(lna, matter, radiation, dark_e, w, wp);
        return w;
    }

    double w_int(double lna){
        double w_int, tf=0.0;
        integrateTz2(lna, tf, w, 1000, w_int);
        return w_int;
    }

    void w_int_spline_set(double vec_t[], double vec_wInt[], double vec_wIntpp[]){
    // # define NUMBER 500

        //double wp, w0, wf;
        double t0 = log(1.0/(1.0+z0));
        double tf=0.0;
        double dt=(tf-t0)/(intEOS_nsplines-1);

        int i;
        int ibcbeg=3;
        int ibcend=10;
        double ybcbeg;
        double ybcend;

        cout<< setprecision(17);
        //Inicializar vec_w
        for ( i = 0; i < intEOS_nsplines; i++ ){
            vec_t[i] = (t0 + (i)*dt);
            vec_wInt[i] =  w_int(vec_t[i]);

        }

        //Boundary conditions at the edges of the interval
        //================================================!
        ybcbeg=w(t0);
        ybcend=w(tf);
        //================================================!

        vec_wIntpp = spline_cubic_set ( intEOS_nsplines, vec_t, vec_wInt, ibcbeg, ybcbeg, ibcend, ybcend ); //Segunda deriada de w_int

        delete [] vec_wIntpp;
    }

    double E2(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(0.0, m_rho, r_rho, de_rho, w, wp);
        double E;
        E = (m_rho/(pow(exp(lna),3))) + de_rho*g(lna);
        return E;
    }

    double Omega_m(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return m_rho;
    }

    double Omega_r(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return r_rho;
    }

    double Omega_w(double lna){ //a = lna
        double m_rho, r_rho, de_rho, w, wp;
        evol_var(lna, m_rho, r_rho, de_rho, w, wp);
        //return (m0_rho) / (pow((a),3)*E2(a));
        return de_rho;
    }

    double g(double lna){
        double gp, tf=0.0, eos_integral;
        integrateTz2(lna, tf, w, 1000, eos_integral);
        gp = exp(-3*lna)*exp(3*eos_integral);
        return gp;
    }

    //------------------------------Interpolación spline----------------------------
    //Asignar el número de splines
    int G_spline_numbe = Number_of_splines;

    double G_m0_density, G_r0_density, G_k0_density, G_de0_density, G_b0_density, w0;

    double* vec_w= new double[G_spline_numbe];
    double* vec_t= new double[G_spline_numbe];
    double* vec_wpp= new double[G_spline_numbe];
    double* vec_w_int= new double[G_spline_numbe];
    double* vec_w_int_pp= new double[G_spline_numbe];
    double* vec_om= new double[G_spline_numbe];
    double* vec_ompp= new double[G_spline_numbe];
    double* vec_ode= new double[G_spline_numbe];
    double* vec_odepp= new double[G_spline_numbe];

    //Función para inicializar el modelo
    void Model_initialization(){

        cout<<"\n \t\t Modelo W dinámica-numérica 2-form \n"<<endl;

        cout << setprecision(15);

        //Llamado a función externa para obtener los datos iniciales
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density,G_b0_density, w0);

        //Imprimir densidades iniciales
        cout << "Densities today...(Intial data)" << endl;
        cout << "Matter: " << G_m0_density << endl;
        cout << "Barions: " << G_b0_density << endl;
        cout << "Radiation: " << G_r0_density << endl;
        cout << "Dark Energy: " << G_de0_density << endl;
        cout << "State of equation: " << w0 << endl;
        cout << "Total: " << G_m0_density+G_r0_density+G_de0_density << endl;        

        cout << endl;
        cout << "==========SPLINE INITIALIZATIONS================" << endl;
            
        //Inicializar el spline de w
        w_spline_set(vec_t, vec_w, vec_wpp);
            
        //Inicializar el spline de w_int_spline
        w_int_spline_set(vec_t, vec_w_int, vec_w_int_pp);

        for (int i=0; i<=20; i++){
            cout<<vec_t[i]<< "  "<<vec_w[i]<< "  "<<vec_w_int[i]<<endl; 
        }

        cout << endl;

        cout << "================================================" << endl;
        

        //Liberar memoria dinámica
        // delete[] vec_w;
        // delete[] vec_w;
        // delete[] vec_wpp;
        // delete[] vec_w_int;
        // delete[] vec_w_int_pp;
        return;
    }

    void g_function(double a, double& ga){
        double y, yp, ypp;

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        //ga=exp(-3.0*a)*exp(-3.0*(y));
        ga=pow(exp(a),-3.0)*exp(3.0*(y));
        
    }

    void cosmologicFuntion(double a, double &E_a, double& eos, double& OmegaM, double& OmegaDE, double& ga) {
        
        present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density,G_b0_density, w0);
        double Omega_M0 = G_m0_density;

        double Omega_DE0 = G_de0_density;
        double wp, wpp;
        double y, yp, ypp;

        double m_rho, r_rho, de_rho, w, wp1;
        evol_var(a, m_rho, r_rho, de_rho, w, wp1);

        y = spline_cubic_val(G_spline_numbe, vec_t, vec_w_int, vec_w_int_pp, a, &yp, &ypp);
        ga=pow(exp(a),-3.0)*exp(3.0*(y));

        E_a = sqrt((Omega_M0/(pow(exp(a),3.0))) + (Omega_DE0 * ga));
        OmegaM = Omega_M0/(E_a*E_a*pow(exp(a),3));
        OmegaDE = (Omega_DE0/(E_a*E_a)) * ga;

        eos = spline_cubic_val(G_spline_numbe, vec_t, vec_w, vec_wpp, a, &wp, &wpp);        
        
    }

} //cierra two_form
