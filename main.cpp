#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cstdlib>
//#include <thread>

#include "lib/scm.cpp"

using namespace scm; //Para poder usar las ecuaciones del SCM, scm.cpp

int main() {
    
    double M1014=1e14;

    double number = 1e3;
    double t0 = log(1e-7);
    //double tf=log(1.0);
    double dmin = 1e-7, dmax = 1e-3;
    ddnli=dmin/t0;
    //double dt=(tf-t0)/number;

    Model_initialization(); //Inicialización del Modelo utilizado en el cálculo.

    cout << setprecision(10);

    if (cosm==1){

    cout<<" Calculando datos cosmológicos ... "<<endl;
    ofstream cosmology;
    cosmology.open("datos/cosmology.dat");
    if (cosmology.fail()){
        cout<<" El archivo no se pudo crear. "<<endl;
    }
    t0 = log(1/(1+z0));
    cosmology << setprecision(17);
    for (int i=1; i<=5000; i++){
        double E_a, EoS, Omega_M, Omega_DE, g;
        cosmologicFuntion(t0, E_a, EoS, Omega_M, Omega_DE, g);
        cosmology << scientific << 1/exp(t0)
        << setw(1) <<" " << E_a 
        << setw(1) <<" " << Omega_M 
        << setw(1) <<" " << Omega_DE 
        << setw(1) <<" " << w(t0) 
        << setw(1) <<" " << g 
        << endl;
        //cout << 1/exp(t0) << setw(1) <<"\t  " << Omega_m(t0) << setw(1) <<"\t  " << Omega_w(t0) << setw(1) <<"\t  " << w(t0) << setw(1) <<"\t  " << g(t0) << endl;
        t0+=0.005;
        if (t0>=5.0) break;
    }
    cout << "Se escribió el archivo cosmology.dat" <<endl;
    cosmology.close();
    //exit(0); //Salir despues de escribir cosmology.dat

    }
    //----------------------------------------------------------------------------------
    cout << " Buscando condiciones iniciales ... " <<endl;
    for(int j=1; j<=1e9; j++){

        dnli = (2.0e-7)*exp(log(dmax/dmin)*((j-1)*1e-3/(number-1)));
        //dnli = (4.0828638e-7)*exp(log(dmax/dmin)*((j-1)*1e-3/(number-1))); //2EXP
        //dnli = (3.7081618e-7)*exp(log(dmax/dmin)*((j-1)*1e-3/(number-1))); //W_de_M
        //dnli = (3.630665049e-7)*exp(log(dmax/dmin)*((j-1)*1e-3/(number-1))); //nonAbelian
        //dnli = (3.619035065e-7)*exp(log(dmax/dmin)*((j-1)*1e-3/(number-1))); //2-form
        t0=log(1.0);
        ci(t0, ddnli, dnli, dnl, d0);
        if(dnl > 1e8) break;
    }
    //-----------------------------------------------------------------------------------
    ofstream datos;
    datos.open("datos/datos-scm-parameters.dat");
    if (datos.fail()){
        cout<<" El archivo no se pudo crear. "<<endl;
    }

    ofstream datos2;
    datos2.open("datos/datos-halos.dat");
    if (datos2.fail()){
        cout<<" El archivo no se pudo crear. "<<endl;
    }

    datos << setprecision(17) <<fixed;
    datos2 << setprecision(17) <<fixed;

    dnli0 = dnli; //Condición inicial de colapso hoy, para calcular D.
    t0 = log(1.0);
    if(ddelta == 1){
        evol_delta(t0, ddnli, dnli, dnl, dl);
    }

    cout << "     z" 
    << setw(3) << "               dc" 
    << setw(3) << "                 D/a" 
    << setw(3) << "                 dnl_ta" 
    << setw(3) << "               y_vir"
    << setw(3) << "             Delta_vir"
    /* << setw(3) << "        dndM(M, z)"
    << setw(3) << "          dV/dz" 
    << setw(3) << "             N_bin1(z)"
    << setw(3) << "           N_bin2(z)"
    << setw(3) << "           N_bin3(z)" 
    << setw(3) << "          n_inf1(z)"
    << setw(3) << "         n_inf2(z)"
    << setw(3) << "         n_inf3(z)"  */ 
    << endl;

    for(t0=log(1.0); t0>=log(1.0/(11.4)); t0-=0.01){
        
        evol_D(t0, D);
        scm_parameters (t0, ddnli, dnli, dnl, dc, dnl_taA, y_vir, Delta_vir);
        /* N_bin(t0, n_bin, n_bin2, n_bin3);
        N_inf(t0, n_inf, n_inf2, n_inf3); */
        
        datos << scientific << (1/exp(t0))-1 
        << setw(3) << "  " << dnli 
        << setw(3) << "  " << dc 
        << setw(3) << "  " << D/exp(t0)
        << setw(3) << "  " << dnl_taA
        << setw(3) << "  " << y_vir
        << setw(3) << "  " << Delta_vir 
        << endl;

        /* datos2 << scientific << (1/exp(t0))-1
        << setw(3) << "  " << dndM(M1014) 
        << setw(3) << "  " << evol_vol(t0) 
        << setw(3) << "  " << n_bin
        << setw(3) << "  " << n_bin2
        << setw(3) << "  " << n_bin3 
        << setw(3) << "  " << n_inf
        << setw(3) << "  " << n_inf2
        << setw(3) << "  " << n_inf3  
        << endl; */

        cout << scientific << (1/exp(t0))-1 
        << setw(3) << "  " << dc 
        << setw(3) << "  " << D/exp(t0)
        << setw(3) << "  " << dnl_taA
        << setw(3) << "  " << y_vir
        << setw(3) << "  " << Delta_vir
        /* << setw(3) << "  " << dndM(M1014)
        << setw(3) << "  " << evol_vol(t0) 
        << setw(3) << "  " << n_bin
        << setw(3) << "  " << n_bin2
        << setw(3) << "  " << n_bin3 
        << setw(3) << "  " << n_inf 
        << setw(3) << "  " << n_inf2
        << setw(3) << "  " << n_inf3 */
        << endl;
    
    }
    
    datos.close();
    datos2.close();

    //-------------Tiempo de ejecución del programa------------------------------//
    cout<<"\n"<<clock()/CLOCKS_PER_SEC<<" segundos"<<endl;
    //---------------------------------------------------------------------------//

    return 0;
}