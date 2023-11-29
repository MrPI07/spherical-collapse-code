#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include<fstream>
#include <ctime>
# include <cfloat>

//#include "rkf45.hpp"
//#include "DE-models.cpp"

#include "../modelSelect.cpp"

namespace scm {

    void fnl(double lna, double dnl[], double ddnl_out[]);
    void evol_fnl(double lna, double& dnli, double& ddnli, double& dnl);
    void fl(double lna, double d[], double dd_out[]);
    void evol_fl(double lna, double& dnli, double& ddnli, double& d);
    void scm_parameters (double t0, double ddnli, double& dnli, double& dnl, double& dc, double& dnl_ta, double& y_vir, double& Delta_vir);
    //void evol_vol(double lna, double& vol);

    //=========================Forma no lineal=======================================
    void fnl(double lna, double dnl[], double ddnl_out[]){

        double E_a, Eos, Omega_M, Omega_DE, ga;
        cosmologicFuntion(lna, E_a, Eos, Omega_M, Omega_DE, ga);

        ddnl_out[0] = dnl[1];
        ddnl_out[1] = -((1.0/2.0)*(1.0-3*Eos*Omega_DE))*dnl[1] + (4.0/3.0)*(pow(dnl[1],2)/(1.0 + dnl[0])) + 
                        (3.0/2.0)*(1+dnl[0])*dnl[0]*(Omega_M);

    }

    void evol_fnl(double lna, double& dnli, double& ddnli, double& dnl){

        int flag = 1;
        //int neqn = 2;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[2];
        dz = new double[2];
        //work = new double[100+21*neqn];

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = 1e-9;//sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        double r_error2 = 7.0e-9;

        //int flag = 1;
        double t0 = log(1.0e-7); 
        double tf = lna;

        z[0] = dnli;
        z[1] = ddnli;//z[0]/t0;

        fnl(t0, z, dz);
        if (modo ==0){
            flag = r8_rkf45(fnl, 2, z, dz, &t0, tf, &r_error, abs_error, flag);
            dnl = z[0];
        }
        else if (modo==1){
            flag = r4_rkf45(fnl, 2, z, dz, &t0, tf, &r_error2, r_error2, flag);
            //ode ( fnl, neqn, z, t0, tf, r_error, 5e-9, flag, work, iwork );
            dnl = z[0];
            // cout<<"Modo numérico en formación"<<endl;
            // exit(1);
        }
        //dnli = z[1];

        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }
    //==========================================================================

    //========================Forma lineal======================================
    void fl(double lna, double d[], double dd_out[]){
        
        double E_a, Eos, Omega_M, Omega_DE, ga;
        cosmologicFuntion(lna, E_a, Eos, Omega_M, Omega_DE, ga);

        dd_out[0] = d[1];
        dd_out[1] = -((1.0/2.0)*(1.0-3*Eos*Omega_DE))*d[1] + (3.0/2.0)*d[0]*(Omega_M);

    }

    void evol_fl(double lna, double& dnli, double& ddnli, double& dl){

        int flag = 1;
        //int neqn = 2;
        //int iwork[5];

        double* z;
        double* dz;
        //double *work;

        z = new double[2];
        dz = new double[2];
        //work = new double[100+21*neqn];

        double r_error = sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root
        double abs_error = 1e-9;//sqrt(numeric_limits<double>::epsilon()); // Use sqrt to calculate the square root

        double r_error2 = 7.0e-9;

        //int flag = 1;
        double t0 = log(1.0e-7); 
        double tf = lna;

        z[0] = dnli;
        z[1] = ddnli;//z[0]/t0;

        fl(t0, z, dz);
        if (modo==0){
            flag = r8_rkf45(fl, 2, z, dz, &t0, tf, &r_error, abs_error, flag);
            dl = z[0];
        }
        else if(modo==1){
            flag = r4_rkf45(fl, 2, z, dz, &t0, tf, &r_error2, r_error2, flag);
            //ode ( fl, neqn, z, t0, tf, r_error, 5e-9, flag, work, iwork );
            dl = z[0];
        }
        //dnli = z[1];

        //delete [] work;
        delete [] z;
        delete [] dz;

        return;
    }
    //==========================================================================
    //Buscar las condiciones iniciales
    void ci(double t0, double ddnli, double& dnli, double& dnl, double& dl0){
        evol_fnl(t0, dnli, ddnli, dnl);
        //cout << exp(t0) <<"  "<< dnli <<"\t"<< dnl << "\t" <<clock()/CLOCKS_PER_SEC<<" segundos "<<endl;
        if (dnl >= 1e8){
                evol_fl(t0, dnli, ddnli, dl0);
                cout << "\nCondiciones iniciales encontradas: (a_c, delta_i, delta_nl(a_c), delta_l(a_c) ) \n" 
                    << exp(t0) <<setw(1)<<"\t"<< dnli <<setw(1)<<"\t"<< dnl <<setw(1)<<"\t"<< dl0 <<endl;
                cout <<endl;
                //-------------Tiempo de ejecución del programa------------------------------//
                cout<<clock()/CLOCKS_PER_SEC<<" segundos"<<endl;
                //---------------------------------------------------------------------------//
        }
    }

    //Evolución de la sobredensidad no lineal y lineal, colapso hoy en día
    void evol_delta(double t0, double ddnli, double dnli, double& dnl, double& dl){

        ofstream evol_delta;
        evol_delta.open("datos/datos-delta.dat");
        if (evol_delta.fail()){
            cout<<" El archivo no se pudo crear. "<<endl;
        }
        evol_delta << setprecision(17);
        cout <<" \nCalculando y escribiendo datos-delta.dat ..." <<endl;
        for(int i = 1; i <= 1e5; i++){
            //t0-=dt;
            evol_fnl(t0, dnli, ddnli, dnl);
            evol_fl(t0, dnli, ddnli, dl);
            evol_delta << scientific << exp(t0) <<setw(3)<<"  "<< dnli <<setw(3)<<"  "<< dnl <<setw(3)<<"  "<< dl <<endl ;
            //cout << scientific << exp(t0) <<setw(2)<<"  "<< dnli <<setw(2)<<"  "<< dnl <<setw(2)<<"  "<< dl <<endl ;
            t0-=((log(1-0)-log(1e-7))/1e5);
        }
        cout << "Se escribió el archivo datos-delta.dat" <<endl;
        evol_delta.close();
        //-------------------------------------------------------------------//
        //-------------Tiempo de ejecución del programa------------------------------//
        cout<<clock()/CLOCKS_PER_SEC<<" segundos"<<endl;
        cout<<"\n";
        //---------------------------------------------------------------------------//
    }

    //Evolución del factor de crecimiento
    void evol_D(double t0, double& D){
    //double evol_D(double t0){
        double d0, dnl;
        //Calculo de la función de crecimiento D+(a)=d(a)/d(a=1)
        scm_parameters (log(1.0), ddnli,  dnli0,  dnl,  d0, dnl_taA, y_vir, Delta_vir);
        evol_fl(t0, dnli0, ddnli, D);
        D = D/d0;
        //return D;
    }

    //Evolución de la sobredensidad crítica
    void scm_parameters (double t0, double ddnli, double& dnli, double& dnl, double& dc, double& dnl_taA, double& y_vir, double& Delta_vir){

        /*
        Calculo de la sobredensidad no lineal dnl.
        Calculo de la sobredensidad crítica, dc(ac) -> dnl >= 1e8 para dnli dado.
        Calculo de la densidad en giro dnl_ta, t_c = 2t_ta.
        Calculo del radio virializado y_vir.
        Calculo de la sobredensidad virializada Delta_vir.
        */

        double var1, var2, varAux;
        double dmin = 1e-7, dmax = 1e-3;

        double y, varAuxR, vary1, vary2, a_ta;
        
        double number = 1e3;

        double Om_m0,Om_r0,Om_w0,Om_k0, Om_b0, w0;
        double E_a_ta, Eos_ta, Omega_M_ta, Omega_DE_ta, ga_ta, n2;
        double E_a_c, Eos_c, Omega_M_c, Omega_DE_c, ga_c, n1;
        //double y_vir;

        present_data(Om_m0,Om_r0,Om_w0,Om_k0, Om_b0, w0);

        var1 = dnli;
        var2 = dnli;
                    
        for(int j=1; j<=1e20; j++){

            evol_fnl(t0, dnli, ddnli, dnl);

            if((dnl >= 1e8)){
                evol_fl(t0, dnli, ddnli, dc);

                //------Implementar rutina para calcular a_ta--------
                y = log((dnl + 1)/pow(exp(t0),3));
                vary1 = y;
                vary2 = y;
                //t_ta = t0+log(pow(2.0,-2.0/3.0));
                for(a_ta = t_ta; a_ta>=log(1.0e-4); a_ta-=0.00005){
                    varAuxR = y;
                    vary1 = vary2;
                    vary2 = varAuxR;
                    if(vary1 < vary2){
                        //dnl_taA = dnl_taA + 1.0;
                        //cout<< exp(a_ta) << "   " <<dnl_taA<<endl;
                        //exit(0);
                        break;
                    }
                    
                    evol_fnl(a_ta, dnli, ddnli, dnl_taA);
                    t_ta=a_ta;
                    y = log((dnl_taA + 1)/pow(exp(a_ta),3));
                }
            //-----------------------------------------------------

                cosmologicFuntion(t0, E_a_c, Eos_c, Omega_M_c, Omega_DE_c, ga_c);
                cosmologicFuntion(a_ta, E_a_ta, Eos_ta, Omega_M_ta, Omega_DE_ta, ga_ta);
                dnl_taA = (dnl_taA + 1.0);
                //n1 = (1.0+3.0*Eos_c)*(Om_w0*ga_c)/(dnl_taA*Om_m0*pow(exp(a_ta),-3));
                //n2 = 2.0*(1.0+3.0*Eos_ta)*(Om_w0*ga_ta)/(dnl_taA*Om_m0*pow(exp(a_ta),-3));
                n1 = -(1.0+3.0*Eos_c)*(Om_w0*ga_c)/(dnl_taA*Om_m0*pow(exp(a_ta),-3));
                n2 = -(1.0+3.0*Eos_ta)*(Om_w0*ga_ta)/(dnl_taA*Om_m0*pow(exp(a_ta),-3));
                //y_vir = (1.0-(n1/4.0))/(2.0+n2-(3.0*n1/4.0));
                y_vir = (1.0-(n1/2.0))/(2.0+n2-(3.0*n1/2.0));
                Delta_vir = (dnl_taA)*(pow((exp(t0)/exp(a_ta)),3))*pow(y_vir,-3.0);
                //------------------------------------------------
                //Auxiliar, colapso tarda menos
                //------------------------------------------------
                varAux = dnli;
                var1 = var2;
                var2 = varAux; 
                dnli += (var2-var1);//*0.97;
                //------------------------------------------------
                break;
            } 

            dnli += (2.5e-9)*exp(log(dmax)+log(dmax/dmin)*((j-1)*1e-10/(number-1)));    
            //dnli += (1.0e-9)*exp(log(dmax)+log(dmax/dmin)*((j-1)*1e-7/(number-1)));                

        }

    }

//--------------Halos-----------------------------------------------------

    double gm(double M){ //Función gamma(M)
        double m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0;
        double h=0.6766, r;
        present_data(m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0);
        double M8=(6e14)*m0_rho, Gm;
        Gm = m0_rho*h*exp(-b0_rho-b0_rho/m0_rho);
        r = (0.3*Gm+0.2)*(2.92+(1.0/3.0)*log(M/M8));
        return r;
    }

    double sigmaM0(double M){ //Sigma(M,z=0)
        double dc0L=1.67575244128603;
        double sigma8 = (d0/dc0L)*0.8083, sgM, M8;
        double m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0;
        present_data(m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0);
        M8=(6e14)*m0_rho;
        sgM = sigma8*pow(M/M8, -gm(M)/3.0);
        return sgM;
    }

    double sigmaMz(double M){ //Sigma(M,z)
        double smz;
        smz = sigmaM0(M)*D; //D(t0)
        return smz;
    }

    double dndM(double M){ //Densidad crítica
        double m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0;
        double dndm, Gm, rho_c0=8.53196555e10, h=0.6766, M8; //rho_c0 (densidad crítica) h^-1 M_solares Mpc^-3
        present_data(m0_rho, r0_rho,de0_rho,k0_rho,b0_rho,w0);
        Gm = m0_rho*h*exp(-b0_rho-b0_rho/m0_rho);
        M8=(6e14)*m0_rho;
        dndm=sqrt(2.0/M_PI)*(rho_c0/(3.0*pow(M,2)))*(dc/sigmaMz(M))*(gm(M)+(log(M/M8)/3.0)*(0.2+0.3*Gm))*exp(-pow(dc,2.0)/(2.0*pow(sigmaMz(M),2)));
        return dndm;
    }

    double n(double M){ //Densidad crítica integrada
        double n;
        n=integrateTz0(M, M_inf, dndM, 1000);
        return n;
    }

    double E_ia(double lna){
        double E_a, Eos, Omega_M, Omega_DE, ga, E;
        cosmologicFuntion(lna, E_a, Eos, Omega_M, Omega_DE, ga);
        E = 1.0/(exp(lna)*E_a);
        return E;
    }

    double evol_vol(double lna){ //dv/dzdOmega(H0/c)^3
        double v_int, tf=log(1.0), vol, cH3=8.80275955e10; //cH3 = (c/H0)^3
        integrateTz2(lna, tf, E_ia, 1000, v_int);
        vol = (cH3)*(E_ia(lna)*exp(lna))*pow(v_int,2); //*pow(integrateSmp(tf, lna, E_ia, 100),2);//
        return vol;
    }

    //double N_bin(double t0){//Número de halos N_bin
    void N_bin(double t0, double& n_bin, double& n_bin2, double& n_bin3){//Número de halos N_bin
 
        n_bin = 4.0*M_PI*integrateTz0(Minf1, Msup1, dndM, 1000)*evol_vol(t0);
        n_bin2 = 4.0*M_PI*integrateTz0(Msup1, Msup2, dndM, 1000)*evol_vol(t0);
        n_bin3 = 4.0*M_PI*integrateTz0(Msup2, Msup3, dndM, 1000)*evol_vol(t0);

        //return n_bin;
    }

    //double N_inf(double t0){//Número de halos integrado N_bin [M_0, M_infinita]
    void N_inf(double t0, double& n_inf, double& n_inf2, double& n_inf3){
    
        n_inf = 4.0*M_PI*n(M_0)*evol_vol(t0);
        n_inf2 = 4.0*M_PI*n(M_1)*evol_vol(t0);
        n_inf3 = 4.0*M_PI*n(M_2)*evol_vol(t0);

        //return n_inf;
    }

    /* void N(double t0, double& N1, double& N2, double& N3){
        double tf = log(1.0);
        N1 = integrateTz0(t0, tf, n_inf, 1000);
        N2 = integrateTz0(t0, tf, n_inf2, 1000);
        N3 = integrateTz0(t0, tf, n_inf3, 1000);
    } */

    /* double N(double t0){//Número de halos integrado N
        double n_bi, tf=log(1.0);
        //integrateTz(M_0, M_inf, dndM, 1000);
        n_bi =integrateTz0(t0, tf, N_inf, 1000);
        return n_bi;
    } */

    //------------------------------------------------------------------------

} //cierra scm


