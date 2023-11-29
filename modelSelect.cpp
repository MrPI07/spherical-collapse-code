#include "lib/DE-models.cpp"

//Modelos descritos en la libreria DE-models.cpp

using namespace EdS;
//using namespace LCDM;
//using namespace CPL;
//using namespace CPLF;
//using namespace neDE;
//using namespace EXP2;
//using namespace devgf;
//using namespace nonabelian;
//using namespace twoform;

int cosm = 0; //Imprimir datos cosmologicos, 1 si -- 0 no.
int ddelta = 0; //Imprimir evolución de perturbaciones: datos-delta.dat, 1 si -- 0 n0.
double dc, d0, D; //Variables globales del ámbito scm
double ddnli, dnli, dnl, dl, dnli0;
double y_vir, Delta_vir, dnl_taA, a_ta, t_ta=log(1.0);
double Minf1=1e13, Msup1=1e14, Msup2=1e15, Msup3=1e16; // 1:[Minf1, Msup1], 2:[Msup1, Msup2], 3:[Msup2, Msup3] Masa para el numero de halos
double n_bin, n_bin2, n_bin3;
double M_0=1e13,M_1=1e14,M_2=1e15, M_inf = 1e16; //Rango de masa para el número integrado de halos
double n_inf, n_inf2, n_inf3;