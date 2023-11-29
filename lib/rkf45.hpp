void r4_fehl ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] );
int r4_rkf45 ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double yp[], double *t, double tout, double *relerr, double abserr, 
  int flag );
double r4_sign ( double x );

void r8_fehl ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] );
int r8_rkf45 ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double yp[], double *t, double tout, double *relerr, double abserr, 
  int flag );
double r8_sign1 ( double x );

void timestamp ( );
