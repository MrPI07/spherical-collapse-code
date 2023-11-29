Programa que calcula la evolución de las perturbaciones esféricas de densidad
usando el modelo de colapso esférico. Los modelos cosmológicos se guardan en
el archivo lib/DE-models.cpp y se enlazan con las ecuaciones del modelo de 
colapso esférico (SCM) mediante los archivos lib/scm.cpp y modelSelect.cpp, en 
este último se debe seleccionar el modelo al que se va a calcular las cantidades 
del SCM. 

En el archivo main.cpp se han escrito las rutinas que realizan el calculo de las
cantidades delta(a) y delta_c(a) para el modelo cosmológico seleccionado. El
programa funciona a una velocidad moderada debido a la precisión con la que se 
hacen los calculos.

Modelos cosmológicos en la carpeta lib/DE-models.cpp y se deben seleccionar usnado
el archivo modelSelect.cpp;

Modelos EdS.
Modelos LCDM.
Modelo CPL.
MOdelo CPLF (CPL Phanton).
Modelo neDE: Energía oscura temprana.
Modelo 2EXP: doble exponecial.
Modelo con ecuación de estado dinámica (numérica) W. "Modelo Melissa"; devgf.
Modelo con ecuación de estado dinámica (numérica) W. "Non Abelian"; nonabelian.
Modelo con ecuación de estado dinámica (numérica) W. "Non Abelian"; twoform.

File compile main.cpp

Compiling: 

make comp (or make):
    build rk45f.o
    compile main.cpp
    link main.o with rk45f.o

make run:
    Run a program: main


make clean:
    Remove file: *.o main

