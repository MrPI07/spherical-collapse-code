CC = g++
CFLAGS = -Wall #-I/usr/local/include
LINKS = -L/usr/local/lib
#LIBS = -lgsl -lgslcblas -lm

comp: #DE-models.cpp
	$(CC) -c -o lib/rkf45.o lib/rkf45.cpp
	$(CC) -c -o lib/ode.o lib/ode.cpp
	$(CC) -c -o lib/integral.o lib/integral.cpp
	$(CC) -c -o lib/spline.o lib/spline.cpp
	$(CC) $(CFLAGS) -c main.cpp
	$(CC) -o main main.o lib/spline.o lib/ode.o lib/rkf45.o
	
run: 
	./main

clean:
	rm -f *.o main


