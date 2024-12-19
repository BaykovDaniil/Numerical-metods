#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>


#define EPS 0.001 // погрешность сравнения вещественных чисел

using namespace std;


double func(double x);
double D_func(double x);
double* get_uniform_mesh(double a, double b, size_t M);
double* Derivative(double* &x, size_t M);
double* Runge(double* &D1, double* &D2, size_t M);
double* TheoreticalRungeError(double* &D1, double* &D2, size_t M);
double l_1(double* &D, size_t M);
double l_2(double* &D, size_t M);
double l_inf(double* &D, size_t M);