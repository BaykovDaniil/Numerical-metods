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
double I_func(double x);
double True_Integral(double a, double b);
double* get_uniform_mesh(double a, double b, size_t M);
double* Int_Trap(double* &x, size_t M);
double* Int_Rect(double* &x, size_t M);
double* Int_Simpson(double* &x, size_t M);
double* Int_Newton(double* &x, size_t M);
double* Int_Gauss(double* &x, size_t M);
double relative_error(double x, double y);