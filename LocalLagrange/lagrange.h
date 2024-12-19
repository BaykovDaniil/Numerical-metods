#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<fstream>
#include<cmath>
#define EPS 0.001 // погрешность сравнения вещественных чисел

using namespace std;


double func(double x);
double* get_uniform_mesh(double a, double b, size_t M);
double* get_denominator(double* &x, size_t N);
double lagrange(double _x, double* &x, double* &y, size_t N, double* &denominator);
double l_1(double* D, size_t M);
double l_2(double* D, size_t M);
double l_inf(double* D, size_t M);