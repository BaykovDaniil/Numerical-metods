#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<fstream>
#define EPS 0.001 // погрешность сравнения вещественных чисел

using namespace std;


double func(double x);
double* get_random_massive(double a, double b, size_t N);
double* get_denominator(double* &x, size_t N);
double lagrange(double _x, double* &x, double* &y, size_t N, double* &denominator);