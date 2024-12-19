#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<ctime>
#include<fstream>
#include<cmath>

#define EPS 0.00001 // погрешность сравнения вещественных чисел

using namespace std;


double func(double x);
double* get_uniform_mesh(double a, double b, size_t M);
double* get_random_massive(double a, double b, size_t M);
double* get_denominator_old(double* &x, size_t N);
double scalar(double* &x, double* &y, size_t N);
double lagrange(double _x, double* &x, double* &y, size_t N, double* denominator);
double basis_func(double _x, size_t i, double* &x, size_t N, size_t K, double* &denominator);
void get_sistem(double** &сoef, size_t N, size_t K, size_t L, double* &B, double* &x, double* &denominator);
double linear_comb(double* coef, size_t L, double _x, double* &x, size_t N, size_t K, double* &denominator);
double l_1(double* D, size_t M);
double l_2(double* D, size_t M);
double l_inf(double* D, size_t M);
int SLAU(double **matrica_a,int n,double *massiv_b, double *x);
