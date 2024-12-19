#include"derivative.h"

double func(double x){

	return sin(x);

}

double D_func(double x){

    return cos(x);

}

// создает равномерный массив из M вещественных чисел от a до b
// 1-е и M-е всегда a и b
// требует отчистки памяти
double* get_uniform_mesh(double a, double b, size_t M){

	double* x = new double[M];
    for(size_t i=0; i<M; ++i){
    	x[i] = a + (b-a) * i / (M-1);
    }
    return x;

}

// создает массив производных функции func в точках сетки x
// сетка долждна быть равномерной
// 2 порядок точности
// требует отчистки памяти
double* Derivative(double* &x, size_t M){

    double h        = x[1] - x[0];
    double* result  = new double[M];

    result[0] = (-3.0*func(x[0]) + 4.0*func(x[1]) - func(x[2]) ) / (2.0*h);
    for(size_t i=1; i<M-1; ++i){
        result[i] = (func(x[i+1]) - func(x[i-1])) / (2.0*h);
    }
    result[M-1] = (3.0*func(x[M-1]) - 4.0*func(x[M-2]) + func(x[M-3]) ) / (2.0*h);
    return result;

}

// уточняет D1 по D2 по второй формуле Рунге
// D2 должен быть в 2 раза мельче D1 (иначе изменить pow)
// требует отчистки памяти
double* Runge(double* &D1, double* &D2, size_t M){

    double* result = new double[M];
    for(size_t i=0; i<M; ++i){
        result[i] = D1[i] + (D1[i] - D2[2*i]) / (pow(1.0/2.0,2) - 1);
    }
    return result;

}

// возвращает массив главных членов погрешности по первой формуле Рунге
// D2 должен быть в 2 раза мельче D1 (иначе изменить pow)
// требует отчистки памяти
double* TheoreticalRungeError(double* &D1, double* &D2, size_t M){

    double* result = new double[M];
    for(size_t i=0; i<M; ++i){
        result[i] = (D1[i] - D2[2*i]) / (pow(1.0/2.0,2) - 1);
    }
    return result;

}


double l_1(double* &D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += abs(D[i]);
    }
    return result;

}


double l_2(double* &D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += pow(D[i], 2);
    }
    return pow(result, 1.0 / 2.0);

}


double l_inf(double* &D, size_t M){

    double result = abs(D[0]);
    for(size_t i=1; i<M; ++i){
        if(abs(D[i]) > result)
            result = abs(D[i]);
    }
    return result;

}


