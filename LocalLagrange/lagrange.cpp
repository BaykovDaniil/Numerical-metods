
#include"lagrange.h"

double func(double x){

	return sin(5.0*x);

}

// создает массив из M равномерно распределенных вещественных чисел от a до b
// 1-е и M-е всегда a и b
// требует отчистки памяти
double* get_uniform_mesh(double a, double b, size_t M){

	double* x = new double[M];
    for(size_t i=0; i<M; ++i){
    	x[i] = a + (b-a) * i / (M-1);
    }
    return x;

}


double* get_denominator(double* &x, size_t N){

    double* denominator = new double[N];
    for(size_t i=0; i<N; ++i){
        double product = 1.0;
        for(size_t j=0; j<N; ++j){
            if(i==j){
                continue;
            }
            product *= x[i] - x[j];
        }
        denominator[i] = product;
    }
    return denominator;

}

// возвращает значение аппрокcимации в точке _x
// (x, y) - точки, по которым идет аппрокcимация
double lagrange(double _x, double* &x, double* &y, size_t N, double* &denominator){
    /* 
    N - степень многочлена
    */
	double result = 0.0, product;
    size_t num = floor( (_x - x[0]) / (x[N-1]-x[0]) ); // номер интервала
    double x_new = _x - (x[N-1]-x[0])*num; // локальные координаты - первый подотрезок
	for(size_t i=0; i<N; ++i){
		product = 1.0;
		for(size_t j=0; j<N; ++j){
			if(i==j){
    			continue;
    		}
    		product *= (x_new - x[j]);
		}
		result += y[i + (N-1)*num] * product / denominator[i];
	}
	return result;

}


double l_1(double* D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += abs(D[i]);
    }
    return result;

}


double l_2(double* D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += pow(D[i], 2);
    }
    return pow(result, 1.0 / 2.0);

}


double l_inf(double* D, size_t M){

    double result = abs(D[0]);
    for(size_t i=1; i<M; ++i){
        if(abs(D[i]) > result)
            result = abs(D[i]);
    }
    return result;

}
