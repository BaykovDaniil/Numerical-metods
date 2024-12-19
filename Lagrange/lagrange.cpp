#include"lagrange.h"

double func(double x){

	return cos(x);

}

// создает массив из N случайных вещественных чисел от a до b
// 1-е и N-е всегда a и b
// требует отчистки памяти
double* get_random_massive(double a, double b, size_t N){

	double* x = new double[N];
	double  r;
    bool    flag = 0; // хранит информацию о совпадении нового случайного числа с уже сгенерированными
    x[0] = a;
    for(size_t i=1; i<N-1; ++i){
    	while(true){
    		r = a + (b-a) * rand() / RAND_MAX;
    		if( abs(r-a) < EPS*(b-a) ){
    				flag = 1;
    			}
    		for(size_t j=0; j<i; ++j){
    			if( abs(r-x[j]) < EPS*(b-a) ){
    				flag = 1;
    			}
    		}
    		if( abs(r-b) < EPS*(b-a) ){
    				flag = 1;
    			}
    		if(flag == 0){
    			x[i] = r;
    			break;
    		}
    		else{
    			flag = 0;
    		}
    	}
    }
    x[N-1] = b;
    return x;

}

// создает массив повторяющихся в многочлене Лагранжа знаменателей
// требует отчисти памяти
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

	double result = 0, product;
	for(size_t i=0; i<N; ++i){
		product = 1.0;
		for(size_t j=0; j<N; ++j){
			if(i==j){
    			continue;
    		}
    		product *= (_x-x[j]);
		}
		result += y[i] * product / denominator[i];
	}
	return result;

}
