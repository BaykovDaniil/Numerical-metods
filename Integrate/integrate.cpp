#include"integrate.h"

double func(double x){

	return pow(x, 7);

}

double I_func(double x){

    return pow(x, 8) / 8.0;

}

double True_Integral(double a, double b){

    return I_func(b) - I_func(a);

}

double* Int_Trap(double* &x, size_t M){

    double  sum_1 = 0, 
            sum_2 = 0,
            h     = x[1] - x[0];
    sum_1 += func(x[0]) / 2.0;
    for(size_t i=1; i<M-1; ++i){
        sum_1 += func(x[i]);
        sum_2 += func(x[i] - h/2.0);
    }
    sum_1 += func(x[M - 1]) / 2.0;
    sum_2 += func(x[M-1] - h/2.0);
    sum_2 += sum_1;

    double* answer = new double[2];
    answer[0] = h * sum_1;
    answer[1] = h/2.0 * sum_2;
    return answer;

}

double* Int_Rect(double* &x, size_t M){

    double  sum_1 = 0, 
            sum_2 = 0,
            h     = x[1] - x[0];
    for(size_t i=0; i<M-1; ++i){
        sum_1 += func(x[i] + h/2.0);
        sum_2 += func(x[i] + h/4.0);
        sum_2 += func(x[i] + 3.0*h/4.0);
    }

    double* answer = new double[2];
    answer[0] = h * sum_1;
    answer[1] = h/2.0 * sum_2;
    return answer;

}

double* Int_Simpson(double* &x, size_t M){

    double  sum_1           = 0,
            sum_2           = 0,
            h               = x[1] - x[0],
            half_h          = h / 2.0,
            quarter_h       = h / 4.0,
            subsum_1_c_is_1 = 0,
            subsum_1_c_is_4 = 0,
            subsum_1_c_is_2 = 0,
            subsum_2_c_is_1 = 0,
            subsum_2_c_is_4 = 0,
            subsum_2_c_is_2 = 0;
    subsum_1_c_is_1 += func(x[0]);
    subsum_2_c_is_1 += func(x[0]);
    // sum_1 += func(x[0]);
    // sum_2 += func(x[0]);
    for(size_t i=0; i<M-2; ++i){
        subsum_1_c_is_4 += func(x[i] + half_h);
        subsum_1_c_is_2 += func(x[i] + h);

        subsum_2_c_is_4 += func(x[i] + quarter_h);
        subsum_2_c_is_2 += func(x[i] + half_h);
        subsum_2_c_is_4 += func(x[i] + 3.0 * quarter_h);
        subsum_2_c_is_2 += func(x[i] + h);

        // sum_1 += 4.0 * func(x[i] + h/2.0);
        // sum_1 += 2.0 * func(x[i] + h);

        // sum_2 += 4.0 * func(x[i] + h/4.0);
        // sum_2 += 2.0 * func(x[i] + h/2.0);
        // sum_2 += 4.0 * func(x[i] + 3.0*h/4.0);
        // sum_2 += 2.0 * func(x[i] + h);
    }
    subsum_1_c_is_4 += func(x[M-1] - half_h);
    subsum_1_c_is_1 += func(x[M-1]);
    sum_1 = subsum_1_c_is_1 +
            4.0 * subsum_1_c_is_4 +
            2.0 * subsum_1_c_is_2;

    subsum_2_c_is_4 += func(x[M-2] + quarter_h);
    subsum_2_c_is_2 += func(x[M-2] + half_h);
    subsum_2_c_is_4 += func(x[M-2] + 3.0 * quarter_h);
    subsum_2_c_is_1 += func(x[M-2] + h);
    sum_2 = subsum_2_c_is_1 +
            4.0 * subsum_2_c_is_4 +
            2.0 * subsum_2_c_is_2;

    // sum_1 += 4.0 * func(x[M - 1] - h/2.0);
    // sum_1 += func(x[M - 1]);

    // sum_2 += 4.0 * func(x[M-2] + h/4.0);
    // sum_2 += 2.0 * func(x[M-2] + h/2.0);
    // sum_2 += 4.0 * func(x[M-2] + 3.0*h/4.0);
    // sum_2 += func(x[M-2] + h);

    double* answer = new double[2];
    answer[0] = half_h * sum_1 / 3.0;
    answer[1] = quarter_h * sum_2 / 3.0;
    return answer;

}

double* Int_Newton(double* &x, size_t M){

    double  sum_1               = 0, 
            sum_2               = 0,
            h                   = x[1] - x[0],
            quarter_h           = h / 4.0,
            quaver_h            = h / 8.0,
            subsum_1_c_is_7     = 0,
            subsum_1_c_is_32    = 0,
            subsum_1_c_is_12    = 0,
            subsum_1_c_is_14    = 0,
            subsum_2_c_is_7     = 0,
            subsum_2_c_is_32    = 0,
            subsum_2_c_is_12    = 0,
            subsum_2_c_is_14    = 0;
    subsum_1_c_is_7 += func(x[0]);
    subsum_2_c_is_7 += func(x[0]);
    // sum_1 += 7.0 * func(x[0]);
    // sum_2 += 7.0 * func(x[0]);
    for(size_t i=0; i<M-2; ++i){
        subsum_1_c_is_32 += func(x[i] + quarter_h);
        subsum_1_c_is_12 += func(x[i] + 2.0 * quarter_h);
        subsum_1_c_is_32 += func(x[i] + 3.0 * quarter_h);
        subsum_1_c_is_14 += func(x[i] + h);

        subsum_2_c_is_32 += func(x[i] + quaver_h);
        subsum_2_c_is_12 += func(x[i] + quarter_h);
        subsum_2_c_is_32 += func(x[i] + 3.0 * quaver_h);
        subsum_2_c_is_14 += func(x[i] + 4.0 * quaver_h);
        subsum_2_c_is_32 += func(x[i] + 5.0 * quaver_h);
        subsum_2_c_is_12 += func(x[i] + 6.0 * quaver_h);
        subsum_2_c_is_32 += func(x[i] + 7.0 * quaver_h);
        subsum_2_c_is_14 += func(x[i] + h);

        // sum_1 += 32.0 * func(x[i] + h/4.0);
        // sum_1 += 12.0 * func(x[i] + h/2.0);
        // sum_1 += 32.0 * func(x[i] + 3.0*h/4.0);
        // sum_1 += 14.0 * func(x[i] + h);

        // sum_2 += 32.0 * func(x[i] + h/8.0);
        // sum_2 += 12.0 * func(x[i] + h/4.0);
        // sum_2 += 32.0 * func(x[i] + 3.0*h/8.0);
        // sum_2 += 14.0 * func(x[i] + h/2.0);
        // sum_2 += 32.0 * func(x[i] + 5.0*h/8.0);
        // sum_2 += 12.0 * func(x[i] + 3.0*h/4.0);
        // sum_2 += 32.0 * func(x[i] + 7.0*h/8.0);
        // sum_2 += 14.0 * func(x[i] + h);
    }
    subsum_1_c_is_32 += func(x[M-2] + quarter_h);
    subsum_1_c_is_12 += func(x[M-2] + 2.0 * quarter_h);
    subsum_1_c_is_32 += func(x[M-2] + 3.0 * quarter_h);
    subsum_1_c_is_7  += func(x[M-2] + h);
    sum_1 = 7.0 * subsum_1_c_is_7 +
            32.0 * subsum_1_c_is_32 +
            12.0 * subsum_1_c_is_12 +
            14.0 * subsum_1_c_is_14;

    subsum_2_c_is_32 += func(x[M-2] + quaver_h);
    subsum_2_c_is_12 += func(x[M-2] + quarter_h);
    subsum_2_c_is_32 += func(x[M-2] + 3.0 * quaver_h);
    subsum_2_c_is_14 += func(x[M-2] + 4.0 * quaver_h);
    subsum_2_c_is_32 += func(x[M-2] + 5.0 * quaver_h);
    subsum_2_c_is_12 += func(x[M-2] + 6.0 * quaver_h);
    subsum_2_c_is_32 += func(x[M-2] + 7.0 * quaver_h);
    subsum_2_c_is_7  += func(x[M-2] + h);
    sum_2 = 7.0 * subsum_2_c_is_7 +
            32.0 * subsum_2_c_is_32 +
            12.0 * subsum_2_c_is_12 +
            14.0 * subsum_2_c_is_14;
    // sum_1 += 32.0 * func(x[M-2] + h/4.0);
    // sum_1 += 12.0 * func(x[M-2] + h/2.0);
    // sum_1 += 32.0 * func(x[M-2] + 3.0*h/4.0);
    // sum_1 += 7.0 * func(x[M-2] + h);

    // sum_2 += 32.0 * func(x[M-2] + h/8.0);
    // sum_2 += 12.0 * func(x[M-2] + h/4.0);
    // sum_2 += 32.0 * func(x[M-2] + 3.0*h/8.0);
    // sum_2 += 14.0 * func(x[M-2] + h/2.0);
    // sum_2 += 32.0 * func(x[M-2] + 5.0*h/8.0);
    // sum_2 += 12.0 * func(x[M-2] + 3.0*h/4.0);
    // sum_2 += 32.0 * func(x[M-2] + 7.0*h/8.0);
    // sum_2 += 7.0 * func(x[M-2] + h);

    double* answer = new double[2];
    answer[0] = 2.0 * quarter_h * sum_1 / 45.0;
    answer[1] = 2.0 * quaver_h * sum_2 / 45.0;
    return answer;

}

double* Int_Gauss(double* &x, size_t M){

    double  sum_1       = 0, 
            sum_2       = 0,
            subsum_1_c1 = 0,
            subsum_1_c2 = 0,
            subsum_2_c1 = 0,
            subsum_2_c2 = 0,
            h           = x[1] - x[0],
            half_h      = h / 2.0,
            quarter_h   = h / 4.0,
            coef1       = 8.0 / 9.0,
            coef2       = 5.0 / 9.0,
            shift1      = pow(3.0 / 5.0, 1.0 / 2.0) * h / 2.0,
            shift2      = shift1 / 2.0;
    for(size_t i=0; i<M-1; ++i){

        subsum_1_c1 += func(x[i] + half_h);
        subsum_1_c2 += func(x[i] + half_h - shift1);
        subsum_1_c2 += func(x[i] + half_h + shift1);

        subsum_2_c1 += func(x[i] + quarter_h);
        subsum_2_c1 += func(x[i] + 3.0*quarter_h);
        subsum_2_c2 += func(x[i] + quarter_h - shift2);
        subsum_2_c2 += func(x[i] + quarter_h + shift2);
        subsum_2_c2 += func(x[i] + 3.0*quarter_h - shift2);
        subsum_2_c2 += func(x[i] + 3.0*quarter_h + shift2);

        // sum_1 += coef1 * func(x[i] + h/2.0);
        // sum_1 += coef2 * func(x[i] + h/2.0 - shift1);
        // sum_1 += coef2 * func(x[i] + h/2.0 + shift1);

        // sum_2 += coef1 * func(x[i] + h/4.0);
        // sum_2 += coef1 * func(x[i] + 3.0*h/4.0);
        // sum_2 += coef2 * func(x[i] + h/4.0 - shift2);
        // sum_2 += coef2 * func(x[i] + h/4.0 + shift2);
        // sum_2 += coef2 * func(x[i] + 3.0*h/4.0 - shift2);
        // sum_2 += coef2 * func(x[i] + 3.0*h/4.0 + shift2);
    }
    sum_1 = coef1 * subsum_1_c1 +
            coef2 * subsum_1_c2;
    sum_2 = coef1 * subsum_2_c1 +
            coef2 * subsum_2_c2;

    double* answer = new double[2];
    answer[0] = half_h * sum_1;
    answer[1] = quarter_h *sum_2;
    return answer;

}

// создает равномерный массив из M вещественных чисел от a до b
// 1-е и M-е всегда a и b
// требует отчистки памяти
double* get_uniform_mesh(double a, double b, size_t M){

	double* x = new double[M];
    double h = (b-a) / (M-1);
    for(size_t i=0; i<M; ++i){
    	x[i] = a + h * i ;
    }
    return x;

}

// относительная погрешность в процентах
// x - истинное значение
// y - приближенное
double relative_error(double x, double y){

    return abs(x-y) / abs(y) * 100.0;

}

