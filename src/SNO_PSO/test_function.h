/* 
 * There are 10 basic functions (minimization) in this file for testing the metaheuristic algorithms.
 * The number before the name is the index for calling that function.
 * 1) Ackley 
 * 2) Griewank
 * 3) Bent Cigar
 * 4) Michalewicz
 * 5) Rosenbrock
 * 6) Schwefel 2.26
 * 7) Zakharov
 * 8) Happy Cat
 * 9) Rastrigin
 * 10) HGBat
 */

#ifndef __TEST_FUNC__
#define __TEST_FUNC__
#include <math.h>
#include <float.h>
#include <stdio.h>

double cal_test_function(const double *x, const int d, const int func_num);
void set_search_range(double *range_max, double *range_min, const int func_num);

double Ackley(const double *x, const int d);
double Griewank(const double *x, const int d);
double BentCigar(const double *x, const int d);
double Michalewicz(const double *x, const int d);
double Rosenbrock(const double *x, const int d);
double Schwefel(const double *x, const int d);
double Zakharov(const double *x, const int d);
double HappyCat(const double *x, const int d);
double Rastrigin(const double *x, const int d);
double HGBat(const double *x, const int d);



void set_search_range(double *range_max, double *range_min, const int func_num)
{
    switch (func_num) {
    case 1:
	*range_max =  32.768;
	*range_min = -32.768;
	printf("Initialize Ackley function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 2:
	*range_max =  600.0;
	*range_min = -600.0;
	printf("Initialize Griewank function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 3:
	*range_max =  100.0;
	*range_min = -100.0;
	printf("Initialize BentCigar function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 4:
	*range_max =  M_PI;
	*range_min =  0.0;
	printf("Initialize Michalewicz function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 5:
	*range_max =  10.0;
	*range_min = -10.0;
	printf("Initialize Rosenbrock function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 6:
	*range_max =  500.0;
	*range_min = -500.0;
	printf("Initialize Schwefel function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 7:
	*range_max =  10.0;
	*range_min = -10.0;
	printf("Initialize Zakharov function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 8:
	*range_max =  20.0;
	*range_min = -20.0;
	printf("Initialize HappyCat function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 9:
	*range_max =  5.12;
	*range_min = -5.12;
	printf("Initialize Rastrigin function with range [%f, %f]\n", *range_min, *range_max);
	break;
    case 10:
	*range_max =  15.0;
	*range_min = -15.0;
	printf("Initialize HGBat function with range [%f, %f]\n", *range_min, *range_max);
	break;
    default:
	printf("invalid function.\n");
	abort();
	break;
    }
}

double cal_test_function(const double *x, const int d, const int func_num)
{
    double f = DBL_MAX;
    switch (func_num) {
    case 1:
	f = Ackley(x, d);
	break;
    case 2:
	f = Griewank(x, d);
	break;
    case 3:
	f = BentCigar(x, d);
	break;
    case 4:
	f = Michalewicz(x, d);
	break;
    case 5:
	f = Rosenbrock(x, d);
	break;
    case 6:
	f = Schwefel(x, d);
	break;
    case 7:
	f = Zakharov(x, d);
	break;
    case 8:
	f = HappyCat(x, d);
	break;
    case 9:
	f = Rastrigin(x, d);
	break;
    case 10:
	f = HGBat(x, d);
	break;
    default:
	printf("invalid function.\n");
	break;
    }
    return f;
}

double Ackley(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * x[i];
        sum2 += cos(2.0 * M_PI * x[i]);
    }
    return (-20.0) * exp(-0.2 * sqrt(sum1 / double(d))) - exp(sum2 / double(d)) + 20.0 + exp(1.0);
}

double Griewank(const double *x, const int d)
{
    double sum1 = 0.0, product1 = 1.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * x[i];
        product1 *= cos(x[i] / sqrt(i + 1));
    }
    return (sum1 / 4000.0) - product1 + 1.0;
}

double BentCigar(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 1; i < d; ++i) {
        sum1 += x[i] * x[i];
    }
    return x[0] * x[0] + pow(10.0, 6) * sum1;
}

double Michalewicz(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += sin(x[i]) * pow(sin((double(i + 1) * x[i] * x[i]) / M_PI), 20.0);
    }
    return sum1 * (-1);
}

double Rosenbrock(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d - 1; ++i) {
        sum1 += 100.0 * (x[i + 1] - (x[i] * x[i])) * (x[i + 1] - (x[i] * x[i])) + ((x[i] - 1.0) * (x[i] - 1.0));
    }
    return sum1;
}

double Schwefel(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * double(d) - sum1;
}

double Zakharov(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * x[i];
        sum2 += 0.5 * (i + 1) * x[i];
    }
    return sum1 + pow(sum2, 2) + pow(sum2, 4);
}

double HappyCat(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * x[i];
        sum2 += x[i];
    }
    return pow(fabs(sum1 - double(d)), 0.25) + (0.5 * sum1 + sum2) / double(d) + 0.5;
}

double Rastrigin(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += (x[i] * x[i]) - (10.0 * cos(2.0 * M_PI * x[i]));
    }
    return sum1 + 10.0 * double(d);
}

double HGBat(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i) {
        sum1 += x[i] * x[i];
        sum2 += x[i];
    }
    return sqrt(fabs((sum1 * sum1) - (sum2 * sum2))) + (0.5 * sum1 + sum2) / double(d) + 0.5;
}


#endif
