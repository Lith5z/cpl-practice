#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM 10
//NUM should >= 3

int main() {
    double k,b;
    printf("[slope] [intercept]\nEnter the standard result:");
    scanf("%lf %lf",&k,&b);

    printf("The program will generate %d lines of x-y data in the output.txt file\n",NUM);
    printf("x: from 0 to 10\n");
    printf("and will produce a linear regression result at line %d\n",NUM+2);

    double x[NUM], y[NUM];

    fopen("output.txt","w");
    freopen("output.txt","w",stdout);
    srand((unsigned int)time(NULL));
    for (int i = 0; i < NUM; i++) {
        x[i] = rand() / 3276.7;
        printf("%.3lf",x[i]);
        srand((unsigned int)rand());
        y[i] = x[i] * k + b + rand() / 37267.0;
        printf(" %.3lf\n",y[i]);
        srand((unsigned int)rand());
    }

    //using the ordinary least squares method to calculate the best result
    double x_avg = 0, y_avg = 0, frac1 = 0, frac2 = 0;
    for (int i = 0; i < NUM; i++) {
        x_avg += x[i];
        y_avg += y[i];
        frac1 += x[i] * y[i];
        frac2 += x[i] * x[i];
    }
    x_avg /= NUM;
    y_avg /= NUM;
    double k_res = (frac1 - NUM * x_avg * y_avg) / (frac2 - NUM * x_avg * x_avg);
    double b_res = y_avg - k_res * x_avg;
    printf("\n%.3f %.3f",k_res,b_res);
    return 0;
}