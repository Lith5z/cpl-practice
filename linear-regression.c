#include <stdio.h>
#include <stdlib.h>

#define DELTA 0.001

int num,steps;
double *x_data,*y_data;

//sum of the squares of the vertical (plumb) distances from points to the line
double cost(double k,double b) {
    double res = 0;
    for (int i=0; i < num; i++) {
        res += (y_data[i] - (k * x_data[i] + b)) * (y_data[i] - (k * x_data[i] + b));
    }
    return res;
}
double k_derivative(double (*func)(double,double),double k,double b) {
    return (func(k + DELTA,b) - func(k,b)) / DELTA;
}
double b_derivative(double (*func)(double,double),double k,double b) {
    return (func(k,b + DELTA) - func(k,b)) / DELTA;
}

int main() {
    printf("Enter the data amount and steps:");
    scanf("%d %d",&num,&steps);
    x_data = malloc(num * sizeof(double));
    y_data = malloc(num * sizeof(double));
    double k_standard,b_standard;
    if (x_data == NULL || y_data == NULL) return -1;

    freopen("output.txt","r",stdin);
    for (int i=0; i < num; i++) {
        scanf("%lf %lf",&x_data[i],&y_data[i]);
    }
    scanf("%lf %lf",&k_standard,&b_standard);

    //initialize the k and b, using the first two points
    double k = (y_data[1] - y_data[0]) / (x_data[1] - x_data[0]);
    double b = y_data[0] - k * x_data[0];

    printf("d(cost) / d(k):%lf\n",k_derivative(cost,k,b));
    printf("d(cost) / d(b):%lf\n",b_derivative(cost,k,b));
    return 0;
}