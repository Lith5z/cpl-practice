//stochastic gradient descent

#include <stdio.h>
#include <stdlib.h>

//for derivative
#define DELTA 0.0001
//exp: rate ~ 1 / (x_max * y_max)
#define RATE 0.0001
#define LIMIT 0.001

int num,steps,batch_size;
double *x_data,*y_data;

//sum of the squares of the vertical (plumb) distances from points to the line
double total_cost(double k,double b) {
    double res = 0;
    for (int i=0; i < num; i++) {
        res += (y_data[i] - (k * x_data[i] + b)) * (y_data[i] - (k * x_data[i] + b));
    }
    return res;
}
double single_cost(double k,double b,int index) {
    return (y_data[index] - (k * x_data[index] + b)) * (y_data[index] - (k * x_data[index] + b));
}
double k_derivative(double (*func)(double,double,int),double k,double b,int start_index) {
    double res = 0;
    for (int i = start_index; i < start_index + batch_size; i++) {
        double tmp = (single_cost(k+DELTA,b,i) - single_cost(k,b,i)) / DELTA;
        res += tmp;
    }
    return res;
}
double b_derivative(double (*func)(double,double,int),double k,double b,int start_index) {
    double res = 0;
    for (int i = start_index; i < start_index + batch_size; i++) {
        double tmp = (single_cost(k,b+DELTA,i) - single_cost(k,b,i)) / DELTA;
        res += tmp;
    }
    return res;
}
double diff_acc(double accurate,double res) {
    return (res - accurate) / accurate * 100;
}

int main() {
    printf("Enter the data amount and steps:");
    scanf("%d %d",&num,&steps);
    printf("Enter the mini batch_size:");
    scanf("%d",&batch_size);
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
    double diff = (diff_acc(k_standard,k) + diff_acc(b_standard,b)) / 2;

    for (int i=1; i <= steps; i++) {
        printf("Steps:%d, batch_size:%d...\n",i,batch_size);
        double avg_cost = total_cost(k,b) / num;
        for (int p=0; p * batch_size < num; p++) {
            avg_cost = total_cost(k,b) / num;
            double k_d = k_derivative(single_cost,k,b,p * batch_size);
            double b_d = b_derivative(single_cost,k,b,p * batch_size);
            //printf("Batch[%d] d(cost) / d(k):%lf\n",p+1,k_d);
            //printf("Batch[%d] d(cost) / d(b):%lf\n",p+1,b_d);
            k -= k_d * RATE;
            b -= b_d * RATE;
        }
        printf("Total Cost:%lf\n",avg_cost);
        printf("\n");

        // if reach the LIMIT, stop at this step
        diff = (diff_acc(k_standard,k) + diff_acc(b_standard,b)) / 2;
        if (-LIMIT < diff && diff < LIMIT) {
            printf("\nReach the limit:%.3lf at step:%d\n",LIMIT,i);
            break;
        }
    }

    printf("The standard result: y=%.4lfx%+.4lf\n",k_standard,b_standard);
    printf("Gradient descent sesult: y=%.4lfx%+.4lf\n",k,b);
    printf("Diff:%.6lf%%\n",diff);

    return 0;
}