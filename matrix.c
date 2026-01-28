#include <stdio.h>

#define LEN 100

int main() {
    puts("Let matrix A_m*n, B_n*p, please enter in the \"m n p\":");
    int m,n,p;
    scanf("%d%d%d",&m,&n,&p);

    int A[LEN][LEN],B[LEN][LEN];
    int C[LEN][LEN] = {0};
    puts("Enter in the A matrix row by row:");
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            scanf("%d",&A[i][j]);
        }
    }
    puts("Enter in the B matrix row by row:");
    for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
            scanf("%d",&B[i][j]);
        }
    }

    //j=column,i=row
    for (int j=0; j<p; j++) {
        for (int i=0; i<m; i++) {
            for (int k=0; k<n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    puts("A * B = C matrix row by row:");
    for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
            printf("%d ",C[i][j]);
        }
        printf("\n");
    }

    return 0;
}