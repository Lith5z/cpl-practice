//Using GBK encoding
//A simple practice for c language programming

#include <stdio.h>
#include <stdlib.h>

#define ACCURACY 0.0001

enum MODE{
    CRAMER_MODE = 1, GAUSSIAN_MODE, DET_MODE, RANK_MODE
};

int debug_mode = 0;
int cal_mode = 1;
//m是矩阵的列数，无论是增广的还是一般的
int n,m;
//增广矩阵的秩
int extend_rank;
//系数矩阵的秩（去掉m列）
int coeff_rank;

int is_equal(double num,double target) {
    if (-ACCURACY <= num - target && num - target <= ACCURACY) return 1;
    else return 0;
}
//判断是否为齐次线性方程组，看常数项构成的向量是否全为0
int is_homogeneous(double **A) {
    int flag = 1;
    for (int i=0; i < n; i++) {
        if (!is_equal(A[i][m],0)) {
            flag = 0;
            break;
        }
    }
    return flag;
}
void swap(double *row1,double *row2);
//r1 = r1 - r2 * times
void elimate(double *row1,double *row2,double times,int start);
//r / divisor
void row_divide(double *row,double divisor,int start);

void show(double **matrix);
double** init_matrix(void);
void matrix_copy(double **org,double **dst);

//det函数会破坏原矩阵，请带入副本
double det(double **A) {
    int flag = 1;
    double res = 1;
    for (int i=0; i < n; i++) {
        if (debug_mode) {
            printf("\n只要求解以下矩阵的行列式：");
            show(A);
        }
        if (is_equal(A[i][i],0.0)) {
            int has_swap = 0;
            for (int j=i+1; j < n; j++) {
                if (!is_equal(A[j][i],0)) {
                    has_swap = 1;
                    flag *= -1;
                    swap(A[i],A[j]);
                    break;
                }
            }
            if (has_swap == 0) return 0;
        }
        for (int i_pos=i+1; i_pos < n; i_pos++) {
            if (!is_equal(A[i_pos][i],0)) elimate(A[i_pos],A[i],A[i_pos][i] / A[i][i],i);
        }
    }
    for (int i=0; i < n; i++) res *= A[i][i];
    if (debug_mode) printf("所以最后行列式的值为%lf\n",res * flag);
    return res * flag;
}
//A本身就是增广矩阵
void replace_column(double **A,int col) {
    for (int i=0; i < n; i++) {
        A[i][col] = A[i][m-1];
    }
}
void cramer_solve(double **A,double *res) {
    double **tempA = init_matrix();
    matrix_copy(A,tempA);
    double detAb = det(tempA);
    free(tempA);
    if (is_equal(detAb,0)) {
        printf("\n系数矩阵的行列式为0，克莱姆法则无法判断解的情况！\n");
        system("pause");
        exit(0);
    }

    for (int i=0; i < n; i++) {
        double **replaceA = init_matrix();
        matrix_copy(A,replaceA);
        replace_column(replaceA,i);
        res[i] = det(replaceA) / detAb;
        if (debug_mode) printf("\n第%d个未知数为%lf\n",i+1,res[i]);
        free(replaceA);
    }
}
//将增广矩阵化为行简化 同时维护两个rank，会修改带入的矩阵
void row_simplify(double **A) {
    //下面一行行的简化，i_cur和j_cur代表非零元的位置
    if (debug_mode) show(A);
    int i_cur = 0,j_cur = 0;
    while (i_cur < n && j_cur < m) {
        //增广矩阵 m 至少为 2
        if (j_cur == m - 2) coeff_rank = i_cur;
        if (is_equal(A[i_cur][j_cur],0.0)) {
            int has_swap = 0;
            for (int i=i_cur+1; i < n; i++) {
                if (!is_equal(A[i][j_cur],0)) {
                    has_swap = 1;
                    swap(A[i],A[i_cur]);
                    break;
                }
            }
            if (has_swap == 0) {
                //说明整个下方的列都是零
                j_cur++;
                continue;
            }
        }
        //assert 目前位置是非零元
        if (!is_equal(A[i_cur][j_cur],1.0)) row_divide(A[i_cur],A[i_cur][j_cur],j_cur);
        for (int i=0; i < i_cur; i++) {
            if (!is_equal(A[i][j_cur],0.0)) {
                elimate(A[i],A[i_cur],A[i][j_cur],j_cur);
            }
        }
        for (int i=i_cur+1; i < n; i++) {
            if (!is_equal(A[i][j_cur],0.0)) {
                elimate(A[i],A[i_cur],A[i][j_cur],j_cur);
            }
        }
        i_cur++;
        j_cur++;
    }
    if (j_cur < m) coeff_rank = i_cur;
    extend_rank = i_cur;
}
void gaussian_solve(double **A) {
    int is_homogen = is_homogeneous(A);
    row_simplify(A);

    if (is_homogen) {

    }
    else {

    }
}

int main() {
    //freopen("input.txt","r",stdin);
    int temp;
    printf("本计算器适用于 1.求解线性方程组Ax = b 2.计算方阵的行列式 3.将矩阵化为行简化阶梯型并求矩阵的秩 \n此外，本计算器的精度为 0.0001 ，如果数值过小，结果可能是错误的！\n如果你已经知晓以上内容，请输入 1 ：");
    scanf("%d",&temp);
    while (temp != 1) {
        printf("如果你已经知晓以上内容，请输入 1 ：");
        scanf("%d",&temp);
    }
    scanf("%*[^\n]");

    printf("\n计算模式选择：\n使用克莱姆法则，输入 1 ；使用高斯消元法，输入 2 ；\n");
    printf("计算方阵的行列式，输入 3 ； 将矩阵化为行简化阶梯型 并求矩阵的秩，输入 4 \n请输入：");
    scanf("%d",&cal_mode);
    while (cal_mode < 1 || cal_mode > 4) {
        printf("计算模式范围错误！请输入 1-4 内的整数：");
        scanf("%d",&cal_mode);
    }
    scanf("%*[^\n]");

    printf("\n是否开启监视模式？\n此模式下，计算过程会尽可能显示，请输入 0 或 1 \n请输入：");
    scanf("%d",&debug_mode);
    while (debug_mode != 0 && debug_mode != 1) {
        printf("输入范围错误！请输入 0 或 1 ：");
        scanf("%d",&debug_mode);
    }
    scanf("%*[^\n]");

    switch (cal_mode) {
        case CRAMER_MODE: case DET_MODE: {
            if (cal_mode == CRAMER_MODE) printf("\n请输入 系数 矩阵的阶数：");
            else printf("\n请输入方阵的阶数：");
            scanf("%d",&n);
            scanf("%*[^\n]");
            while (n < 1) {
                printf("输入范围错误！请输入一个正数：");
                scanf("%d",&n);
                scanf("%*[^\n]");
            }
            if (cal_mode == CRAMER_MODE) m = n + 1;
            else m = n;
            break;
        }
        case 2: case 4: {
            if (cal_mode == GAUSSIAN_MODE) printf("\n请输入 增广 矩阵 A_n*m的行数n：");
            else printf("\n请输入矩阵 A_n*m的行数n：");
            scanf("%d",&n);
            scanf("%*[^\n]");
            while (n < 1) {
                printf("输入范围错误！请输入一个正数：");
                scanf("%d",&n);
                scanf("%*[^\n]");
            }
            if (cal_mode == GAUSSIAN_MODE) printf("请输入 增广 矩阵 A_n*m的列数m：");
            else printf("请输入矩阵 A_n*m的列数m：");
            scanf("%d",&m);
            scanf("%*[^\n]");
            while (m < 1) {
                printf("输入范围错误！请输入一个正数：");
                scanf("%d",&m);
                scanf("%*[^\n]");
            }
            break;
        }
        default: {
            printf("\n这怎么可能？？\n");
            system("pause");
            return 1;
        }
    }

    double **matrix = init_matrix();

    switch (cal_mode) {
        case CRAMER_MODE: case GAUSSIAN_MODE: printf("\n请逐行输入矩阵(A b)：\n");break;
        case DET_MODE: printf("\n请逐行输入方阵：\n");break;
        case RANK_MODE: printf("\n请逐行输入矩阵：\n");break;
        default: {
            printf("\n这怎么可能？？\n");
            free(matrix);
            system("pause");
            return 1;
        }
    }
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            scanf("%lf",&matrix[i][j]);
        }
    }

    show(matrix);

    switch (cal_mode) {
        case CRAMER_MODE: {
            double *res_vector = malloc(n * sizeof(double));
            if (res_vector == NULL) {
                printf("\n内存分配失败！\n");
                system("pause");
                exit(0);
            }
            cramer_solve(matrix,res_vector);
            printf("\n所求解向量为\n");
            for (int i=0; i < n; i++) {
                printf("%lf\n",res_vector[i]);
            }
            free(res_vector);
            break;
        }
        case GAUSSIAN_MODE: {
            gaussian_solve(matrix);
            show(matrix);
            break;
        }
        case DET_MODE: {
            printf("\n所求方阵的行列式为%lf\n",det(matrix));
            break;
        }
        case RANK_MODE: {
            row_simplify(matrix);
            show(matrix);
            printf("该矩阵的秩为%d %d\n",coeff_rank,extend_rank);
            break;
        }
        default: {
            printf("\n这怎么可能？？\n");
            free(matrix);
            system("pause");
            return 1;
        }
    }

    free(matrix);
    system("pause");
    return 0;
}



void show(double **matrix) {
    if (cal_mode == CRAMER_MODE || cal_mode == GAUSSIAN_MODE) printf("\n当前增广矩阵为\n");
    else printf("\n当前矩阵为\n");
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            printf("%.3lf ",matrix[i][j]);
        }
        printf("\n");
    }
}
double** init_matrix(void) {
    double **matrix = malloc(n * sizeof(double*));
    if (matrix == NULL) {
        printf("\n内存分配失败！\n");
        system("pause");
        exit(0);
    }
    for (int i=0; i < n; i++) {
        matrix[i] = malloc(m * sizeof(double));
        if (matrix[i] == NULL) {
            printf("\n内存分配失败！\n");
            system("pause");
            exit(0);
        }
    }
    return matrix;
}
void matrix_copy(double **org,double **dst) {
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            dst[i][j] = org[i][j];
        }
    }
}
void swap(double *row1,double *row2) {
    for (int j=0; j < m; j++) {
        double temp = row1[j];
        row1[j] = row2[j];
        row2[j] = temp;
    }
}
void elimate(double *row1,double *row2,double times,int start) {
    for (int j=start; j < m; j++) {
        row1[j] -= row2[j] * times;
    }
}
void row_divide(double *row,double divisor,int start) {
    for (int j=start; j < m; j++) {
        row[j] /= divisor;
    }
}