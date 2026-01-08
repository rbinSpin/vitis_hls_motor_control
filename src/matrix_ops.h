#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

void Gaussian_Elimination_3d(float MatrixA[3][3], float MatrixB[3][3], float Minv[3][3]);
void Gaussian_Elimination_4d(float MatrixA[4][4], float MatrixB[4][4], float Minv[4][4]);
void Gaussian_Elimination_T_3d(float MatrixA[3][3], float MatrixB[3][3], float Minv[3][3]);
void Gaussian_Elimination_T_4d(float MatrixA[4][4], float MatrixB[4][4], float Minv[4][4]);
void Multiply_3d(float MatA[3][3], float MatB[3][3], float MatC[3][3]);
void Multiply_4d(float MatA[4][4], float MatB[4][4], float MatC[4][4]);
void Add_3d(float MatA[3][3], float MatB[3][3], float MatC[3][3]);
void Add_4d(float MatA[4][4], float MatB[4][4], float MatC[4][4]);
void Transpose_3d(float MatA[3][3], float MatB[3][3]);
void Transpose_4d(float MatA[4][4], float MatB[4][4]);
void Swap(float &A,float &B);
void Saturate(float A,  float upper_limit, float lower_limit, float &B);

#endif // MATRIX_OPS_H