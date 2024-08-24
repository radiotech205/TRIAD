#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

double dot_product(double *fv, double *sv);
void cross_product(double *fv, double *sv, double *dest);
double vector_normalazer(double* xyz);
void matrix_multiplication(double *fm, double *sm, double *dest);
void transponse_matrix(double *matrix);
void create_matrix_from_vectors(double *fv, double *sv, double *tv, double *dest);


void RotMatrixX(double angel, double* matrix);
void RotMatrixY(double angel, double* matrix);
void RotMatrixZ(double angel, double* matrix);

void dcm2quaternion(const double* dcm, double* quat);
void MatVectMult(const double* matrix, const double* vector, double* dest);
void PrintMatrix(const double* matrix);

double determinant3x3(double* matrix);
void quat2Euler(const double* quat, double* ra, double* dec, double* phi);
double quat_vector_norm(const double* quat);


#endif // MATRIX_H
