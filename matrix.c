#include "matrix.h"

double dot_product(double *fv, double *sv)
{
    return fv[0]*sv[0] + fv[1]*sv[1] + fv[2]*sv[2];
}

void cross_product(double *fv, double *sv, double *dest)
{
    dest[0] = (fv[1] * sv[2]) - (sv[1] * fv[2]);
    dest[1] = ((fv[0] * sv[2]) - (sv[0] * fv[2])) * -1;
    dest[2] = (fv[0] * sv[1]) - (sv[0] * fv[1]);
}

double vector_normalazer(double* xyz)
{
    double norma = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
    xyz[0] /= norma;
    xyz[1] /= norma;
    xyz[2] /= norma;
    return norma;
}

void matrix_multiplication(double *fm, double *sm, double *dest)
{
    dest[0] = ((fm[0] * sm[0]) + (fm[1] * sm[3]) + (fm[2] * sm[6]));
    dest[1] = ((fm[0] * sm[1]) + (fm[1] * sm[4]) + (fm[2] * sm[7]));
    dest[2] = ((fm[0] * sm[2]) + (fm[1] * sm[5]) + (fm[2] * sm[8]));

    dest[3] = ((fm[3] * sm[0]) + (fm[4] * sm[3]) + (fm[5] * sm[6]));
    dest[4] = ((fm[3] * sm[1]) + (fm[4] * sm[4]) + (fm[5] * sm[7]));
    dest[5] = ((fm[3] * sm[2]) + (fm[4] * sm[5]) + (fm[5] * sm[8]));

    dest[6] = ((fm[6] * sm[0]) + (fm[7] * sm[3]) + (fm[8] * sm[6]));
    dest[7] = ((fm[6] * sm[1]) + (fm[7] * sm[4]) + (fm[8] * sm[7]));
    dest[8] = ((fm[6] * sm[2]) + (fm[7] * sm[5]) + (fm[8] * sm[8]));
}

void transponse_matrix(double *matrix)
{
    double temp;

    temp = matrix[3];
    matrix[3] = matrix[1];
    matrix[1] = temp;

    temp = matrix[6];
    matrix[6] = matrix[2];
    matrix[2] = temp;

    temp = matrix[7];
    matrix[7] = matrix[5];
    matrix[5] = temp;
}

void create_matrix_from_vectors(double *fv, double *sv, double *tv, double *dest)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(j == 0)
                *dest = fv[i];

            else if(j == 1)
                *dest = sv[i];

            else if(j == 2)
                *dest = tv[i];

            dest++;
        }
    }
}

void RotMatrixX(double angel, double* matrix)
{
    matrix[0] = 1.0;
    matrix[1] = 0.0;
    matrix[2] = 0.0;

    matrix[3] = 0.0;
    matrix[4] = cos(angel * M_PI/180.0);
    matrix[5] = -sin(angel * M_PI/180.0);

    matrix[6] = 0.0;
    matrix[7] = sin(angel * M_PI/180.0);
    matrix[8] = cos(angel * M_PI/180.0);
}

void RotMatrixY(double angel, double* matrix)
{
    matrix[0] = cos(angel * M_PI/180.0);
    matrix[1] = 0.0;
    matrix[2] = sin(angel * M_PI/180.0);

    matrix[3] = 0.0;
    matrix[4] = 1.0;
    matrix[5] = 0.0;

    matrix[6] = -sin(angel * M_PI/180.0);
    matrix[7] = 0.0;
    matrix[8] = cos(angel * M_PI/180.0);
}

void RotMatrixZ(double angel, double* matrix)
{
    matrix[0] = cos(angel * M_PI/180.0);
    matrix[1] = -sin(angel * M_PI/180.0);
    matrix[2] = 0.0;

    matrix[3] = sin(angel * M_PI/180.0);
    matrix[4] = cos(angel * M_PI/180.0);
    matrix[5] = 0.0;

    matrix[6] = 0.0;
    matrix[7] = 0.0;
    matrix[8] = 1.0;
}

void dcm2quaternion(const double* dcm, double* quat)
{
    double tol = 1e-5;
    quat[0] = sqrt(0.25*(dcm[0]+dcm[4]+dcm[8]+1));
    if(quat[0] > tol)
    {
        double s = 0.25/quat[0];
        quat[1] = s * (dcm[5] - dcm[7]);
        quat[2] = s * (dcm[6] - dcm[2]);
        quat[3] = s * (dcm[1] - dcm[3]);
    }
    else
    {
        quat[1] = sqrt(0.25 * (1.0 + dcm[0] - dcm[4] - dcm[8]));
        double s = 0.25/quat[1];
        quat[0] = s * (dcm[5] - dcm[7]);
        quat[2] = s * (dcm[1] + dcm[3]);
        quat[3] = s * (dcm[2] + dcm[6]);

    }
}

void MatVectMult(const double* matrix, const double* vector, double* dest)
{
    for(int i = 0; i < 3; i++)
    {
        dest[i] = 0.0;
        for(int j = 0; j < 3; j++)
        {
            dest[i] += matrix[3*i+j] * vector[j];
        }
    }
}

void PrintMatrix(const double* matrix)
{
    for(int j = 0; j < 3; j++)
    {
        for(int i = 0; i < 3; i++)
            printf("%lf\t", matrix[3*j+i]);
        printf("\n");
    }
}


double determinant3x3(double* matrix)
{
    double a, b, c, d, e, f, g, h, i, det;
    a = matrix[0];
    b = matrix[1];
    c = matrix[2];
    d = matrix[3];
    e = matrix[4];
    f = matrix[5];
    g = matrix[6];
    h = matrix[7];
    i = matrix[8];
    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    return det;
}

void quat2Euler(const double* quat, double* ra, double* dec, double* phi)
{
    // 3-2-1//+/-
    *ra = (180.0/M_PI)*atan2(2.0*(quat[0]*quat[1] + quat[2]*quat[3]),
                             1.0 - 2.0*(quat[1]*quat[1] + quat[2]*quat[2]));
    *dec = (180.0/M_PI)*atan2(2.0*(quat[0]*quat[3] + quat[1]*quat[2]),
                              1.0 - 2.0*(quat[2]*quat[2] + quat[3]*quat[3]));
    *phi = (180.0/M_PI)*asin(2.0*(quat[0]*quat[2] - quat[1]*quat[3]));
}
double quat_vector_norm(const double* quat)
{
    return sqrt(quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);
}
