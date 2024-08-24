#include "matrix.h"
#define RADEC

double randfrom(double min, double max);
void TRIAD(const double* sun_body,
           const double* magnetic_field_body,
           const double* sun_inertial,
           const double* magnetic_field_inertial,
           double* Rbi);
void radec2xyzarr(double ra, double dec, double* xyz);
void xyz2radec(double x, double y, double z, double *ra, double *dec);
int main()
{
    //the sun vector respect to the body frame
    double sun_body[3] = {0.8273, 0.5541, -0.09020};
    //the magnetic field vector respect to the body frame
    double magnetic_field_body[3] = {-0.8285, 0.5522, -0.0955};
    //the sun vector respect to the inertial frame
    double sun_inertial[3] = {-0.1517, -0.9669, 0.2050};
    //the magnetic field vector respect to the inertial frame
    double magnetic_field_inertial[3] = {-0.8393, 0.4494, -0.3044};

    double Rbi[9];

    TRIAD(sun_body, magnetic_field_body,
          sun_inertial, magnetic_field_inertial, Rbi);

    printf("Rbi:\n");
    PrintMatrix(Rbi);
    printf("\n");
    /*************************************************/
    printf("*********************************\n");
    srand(time(NULL));
    for(int i = 0; i < 3; i++)
    {
        sun_body[i] = randfrom(-1.0, 1.0);
        magnetic_field_body[i] = randfrom(-1.0, 1.0);
    }
    vector_normalazer(sun_body);
    vector_normalazer(magnetic_field_body);

    RotMatrixX(45, Rbi);

    MatVectMult(Rbi, sun_body, sun_inertial);
    MatVectMult(Rbi, magnetic_field_body, magnetic_field_inertial);

    TRIAD(sun_body, magnetic_field_body,
          sun_inertial, magnetic_field_inertial, Rbi);

    printf("Rbi:\n");
    PrintMatrix(Rbi);
    printf("\n");
    /*************************************************/
    printf("*********************************\n");
#ifndef RADEC
    double ra, dec;
#endif

#ifdef RADEC
    double raSB = 90.0;
    double decSB = 5.0;
    radec2xyzarr(raSB*M_PI/180.0, decSB*M_PI/180.0, sun_body);
#else
    sun_body[0] = 0.0;
    sun_body[1] = 0.0;
    sun_body[2] = 1.0;
    vector_normalazer(sun_body);
    xyz2radec(sun_body[0], sun_body[1], sun_body[2], &ra, &dec);
    printf("sun_body = %lf,%lf\n", ra*180.0/M_PI, dec*180.0/M_PI);
#endif

#ifdef RADEC
    double raMFB = 90.0;
    double decMFB = -5.0;
    radec2xyzarr(raMFB*M_PI/180.0, decMFB*M_PI/180.0, magnetic_field_body);
#else
    magnetic_field_body[0] = 1.0;
    magnetic_field_body[1] = 0.0;
    magnetic_field_body[2] = 0.0;
    vector_normalazer(magnetic_field_body);
    xyz2radec(magnetic_field_body[0], magnetic_field_body[1], magnetic_field_body[2], &ra, &dec);
    printf("magnetic_field_body = %lf,%lf\n", ra*180.0/M_PI, dec*180.0/M_PI);
#endif

#ifdef RADEC
    double raSI = 0.0;
    double decSI = -5.0;
    radec2xyzarr(raSI*M_PI/180.0, decSI*M_PI/180.0, sun_inertial);
#else
    sun_inertial[0] = 1.0;
    sun_inertial[1] = 0.0;
    sun_inertial[2] = 0.0;
    vector_normalazer(sun_inertial);
    xyz2radec(sun_inertial[0], sun_inertial[1], sun_inertial[2], &ra, &dec);
    printf("sun_inertial = %lf,%lf\n", ra*180.0/M_PI, dec*180.0/M_PI);
#endif

#ifdef RADEC
    double raMFI = 0.0;
    double decMFI = 5.0;
    radec2xyzarr(raMFI*M_PI/180.0, decMFI*M_PI/180.0, magnetic_field_inertial);
#else
    magnetic_field_inertial[0] = 0.0;
    magnetic_field_inertial[1] = 0.0;
    magnetic_field_inertial[2] = 1.0;
    vector_normalazer(magnetic_field_inertial);
    xyz2radec(magnetic_field_inertial[0], magnetic_field_inertial[1], magnetic_field_inertial[2], &ra, &dec);
    printf("magnetic_field_inertial = %lf,%lf\n", ra*180.0/M_PI, dec*180.0/M_PI);
#endif

    TRIAD(sun_body, magnetic_field_body,
          sun_inertial, magnetic_field_inertial, Rbi);

    double Qbi[4];
    dcm2quaternion(Rbi, Qbi);

    double raBi, decBi, phiBi;
    quat2Euler(Qbi, &raBi, &decBi, &phiBi);

    printf("Rbi:\n");
    PrintMatrix(Rbi);
    printf("\n");
    printf("det(Rbi):%lf\n", determinant3x3(Rbi));

    printf("Qbi:\n");
    for(int i = 0; i < 4; i++)
        printf("%lf\t", Qbi[i]);
    printf("\n");

    double angel = (180.0/M_PI)*acos(Qbi[0])*2.0;
    double normQuat = quat_vector_norm(Qbi);
    double orth_x = Qbi[1]/asin(normQuat)*/*(M_PI/180.0)**/2.0;
    double orth_y = Qbi[2]/asin(normQuat)*/*(M_PI/180.0)**/2.0;
    double orth_z = Qbi[3]/asin(normQuat)*/*(M_PI/180.0)**/2.0;



//    double normQuat = quat_vector_norm(Qbi);
//    double angel_p = (180.0/M_PI)*acos(normQuat)*2.0;
//    double ex = Qbi[1]/asin(normQuat)*2.0;
//    double ey = Qbi[2]/asin(normQuat)*2.0;
//    double ez = Qbi[3]/asin(normQuat)*2.0;

//    printf("ra = %lf, dec = %lf, phi = %lf\n", raBi, decBi, phiBi);

    printf("angel = %lf, orth_x = %lf, orth_y = %lf, orth_z = %lf\n", angel, orth_x, orth_y, orth_z);
//    printf("angel_p = %lf, ex = %lf, ey = %lf, ez = %lf\n", angel_p, ex, ey, ez);

    printf("Hello World!\n");
    return 0;
}

double randfrom(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void TRIAD(const double* sun_body,
           const double* magnetic_field_body,
           const double* sun_inertial,
           const double* magnetic_field_inertial,
           double* Rbi)
{
    double t1b[3];
    memcpy(t1b, sun_body, sizeof(double) * 3);

    double t1i[3];
    memcpy(t1i, sun_inertial, sizeof(double) * 3);

    double tb2[3];
    cross_product(sun_body, magnetic_field_body, tb2);
    vector_normalazer(tb2);
    double tb3[3];
    cross_product(t1b, tb2, tb3);

    double ti2[3];
    cross_product(sun_inertial, magnetic_field_inertial, ti2);
    vector_normalazer(ti2);
    double ti3[3];
    cross_product(t1i, ti2, ti3);

    double Rbt[9], Rti[9];

    create_matrix_from_vectors(t1b, tb2, tb3, Rbt);
    create_matrix_from_vectors(t1i, ti2, ti3, Rti);

    transponse_matrix(Rti);

    matrix_multiplication(Rbt, Rti, Rbi);
}

void radec2xyzarr(double ra, double dec, double* xyz) {
    double cosdec = cos(dec);
    xyz[0] = cosdec * cos(ra);
    xyz[1] = cosdec * sin(ra);
    xyz[2] = sin(dec);
}

double z2dec(double z) {
    return asin(z);
}

double xy2ra(double x, double y) {
    double a = atan2(y, x);
    if (a < 0)
        a += 2.0 * M_PI;
    return a;
}

void xyz2radec(double x, double y, double z, double *ra, double *dec) {
    if (ra)
        *ra = xy2ra(x, y);
    if (dec)
        *dec = z2dec(z);
}
