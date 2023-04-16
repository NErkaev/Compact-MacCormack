#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <algorithm>
#define M 2000 // кол-во точек по пространству
long N = 0;    // итерации по времени   in main()
long int i, j, it, file_read = 0,
                   savepoints = 10000;
char del_file_L; // интенсивность потери атмосферы

using namespace std;
double r0 = 1., h = 1 / M;
double geo = 0.0, geo_pow = 1.0, pi = 3.1415926, T00 = 240.4,
       kT0 = 1.1 / 2.0, T0d = T00 * kT0;
double T_pl = 1.0 / kT0;
double Rs = 1., mh = 1.67e-24, d_s = 1.0 * 1.5e13, Grav = 6.67e-8, k_B = 1.38e-16;
double M_p = 1.1 * 0.5045 * 5.9676e27, M_st = 1.0 * 1.98892e33, R0d = 2.4119 * Rs * 6.371e8;
double Cs = sqrt(k_B * T0d / mh), beta_0 = Grav * M_p / (R0d * Cs * Cs);
double beta_s = Grav * M_st / (R0d * Cs * Cs);
double density_0 = 3.960 / k_B / T00 * exp(-1.0 * kT0 * beta_0 * Rs * (1.0 - 1.0 / Rs));
/* density_0= 500/k_B/T00*exp(-1.0*kT0*beta0*(1.0-1.0/Rs));*/
double a_T = 118348.0 / T0d, b_T = log(7.5 * R0d * pow(10.0, -19) * density_0 * 1.0 / (mh * Cs * Cs * Cs));
double d_sp = d_s / R0d, L_1 = d_sp * pow(M_p / M_st / 3, 1.0 / 3.0), difc = 1.0 / (pow(10.0, -15) * density_0 * R0d);
double XUV = 8.67e29 / (4.0 * pi * d_s * d_s), XUV1 = 1.27e+30 / (4.0 * pi * d_s * d_s),
       ai1 = 2.68e-5 * R0d * XUV / (Cs * 400 * 4.62);
/*  ai= 4.76e-5*R0d*XUV/(Cs*400*4.62);*/
double ai = 0.59e-7 * R0d * XUV / Cs, Arec = 2.7e-13 * R0d / Cs * density_0 * 1.0 * pow(10000 / T0d, 0.9);
double Arec1 = 8.0e-33 * pow(300 / T0d, 0.6) * R0d / Cs * density_0 * density_0 * 1.0;
double Arec2 = 0.0 * 2.3e-8 * pow(300 / T0d, 0.4) * density_0 * 1.0;
double Rm = 25, T_b = 1.0 / kT0;
double a = density_0 * R0d * 1.89e-18;
double qN = 1.89e-18 * XUV * R0d * 0.15 / (Cs * Cs * Cs * mh);
double aX = 1 * density_0 * R0d * 0.8895e-21;

double qNX = 0.8895e-21 * XUV1 * R0d * 0.32 / (Cs * Cs * Cs * mh);

double Gamma = 5. / 3.;
int bt = beta_0;
int nfile = bt * 10 + M / 1000;

double r[M + 1], ksi[M + 1], X[M + 1], L;
double D_F[4][M + 1], D_B[4][M + 1];
double S0[M + 1], SX0[M + 1], X0[M + 1], G0[M + 1], E0[M + 1], W0[M + 1], tt[M + 1],
    S1[M + 1], SX1[M + 1], G1[M + 1], E1[M + 1], W1[M + 1], X1[M + 1],
    S2[M + 1], SX2[M + 1], G2[M + 1], E2[M + 1], W2[M + 1], X2[M + 1],
    S3[M + 1], SX3[M + 1], G3[M + 1], E3[M + 1], W3[M + 1], X3[M + 1],
    // n - density, v - velocity, T - temperature P - pressure
    n[M + 1], nX[M + 1], S[M + 1], SX[M + 1], E[M + 1], T[M + 1], P[M + 1],
    q[M + 1], v[M + 1], v_[M + 1], qX[M + 1],
    T0[M + 1], n0[M + 1], nX0[M + 1], P0[M + 1], v0[M + 1],
    T1[M + 1], n1[M + 1], nX1[M + 1], P1[M + 1],
    T2[M + 1], n2[M + 1], nX2[M + 1], P2[M + 1],
    T3[M + 1], n3[M + 1], nX3[M + 1], P3[M + 1];
double h1[4][M + 1], h2[4][M + 1], h3[4][M + 1], h4[4][M + 1];

char fname[50], fname_time[50];
FILE *fin_n, *fin_X, *fin_r, *fin_v, *fin_q, *fin_P, *fin_T, *fin_L;
FILE *fout_n, *fout_X, *fout_r, *fout_v, *fout_q, *fout_P, *fout_T, *fout_L;

double Force(double r)
{
    double ss;
    ss = (-beta_0 / r / r + beta_s * (3 * d_sp * d_sp * r - 3 * d_sp * r * r + r * r * r) / (d_sp * d_sp * d_sp * (d_sp - r) * (d_sp - r)) +
          beta_0 * r / (d_sp * d_sp * d_sp));
    /*   ss=-beta_0/r/r + 3*beta_s*r/(d_sp*d_sp*d_sp);*/

    /* ss=-beta_0/r/r;*/
    return ss;
}

void Boundary_stensils_4_4_backward(double *S, double *G, double *W, double *X)
{
    D_B[0][0] = (-(37. / 18) * S[0] + (35. / 9) * S[1] - (17. / 6) * S[2] +
                 (11. / 9) * S[3] - (2. / 9) * S[4]) /
                (ksi[1] - ksi[0]); //(r[1] - r[0]);
    D_B[1][0] = (-(37. / 18) * G[0] + (35. / 9) * G[1] - (17. / 6) * G[2] +
                 (11. / 9) * G[3] - (2. / 9) * G[4]) /
                (ksi[1] - ksi[0]); //(r[1] - r[0]);
    D_B[2][0] = (-(37. / 18) * W[0] + (35. / 9) * W[1] - (17. / 6) * W[2] +
                 (11. / 9) * W[3] - (2. / 9) * W[4]) /
                (ksi[1] - ksi[0]); //(r[1] - r[0]);
    D_B[3][0] = (-(37. / 18) * X[0] + (35. / 9) * X[1] - (17. / 6) * X[2] +
                 (11. / 9) * X[3] - (2. / 9) * X[4]) /
                (ksi[1] - ksi[0]); //(r[1] - r[0]);

    D_B[0][M] = ((19. / 9) * S[M] - (37. / 9) * S[M - 1] + (19. / 6) * S[M - 2] - (13. / 9) * S[M - 3] + (5. / 18) * S[M - 4]) / (ksi[M] - ksi[M - 1]); //(r[M]-r[M-1]);
    D_B[1][M] = ((19. / 9) * G[M] - (37. / 9) * G[M - 1] + (19. / 6) * G[M - 2] - (13. / 9) * G[M - 3] + (5. / 18) * G[M - 4]) / (ksi[M] - ksi[M - 1]); //(r[M]-r[M-1]);
    D_B[2][M] = ((19. / 9) * W[M] - (37. / 9) * W[M - 1] + (19. / 6) * W[M - 2] - (13. / 9) * W[M - 3] +
                 (5. / 18) * W[M - 4]) /
                (ksi[M] - ksi[M - 1]); //(r[M]-r[M-1]);
    D_B[3][M] = ((19. / 9) * X[M] - (37. / 9) * X[M - 1] + (19. / 6) * X[M - 2] - (13. / 9) * X[M - 3] +
                 (5. / 18) * X[M - 4]) /
                (ksi[M] - ksi[M - 1]); //(r[M]-r[M-1]);
}

void Boundary_stensils_4_4_forward(double *S, double *G, double *W, double *SX)
{
    D_F[0][0] = (-(19. / 9) * S[0] + (37. / 9) * S[1] - (19. / 6) * S[2] + (13. / 9) * S[3] - (5. / 18) * S[4]) / (ksi[1] - ksi[0]); //(r[1]-r[0]);
    D_F[1][0] = (-(19. / 9) * G[0] + (37. / 9) * G[1] - (19. / 6) * G[2] + (13. / 9) * G[3] - (5. / 18) * G[4]) / (ksi[1] - ksi[0]); //(r[1]-r[0]);
    D_F[2][0] = (-(19. / 9) * W[0] + (37. / 9) * W[1] - (19. / 6) * W[2] + (13. / 9) * W[3] -
                 (5. / 18) * W[4]) /
                (ksi[1] - ksi[0]); //(r[1]-r[0]);
    D_F[3][0] = (-(19. / 9) * X[0] + (37. / 9) * X[1] - (19. / 6) * X[2] + (13. / 9) * X[3] -
                 (5. / 18) * X[4]) /
                (ksi[1] - ksi[0]); //(r[1]-r[0]);

    D_F[0][M] = ((37. / 18) * S[M] - (35. / 9) * S[M - 1] + (17. / 6) * S[M - 2] - (11. / 9) * S[M - 3] + (2. / 9) * S[M - 4]) / (ksi[M] - ksi[M - 1]); //(r[M] - r[M - 1]);
    D_F[1][M] = ((37. / 18) * G[M] - (35. / 9) * G[M - 1] + (17. / 6) * G[M - 2] - (11. / 9) * G[M - 3] + (2. / 9) * G[M - 4]) / (ksi[M] - ksi[M - 1]); //(r[M] - r[M - 1]);
    D_F[2][M] = ((37. / 18) * W[M] - (35. / 9) * W[M - 1] + (17. / 6) * W[M - 2] -
                 (11. / 9) * W[M - 3] + (2. / 9) * W[M - 4]) /
                (ksi[M] - ksi[M - 1]); //(r[M] - r[M - 1]);
    D_F[3][M] = ((37. / 18) * SX[M] - (35. / 9) * SX[M - 1] + (17. / 6) * X[M - 2] -
                 (11. / 9) * SX[M - 3] + (2. / 9) * SX[M - 4]) /
                (ksi[M] - ksi[M - 1]); //(r[M] - r[M - 1]);
}

void D_backward_4_2(double *S, double *G, double *W, double *SX)
{
    // for D_B(i min)
    /*for (long it = 1; it <= M; it++) { //(r[it] - r[it - 1])
        D_B[0][it] = ((1. / 6.*S[it + 1] + 2. / 3.*S[it] - 5. / 6.*S[it - 1]) / (ksi[it] - ksi[it - 1]) - 1. / 3.*D_B[0][it - 1]) * 3. / 2.;
        D_B[1][it] = ((1. / 6.*G[it + 1] - 2. / 3.*G[it] - 5. / 6.*G[it - 1]) / (ksi[it] - ksi[it - 1]) - 1. / 3.*D_B[1][it - 1]) * 3. / 2.;
        D_B[2][it] = ((1. / 6.*W[it + 1] - 2. / 3.*W[it] - 5. / 6.*W[it - 1]) / (ksi[it] - ksi[it - 1]) - 1. / 3.*D_B[2][it - 1]) * 3. / 2.;
    }*/

    // Компактная схема 4_2
    double a_ = 0.5 * (1 - 1. / sqrt(3));
    double G1_S, G1_G, G1_W, G1_X;
    G1_S = -1.0 / 12.0 * (25 + 17. / sqrt(3)) * S[0] + (4 + 25.0 / 6. / sqrt(3)) * S[1] - 3.0 * (1 + sqrt(3) / 2.) * S[2] + 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * S[3] - 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * S[4];
    D_B[0][0] = G1_S / (ksi[1] - ksi[0]);

    G1_G = -1.0 / 12.0 * (25 + 17. / sqrt(3)) * G[0] + (4 + 25.0 / 6. / sqrt(3)) * G[1] - 3.0 * (1 + sqrt(3) / 2.) * G[2] + 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * G[3] - 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * G[4];
    D_B[1][0] = G1_G / (ksi[1] - ksi[0]);

    G1_W = -1.0 / 12.0 * (25 + 17. / sqrt(3)) * W[0] + (4 + 25.0 / 6. / sqrt(3)) * W[1] - 3.0 * (1 + sqrt(3) / 2.) * W[2] + 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * W[3] - 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * W[4];
    D_B[2][0] = G1_W / (ksi[1] - ksi[0]);

    G1_X = -1.0 / 12.0 * (25 + 17. / sqrt(3)) * SX[0] + (4 + 25.0 / 6. / sqrt(3)) * SX[1] - 3.0 * (1 + sqrt(3) / 2.) * SX[2] + 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * SX[3] - 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * SX[4];
    D_B[3][0] = G1_X / (ksi[1] - ksi[0]);

    for (it = 1; it <= M; it++)
    {
        D_B[0][it] = (S[it] - S[it - 1]) / (ksi[it] - ksi[it - 1]) * 1 / (1 - a_) - D_B[0][it - 1] * a_ / (1 - a_);
        D_B[1][it] = (G[it] - G[it - 1]) / (ksi[it] - ksi[it - 1]) * 1 / (1 - a_) - D_B[1][it - 1] * a_ / (1 - a_);
        D_B[2][it] = (W[it] - W[it - 1]) / (ksi[it] - ksi[it - 1]) * 1 / (1 - a_) -
                     D_B[2][it - 1] * a_ / (1 - a_);
        D_B[3][it] = (SX[it] - SX[it - 1]) / (ksi[it] - ksi[it - 1]) * 1 / (1 - a_) -
                     D_B[3][it - 1] * a_ / (1 - a_);
    }
}

void D_forward_4_2(double *S, double *G, double *W, double *SX)
{
    // for D_F(i max)
    /*for (long k = M - 1; k >= 0; k--) { //(r[k + 1] - r[k])
        D_F[0][k] = ((5. / 6.*S[k + 1] - 2. / 3.*S[k] - 1. / 6.*S[k - 1]) / (ksi[k + 1] - ksi[k]) - 1. / 3.*D_F[0][k + 1]) * 3. / 2.;
        D_F[1][k] = ((5. / 6.*G[k + 1] - 2. / 3.*G[k] - 1. / 6.*G[k - 1]) / (ksi[k + 1] - ksi[k]) - 1. / 3.*D_F[1][k + 1]) * 3. / 2.;
        D_F[2][k] = ((5. / 6.*W[k + 1] - 2. / 3.*W[k] - 1. / 6.*W[k - 1]) / (ksi[k + 1] - ksi[k]) - 1. / 3.*D_F[2][k + 1]) * 3. / 2.;
    }*/

    // Компактная схема 4_2
    double a_ = 0.5 * (1 - 1 / sqrt(3));
    double G1_S, G1_G, G1_W, G1_X;
    G1_S = 1.0 / 12.0 * (25 + 17. / sqrt(3)) * S[M] - (4 + 25.0 / 6. / sqrt(3)) * S[M - 1] + 3.0 * (1 + sqrt(3) / 2) * S[M - 2] - 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * S[M - 3] + 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * S[M - 4];
    D_F[0][M] = G1_S / (ksi[M] - ksi[M - 1]);

    G1_G = 1.0 / 12.0 * (25 + 17 / sqrt(3)) * G[M] - (4 + 25.0 / 6 / sqrt(3)) * G[M - 1] + 3.0 * (1 + sqrt(3) / 2) * G[M - 2] - 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * G[M - 3] + 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * G[M - 4];
    D_F[1][M] = G1_G / (ksi[M] - ksi[M - 1]);

    G1_W = 1.0 / 12.0 * (25 + 17 / sqrt(3)) * W[M] - (4 + 25.0 / 6 / sqrt(3)) * W[M - 1] + 3.0 * (1 + sqrt(3) / 2) * W[M - 2] - 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * W[M - 3] + 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * W[M - 4];
    D_F[2][M] = G1_W / (ksi[M] - ksi[M - 1]);

    G1_X = 1.0 / 12.0 * (25 + 17 / sqrt(3)) * SX[M] - (4 + 25.0 / 6 / sqrt(3)) * SX[M - 1] + 3.0 * (1 + sqrt(3) / 2) * SX[M - 2] - 1.0 / 3.0 * (4 + 13.0 / 2.0 / sqrt(3)) * SX[M - 3] + 1 / 4.0 * (1 + 5.0 / 3.0 / sqrt(3)) * SX[M - 4];
    D_F[3][M] = G1_X / (ksi[M] - ksi[M - 1]);

    for (i = M - 1; i >= 0; i--)
    {
        D_F[0][i] = (S[i + 1] - S[i]) / (ksi[i + 1] - ksi[i]) * 1. / (1 - a_) - D_F[0][i + 1] * a_ / (1 - a_);
        D_F[1][i] = (G[i + 1] - G[i]) / (ksi[i + 1] - ksi[i]) * 1. / (1 - a_) - D_F[1][i + 1] * a_ / (1 - a_);
        D_F[2][i] = (W[i + 1] - W[i]) / (ksi[i + 1] - ksi[i]) * 1. / (1 - a_) -
                    D_F[2][i + 1] * a_ / (1 - a_);
        D_F[3][i] = (SX[i + 1] - SX[i]) / (ksi[i + 1] - ksi[i]) * 1. / (1 - a_) -
                    D_F[3][i + 1] * a_ / (1 - a_);
    }
}

void For_write(double v[], double n[], double P[], double T[], double q[], double X[], long j)
{
    // Non-uniform grid/
    sprintf(fname, "init/Density%d.dat", nfile + file_read + j);
    fout_n = fopen(fname, "w");
    sprintf(fname, "init/Velocity%d.dat", nfile + file_read + j);
    fout_v = fopen(fname, "w");
    sprintf(fname, "init/Heating%d.dat", nfile + file_read + j);
    fout_q = fopen(fname, "w");
    sprintf(fname, "init/Pressure%d.dat", nfile + file_read + j);
    fout_P = fopen(fname, "w");
    sprintf(fname, "init/Temperature%d.dat", nfile + file_read + j);
    fout_T = fopen(fname, "w");
    sprintf(fname, "init/X%d.dat", nfile + file_read + j);
    fout_X = fopen(fname, "w");
    sprintf(fname, "init/r%d.dat", nfile + file_read + j);
    fout_r = fopen(fname, "w");
    for (i = 0; i <= M; i++)
        fprintf(fout_n, "%.10e\n", n[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_v, "%.10e\n", v[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_q, "%.10e\n", q[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_P, "%.10e\n", P[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_T, "%.10e\n", T[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_X, "%.10e\n", X[i]);
    for (i = 0; i <= M; i++)
        fprintf(fout_r, "%.10e\n", r[i]);

    fclose(fout_n);
    fclose(fout_v);
    fclose(fout_q);
    fclose(fout_P);
    fclose(fout_T);
    fclose(fout_X);
    fclose(fout_r);
}

int For_read(double v[], double n[], double P[], double T[], double q[], double time[])
{ // прописать путь до файла
    long int it;
    double tmpvar;
    // Non-uniform grid/
    sprintf(fname, "init%d/Density%d.dat", nfile * 1000, nfile + file_read);
    fin_n = fopen(fname, "r");
    if (fin_n == NULL)
    {
        printf("File init%d/Density%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_n, "%lf", &tmpvar) != EOF)
    {
        n[it] = tmpvar;
        it++;
    }
    fclose(fin_n);

    sprintf(fname, "init%d/X%d.dat", nfile * 1000, nfile + file_read);
    fin_X = fopen(fname, "r");
    if (fin_X == NULL)
    {
        printf("File init%d/X%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_X, "%lf", &tmpvar) != EOF)
    {
        X[it] = tmpvar;
        it++;
    }
    fclose(fin_X);

    sprintf(fname, "init%d/Velocity%d.dat", nfile * 1000, nfile + file_read);
    fin_v = fopen(fname, "r");
    if (fin_v == NULL)
    {
        printf("File init%d/Velocity%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_v, "%lf", &tmpvar) != EOF)
    {
        v[it] = tmpvar;
        it++;
    }
    fclose(fin_v);

    sprintf(fname, "init%d/Heating%d.dat", nfile * 1000, nfile + file_read);
    fin_q = fopen(fname, "r");
    if (fin_q == NULL)
    {
        printf("File init%d/Heating%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_q, "%lf", &tmpvar) != EOF)
    {
        q[it] = tmpvar;
        it++;
    }
    fclose(fin_q);

    sprintf(fname, "init%d/Pressure%d.dat", nfile * 1000, nfile + file_read);
    fin_P = fopen(fname, "r");
    if (fin_P == NULL)
    {
        printf("File init%d/Pressure%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_P, "%lf", &tmpvar) != EOF)
    {
        P[it] = tmpvar;
        it++;
    }
    fclose(fin_P);
    sprintf(fname, "init%d/Temperature%d.dat", nfile * 1000, nfile + file_read);
    fin_T = fopen(fname, "r");
    if (fin_T == NULL)
    {
        printf("File init%d/Temperature%d.dat doesn't exist\n", nfile * 1000, nfile + file_read);
        return -1;
    }
    it = 0;
    while (fscanf(fin_T, "%lf", &tmpvar) != EOF)
    {
        T[it] = tmpvar;
        it++;
    }
    fclose(fin_T);

    sprintf(fname, "init%d/Time%d.dat", nfile * 1000, nfile);
    fin_T = fopen(fname, "r");
    if (fin_T == NULL)
    {
        printf("File init%d/Time%d.dat doesn't exist\n", nfile * 1000, nfile);
        return -1;
    }
    it = 0;
    while (fscanf(fin_T, "%lf", &tmpvar) != EOF)
    {
        time[0] = tmpvar;
    }
    fclose(fin_T);
    return 1;
}

void RK4_for_system(double n0_[], double v0_[], double T0_[])
{
    double alpha[4] = {0, 1. / 2, 1. / 2, 1.},
           beta[4] = {1. / 6, 1. / 3, 1. / 3, 1. / 6};
    double c[M + 1]; // скорость звука
    double K = 1.0, tau = 0;
    double *time = new double[N + 2];

    double max_v, max_c;

    // переобозначения
    if (!file_read)
        for (i = 0; i <= M; i++)
        {
            n0[i] = pow(r[i], 2) * n0_[i];
            T0[i] = 1. / T_b;
            X0[i] = 0;
            v0[i] = v0_[i];
            P0[i] = n0[i] * T0[i];
            S0[i] = n0[i] * v0[i];
            G0[i] = pow(S0[i], 2) / n0[i] + P0[i];
            E0[i] = pow(S0[i], 2) / 2. / n0[i] + P0[i] / (Gamma - 1);
            W0[i] = S0[i] / n0[i] * (Gamma * E0[i] - (Gamma - 1) * pow(S0[i], 2) / 2. / n0[i]);
            time[0] = 0.0;
        }
    else if (For_read(v0, n0, P0, T0, q, time) < 0)
    {
        printf("ERROR!");
        delete[] time;
        exit(1);
    }
    printf("L0=%.14lf\n", v0[M] * n0[M]);

    n0[0] = 1.;
    X0[0] = 0, T0[0] = T_b;
    P0[0] = n0[0] * T0[0];
    v0[0] = v0[1];
    for (i = 0; i <= M; i++)
    {
        S0[i] = v0[i] * n0[i];
        T0[i] = P0[i] / n0[i] / (1 + X0[i]);
        E0[i] = pow(S0[i], 2) / 2. / n0[i] + P0[i] / (Gamma - 1);
    }
    for (j = 0; j <= N; j++)
    {
        for (i = 0; i <= M; i++)
        {
            P0[i] = (E0[i] - pow(S0[i], 2) / 2. / n0[i]) * (Gamma - 1);
            T0[i] = P0[i] / n0[i] / (1 + X0[i]);
            SX0[i] = S0[i] * X0[i];
            nX0[i] = n0[i] * X0[i];
            v0[i] = S0[i] / n0[i];
            G0[i] = (Gamma - 1) * E0[i] + (3 - Gamma) / (2 * n0[i]) * pow(S0[i], 2);
            W0[i] = S0[i] / n0[i] * (Gamma * E0[i] - (Gamma - 1) * pow(S0[i], 2) / 2. / n0[i]);
        }

        for (i = 1; i <= M; i++)
        {
            c[i] = sqrt(Gamma * P0[i] / n0[i]);
            v_[i] = (sqrt(v0[i] * v0[i]) + c[i]) / (r[i] - r[i - 1]);
        }
        max_v = *max_element(v_, v_ + M);
        tau = K * 1.0 / (max_v); // max(|v|+c)
        time[j + 1] = time[j] + tau;

        q[M] = qN;
        // Heating_rate(n0, q);
        tt[M] = 0;
        for (i = 0; i < M; i++)
        {
            tt[M - i - 1] = tt[M - i] + 0.5 * a * (r[M - i] - r[M - i - 1]) * (n0[M - i] * (1. - X0[M - i]) / r[M - i] / r[M - i] +
                                                                               /* tt[M-i-1]= tt[M-i] + 0.5*a*(r[M-i]-r[M-i-1])*(n0[M-i]*(1.)/r[M-i]/r[M-i]+ */
                                                                               n0[M - i - 1] * (1. - X0[M - i - 1]) / r[M - i - 1] / r[M - i - 1]);
            q[M - i - 1] = q[M - 1] * exp(-tt[M - i - 1]) / (1. + geo * pow(tt[M - i - 1], geo_pow));
        }
        qX[M] = qNX;
        tt[M] = 0;
        for (i = 0; i < M; i++)
        {
            tt[M - i - 1] = tt[M - i] + aX * (r[M - i] - r[M - i - 1]) * (n0[M - i] * (1. - X0[M - i]) / r[M - i] / r[M - i] + n0[M - i - 1] * (1. - X0[M - i - 1]) / r[M - i - 1] / r[M - i - 1]);

            qX[M - i - 1] = qX[M] * exp(-tt[M - i - 1]) / (1. + geo * pow(tt[M - i - 1], geo_pow));
        }

        // Boundary_stensils_4_4_forward(S0, G0, W0);
        D_forward_4_2(S0, G0, W0, SX0);
        // FIRST STAGE
        for (i = 0; i <= M; i++)
        {
            /*----->*/
            h1[0][i] = -tau * D_F[0][i] / r[i] / log(Rm);
            h1[1][i] = -tau * D_F[1][i] / r[i] / log(Rm) + tau * n0[i] * Force(r[i]) + tau * 2 * (Gamma - 1) / r[i] * (E0[i] - S0[i] * S0[i] / 2 / n0[i]);
            h1[2][i] = -tau * D_F[2][i] / r[i] / log(Rm) + tau * S0[i] * Force(r[i]) + tau * n0[i] * (1 - X0[i]) * q[i] - tau * X0[i] * (1 - X0[i]) * n0[i] * n0[i] * exp(b_T - a_T / T0[i]) / r[i] / r[i];
            /*    - tau*S0[i]*beta_0/pow(r[i],2) + tau*n0[i]*(1)*q[i];*/

            h1[3][i] = -1.0 * tau * D_F[3][i] / r[i] / log(Rm) +
                       tau * ai * n0[i] * (1 - X0[i]) * q[i] / qN -
                       1. * tau * Arec * n0[i] * n0[i] * X0[i] * X0[i] / r[i] / r[i] * pow(T0[i], -0.9);

            n1[i] = n0[i] + h1[0][i] * alpha[1];
            S1[i] = S0[i] + h1[1][i] * alpha[1];
            E1[i] = E0[i] + h1[2][i] * alpha[1];
            nX1[i] = nX0[i] + h1[3][i] * alpha[1];
            X1[i] = nX1[i] / n1[i];
            if (X1[i] > 1)
                X1[i] = 1;
            if (X1[i] < 0)
                X1[i] = 0;
            SX1[i] = S1[i] * X1[i];

            P1[i] = (E1[i] - pow(S1[i], 2) / 2. / n1[i]) * (Gamma - 1);
            T1[i] = P1[i] / n1[i] / (1 + X1[i]);
            /*   T1[i]=P1[i]/n1[i]/(1);*/
            G1[i] = (Gamma - 1) * E1[i] + (3 - Gamma) / (2 * n1[i]) * pow(S1[i], 2);
            W1[i] = S1[i] / n1[i] * (Gamma * E1[i] - (Gamma - 1) * pow(S1[i], 2) / (2 * n1[i]));
        }

        // Boundary_stensils_4_4_backward(S1, G1, W1);
        D_backward_4_2(S1, G1, W1, SX1);
        // SECOND STAGE
        for (i = 0; i <= M; i++)
        {
            /*----->*/
            h2[0][i] = -tau * D_B[0][i] / r[i] / log(Rm);
            h2[1][i] = -tau * D_B[1][i] / r[i] / log(Rm) + tau * n1[i] * Force(r[i]) + tau * 2 * (Gamma - 1) / r[i] * (E1[i] - S1[i] * S1[i] / 2 / n1[i]);
            h2[2][i] = -tau * D_B[2][i] / r[i] / log(Rm) + tau * S1[i] * Force(r[i]) + tau * n1[i] * (1 - X1[i]) * q[i] - tau * X1[i] * (1 - X1[i]) * n1[i] * n1[i] * exp(b_T - a_T / T1[i]) / r[i] / r[i];
            /*   - tau*S1[i] *beta_0/pow(r[i],2) + tau*n1[i] * (1)*q[i];*/
            h2[3][i] = -1.0 * tau * D_B[3][i] / r[i] / log(Rm) +
                       tau * ai * n1[i] * (1 - X1[i]) * q[i] / qN -
                       1. * tau * Arec * n1[i] * n1[i] * X1[i] * X1[i] / r[i] / r[i] * pow(T1[i], -0.9);
            n2[i] = n0[i] + h2[0][i] * alpha[2];
            S2[i] = S0[i] + h2[1][i] * alpha[2];
            E2[i] = E0[i] + h2[2][i] * alpha[2];
            nX2[i] = nX0[i] + h2[3][i] * alpha[2];
            X2[i] = nX2[i] / n2[i];
            if (X2[i] > 1)
                X2[i] = 1;
            if (X2[i] < 0)
                X2[i] = 0;
            SX2[i] = S2[i] * X2[i];
            P2[i] = (E2[i] - pow(S2[i], 2) / 2. / n2[i]) * (Gamma - 1);
            T2[i] = P2[i] / n2[i] / (1 + X2[i]);
            /*    T2[i]=P2[i]/n2[i]/(1);*/
            G2[i] = (Gamma - 1) * E2[i] + (3 - Gamma) / (2 * n2[i]) * pow(S2[i], 2);
            W2[i] = S2[i] / n2[i] * (Gamma * E2[i] - (Gamma - 1) * pow(S2[i], 2) / (2 * n2[i]));
        }

        // Boundary_stensils_4_4_forward(S2, G2, W);
        D_forward_4_2(S2, G2, W2, SX2);

        // THIRD STAGE
        for (i = 0; i <= M; i++)
        {
            /*----->*/
            h3[0][i] = -tau * D_F[0][i] / r[i] / log(Rm);
            h3[1][i] = -tau * D_F[1][i] / r[i] / log(Rm) + tau * n2[i] * Force(r[i]) +
                       tau * 2 * (Gamma - 1) / r[i] * (E2[i] - S2[i] * S2[i] / 2 / n2[i]);
            h3[2][i] = -tau * D_F[2][i] / r[i] / log(Rm) + tau * S2[i] * Force(r[i]) + tau * n2[i] * (1 - X2[i]) * q[i] - tau * X2[i] * (1 - X2[i]) * n2[i] * n2[i] * exp(b_T - a_T / T2[i]) / r[i] / r[i];
            /* - tau*S2[i] *beta_0/pow(r[i],2) + tau*n2[i] *(1)* q[i];*/
            h3[3][i] = -1.0 * tau * D_F[3][i] / r[i] / log(Rm) +
                       tau * ai * n2[i] * (1 - X2[i]) * q[i] / qN -
                       1. * tau * Arec * n2[i] * n2[i] * X2[i] * X2[i] / r[i] / r[i] * pow(T2[i], -0.9);
            n3[i] = n0[i] + h3[0][i] * alpha[3];
            S3[i] = S0[i] + h3[1][i] * alpha[3];
            E3[i] = E0[i] + h3[2][i] * alpha[3];
            nX3[i] = nX0[i] + h3[3][i] * alpha[3];
            X3[i] = nX3[i] / n3[i];
            if (X3[i] > 1)
                X3[i] = 1;
            if (X3[i] < 0)
                X3[i] = 0;
            SX3[i] = S3[i] * X3[i];
            P3[i] = (E3[i] - pow(S3[i], 2) / 2. / n3[i]) * (Gamma - 1);
            T3[i] = P3[i] / n3[i] / (1 + X3[i]);
            /*  T3[i]=P3[i]/n3[i]/(1);*/
            G3[i] = (Gamma - 1) * E3[i] + (3 - Gamma) / (2 * n3[i]) * pow(S3[i], 2);
            W3[i] = S3[i] / n3[i] * (Gamma * E3[i] - (Gamma - 1) * pow(S3[i], 2) / (2 * n3[i]));
        }
        // Boundary_stensils_4_4_backward(S3, G3, W3);
        D_backward_4_2(S3, G3, W3, SX3);

        // FOURTH STAGE
        for (i = 0; i <= M; i++)
        {
            /*----->*/
            h4[0][i] = -tau * D_B[0][i] / r[i] / log(Rm);
            h4[1][i] = -tau * D_B[1][i] / r[i] / log(Rm) + tau * n3[i] * Force(r[i]) +
                       tau * 2 * (Gamma - 1) / r[i] * (E3[i] - S3[i] * S3[i] / 2 / n3[i]);
            h4[2][i] = -tau * D_B[2][i] / r[i] / log(Rm) + tau * S3[i] * Force(r[i]) + tau * n3[i] * (1 - X3[i]) * q[i] - tau * X3[i] * (1 - X3[i]) * n3[i] * n3[i] * exp(b_T - a_T / T3[i]) / r[i] / r[i];
            /*	 - tau*S3[i] *beta_0/pow(r[i],2) + tau*n3[i]*(1)*q[i];*/
            h4[3][i] = -1.0 * tau * D_B[3][i] / r[i] / log(Rm) +
                       tau * ai * n3[i] * (1 - X3[i]) * q[i] / qN -
                       1. * tau * Arec * n3[i] * n3[i] * X3[i] * X3[i] / r[i] / r[i] * pow(T3[i], -0.9);
        }

        for (i = 1; i < M; i++)
        {
            n[i] = n0[i] + beta[0] * h1[0][i] + beta[1] * h2[0][i] +
                   beta[2] * h3[0][i] + beta[3] * h4[0][i];
            S[i] = S0[i] + beta[0] * h1[1][i] + beta[1] * h2[1][i] +
                   beta[2] * h3[1][i] + beta[3] * h4[1][i];

            if (S[i] < 0)
                S[i] = 0;
            E[i] = E0[i] + beta[0] * h1[2][i] + beta[1] * h2[2][i] +
                   beta[2] * h3[2][i] + beta[3] * h4[2][i] + 0.5 * ((T0[i + 1] - T0[i]) / (r[i + 1] + r[i]) * (n0[i] + n0[i + 1]) - (T0[i] - T0[i - 1]) / (r[i] + r[i - 1]) * (n0[i] + n0[i - 1])) / (r[i]) / log(Rm) / log(Rm);
            nX[i] = nX0[i] + beta[0] * h1[3][i] + beta[1] * h2[3][i] +
                    beta[2] * h3[3][i] + beta[3] * h4[3][i];
            X[i] = nX[i] / n[i];
            if (X[i] > 1)
                X[i] = 1;
            if (X[i] < 0)
                X[i] = 0;
            SX[i] = S[i] * X[i];

            P[i] = (E[i] - pow(S[i], 2) / 2. / n[i]) * (Gamma - 1);
            T[i] = P[i] / n[i] / (1 + X[i]);
            /* T[i] = P[i] / n[i]/(1+X[i]);*/
            /* T[i] = P[i] / n[i]/(1);*/
            v[i] = S[i] / n[i];
        }

        n[0] = 1.;
        X[0] = 0;
        T[0] = T_b;
        P[0] = n[0] * T[0];
        v[0] = v[1];
        S[0] = n[0] * v[0];
        E[0] = n[0] * pow(v[0], 2) / 2. + P[0] / (Gamma - 1);

        n[M] = n[M - 1];
        T[M] = T[M - 1];
        S[M] = S[M - 1];
        E[M] = E[M - 1];
        P[M] = (E[M] - pow(S[M], 2) / 2. / n[M]) * (Gamma - 1);
        v[M] = S[M] / n[M];
        X[M] = X[M - 1];
        SX[M] = S[M] * X[M];

        q[M] = qN;
        // Heating_rate(n0, q);
        tt[M] = 0;
        for (i = 0; i < M; i++)
        {
            tt[M - i - 1] = tt[M - i] + 0.5 * a * (r[M - i] - r[M - i - 1]) * (n[M - i] * (1. - X[M - i]) / r[M - i] / r[M - i] +
                                                                               /*  tt[M-i-1]= tt[M-i] + 0.5*a*(r[M-i]-r[M-i-1])*(n[M-i]*(1.)/r[M-i]/r[M-i]+ */
                                                                               n[M - i - 1] * (1. - X[M - i - 1]) / r[M - i - 1] / r[M - i - 1]);
            q[M - i - 1] = q[M] * exp(-tt[M - i - 1]) / (1. + geo * pow(tt[M - i - 1], geo_pow));
        }

        qX[M] = qNX;
        tt[M] = 0;
        for (i = 0; i < M; i++)
        {
            tt[M - i - 1] = tt[M - i] + aX * (r[M - i] - r[M - i - 1]) * (n[M - i] * (1. - X[M - i]) / r[M - i] / r[M - i] + n[M - i - 1] * (1. - X[M - i - 1]) / r[M - i - 1] / r[M - i - 1]);

            qX[M - i - 1] = qX[M] * exp(-tt[M - i - 1]) / (1. + geo * pow(tt[M - i - 1], geo_pow));
        }

        for (i = 0; i <= M; i++)
        {
            E0[i] = E[i];
            S0[i] = S[i];
            n0[i] = n[i];
            v0[i] = S0[i] / n0[i];
            T0[i] = T[i];
            X0[i] = X[i];
        }

        if (j % 100 == 0)
        {
            printf("\n%d\n", j);
            // printf("%e\t%e\t%e  %e\n", n[M], T[M], v[M], q[M]);
        }
        if (j % savepoints == 0 && j != 0)
        { // && j !=0
            printf("\ntime=%g\n", time[j + 1]);
            for (i = M - 5; i <= M; i++)
                printf("%e\t%e\t%e  %e\n", n[i], T[i], v[i], q[i]); //%e\t%e\t%e\t  S1[i],G[i],E_N[i],
            For_write(v, n, P, T, q, X, j);
        }
        if (j % 1000 == 0 && j != 0)
        { //  && j!=0
            printf("\nv[M]=%.10lf\nC[M]=%.10lf\n", v[M], sqrt(Gamma * P[M] / n[M]));
            printf("    L=%.14lf\n", v0[M] * n0[M]);
            // Non-uniform grid/
            sprintf(fname, "init%d/L%d.dat", nfile * 1000, nfile);
            sprintf(fname_time, "init%d/Time%d.dat", nfile * 1000, nfile);

            if (del_file_L == 'y')
            {
                fout_L = fopen(fname, "w");
                fout_T = fopen(fname_time, "w");
                del_file_L = 'n';
            }
            else
            {
                fout_L = fopen(fname, "a");
                fout_T = fopen(fname_time, "a");
            }
            if (fout_L == NULL || fout_T == NULL)
            {
                fout_L = fopen(fname, "w");
                fout_T = fopen(fname_time, "w");
            }
            L = n[M] * v[M];
            fprintf(fout_L, "%.14lf\n", L);
            fclose(fout_L);
            fprintf(fout_T, "%.14lf\n", time[j + 1]);
            fclose(fout_T);
        }
    }
    printf("\n\n");
    printf("Last iteration\n");
    for (i = M - 5; i <= M; i++)
        printf("%e\t%e\t%e  %e\n", n[i], T[i], v[i], q[i]); //%e\t%e\t%e\t  S1[i],G[i],E_N[i],
    delete[] time;
}

int main(int argc, char *argv[])
{
    long in_it = 0, count_it = 0;
    long whole_part_it = 0, remainder = 0; // Non-uniform grid

    if (argc == 5)
    {
        nfile = atoi(argv[1]);
        in_it = atol(argv[2]);
        count_it = atol(argv[3]);
        del_file_L = argv[4][0];
    }
    else
    {
        printf("\n-----> YOU MUST HAVE THE FOLDER 'init%d\\' \n\n", nfile * 1000);
        printf("It is saved every 10000 iterations\n");
        cout << "Enter initial iteration: ";
        cin >> in_it;
        cout << "Enter count of iterations: N= ";
        cin >> count_it;
        cout << "Rewrite files L" << nfile << " and Time" << nfile << "? (y/n): ";
        cin >> del_file_L;
    }

    cout << nfile << " " << in_it << " " << count_it << " " << del_file_L;
    whole_part_it = int(in_it / savepoints);
    remainder = in_it % savepoints;
    N = remainder + count_it;
    file_read = whole_part_it * savepoints;

    printf("\nCount of iterations=%d \tNumber file for read=%d\n", N, file_read);

    double n0[M + 1], v0[M + 1], T0[M + 1];
    // Initial conditions
    // double L0=4e-13;
    for (i = 0; i <= M; i++)
    {
        ksi[i] = double(i) / double(M);
        r[i] = pow(Rm, ksi[i]);
        n0[i] = 1.0 * exp(beta_0 * (-1. + 1. / r[i]));
        v0[i] = 0.1 * (r[i] - 1.);
        // v0[i] = L0/(pow(r[i], 2)*n0[i]);
        T0[i] = 1.;
        P0[i] = n0[i] * T0[i];
    }
    printf("v0[M]=%g\n C[M]=%g\n", v0[M], sqrt(Gamma * P0[M] / n0[M]));

    /*for (i=0; i<5; i++) {
        //printf("r=%g\n", r[i]);
        printf("r_nonUn=%g\t n0=%g\t v0=%g\n", r[i],n0[i],v0[i]);
    }
    printf("\n");
    for (i=M-5; i<M+1; i++) {
        //printf("r=%g\n", r[i]);
        printf("r_nonUn=%g\t n0=%g\t v0=%g\n", r[i],n0[i],v0[i]);
    }*/
    printf("beta_s = %f\n", beta_s);
    printf("beta = %f\n", beta_0);
    printf("d_sp = %f\n", d_sp);
    printf("a = %f\n", a);
    printf("qN = %f\n", qN);
    printf("R0d = %e\n", R0d);
    printf("N0 = %e\n", density_0);
    printf("a_T = %e\n", a_T);
    printf("b_T = %e\n", b_T);
    printf("ai = %e\n", ai);
    printf("ai1 = %e\n", ai1);

    printf("Arec = %e\n", Arec);
    printf("Arec1 = %e\n", Arec1);
    printf("Arec2 = %e\n", Arec2);

    printf("Mp = %e\n", M_p);

    RK4_for_system(n0, v0, T0);
    return 0;
}
