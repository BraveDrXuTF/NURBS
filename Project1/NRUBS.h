#pragma once
class NRUBS
{
    
private:
    double** Bin;
    double w;
    
public:
    int FindSpan(int n, int p, double u, double* U);
    void RatSurfaceDerivs(double** Aders, double** wders, int d, double** SKL);
    void binCalculate(int n);
    void BasisFuns(int i, double u, double p, double* U, double* N);
    void SurfacePoint(int n, int p, double* U, int m, int q, double* V, double** Pw, double u, double v, double S, double* Nu, double* Nv);
    void SurfacePoint_B(int n, int p, double* U, int m, int q, double* V, double** Pw, double u, double v, double S, double* Nu, double* Nv);
    void DersBasisFuns(int i, double u, int p, int n, double* U, double** ders);
    void SurfaceDerivsAlg1(int n, int p, double* U, int m, int q, double* V, double** P, double u, double v, int d, double** SKL, double** Nu, double** Nv);
    NRUBS(int n);
};

