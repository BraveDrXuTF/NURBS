#include "NRUBS.h"
#include "string.h"
NRUBS::NRUBS(int n)
{
    binCalculate(n);
}
int NRUBS::FindSpan(int n, int p, double u, double* U)
{
    /*A2.1查找u所在的节点区间，返回下标*/
    if (u == U[n + 1])  return(n);
    int low = p;  int high = n + 1;
    int mid = (low + high) / 2;
    while (u < U[mid] || u >= U[mid + 1])
    {
        if (u < U[mid])  high = mid;
        else   low = mid;
        mid = (low + high) / 2;
    }
    return(mid);
}
void NRUBS::BasisFuns(int i, double u, double p, double* U, double* N)
{
    /*A2.2计算所有非0B样条函数的值*/
    N[0] = 1.0;
    double* left = new double[p + 1];
    double* right = new double[p + 1];
    for (int j = 1; j <= p; j++)
    {
        left[j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; r++)
        {
            double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}
void NRUBS::DersBasisFuns(int i,double u,int p,int n,double* U,double** ders)
{ 
    /*A2.3计算非0b样条函数及其导数*/
    double** ndu = new double* [p + 1];
    for (int i = 0; i < p + 1; i++) {
        ndu[i] = new double[p + 1];
    }
    ndu[0][0] = 1.0;
    double* left = new double[p + 1];
    double* right = new double[p + 1];
    for (int j = 1; j <= p; j++)
    {
        left[j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; r++)
        { 
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= p; j++)  
        ders[0][j] = ndu[j][p];
    double** a = new double* [2];
    a[0] = new double[p + 1];
    a[1] = new double[p + 1];
    for (int r = 0; r <= p; r++)  
    {
        int s1, s2, j1, j2;
        s1 = 0; s2 = 1;
        a[0][0] = 1.0;

        for (int k = 1; k <= n; k++)
        {
            double d = 0.0;
            int rk = r - k;  int pk = p - k;
            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                double d = a[s2][0] * ndu[rk][pk];
            }
            if (rk >= -1)    j1 = 1;
            else         j1 = -rk;
            if (r - 1 <= pk) j2 = k - 1;
            else           j2 = p - r;
            for (int j = j1; j <= j2; j++)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][k] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            int j;
            j = s1;  s1 = s2; s2 = j;  
        }
    }
    int r;
    r = p;
    for (int k = 1; k <= n; k++)
    {
        for (int j = 0; j <= p; j++)  ders[k][j] *= r;
        r *= (p - k);
    }
}

void NRUBS::SurfacePoint_B(int n, int p, double* U, int m, int q, double* V, double** P, double u, double v, double S, double* Nu, double* Nv)
{ 
    int uspan = FindSpan(n, p, u, U);
    BasisFuns(uspan, u, p, U, Nu);
    int vspan = FindSpan(m, q, v, V);
    BasisFuns(vspan, v, q, V, Nv);
    int uind = uspan - p;
    double S = 0.0;
    for (int l = 0; l <= q; l++)
    {
        double temp = 0.0;
        int vind = vspan - q - l;
        for (int k = 0; k <= p; k++)
            temp = temp + Nu[k] * P[uind + k][vind];
        S = S + Nv[l] * temp;
    }
}
void NRUBS::SurfaceDerivsAlg1(int n,int p,double* U,int m,int q,double* V,double** P,double u,double v,int d,double** SKL,double** Nu,double** Nv)
{
    /*A3.6计算B样条曲面上的点及其所有直到d阶偏导矢*/
    int du = d > p ? p : d;
    for (int k = p + 1; k <= d; k++)
        for (int l = 0; l <= d - k; l++)   SKL[k][l] = 0.0;
    int dv = d > q ? q : d;
    for (int l = q + 1; l <= d; l++)
        for (int k = 0; k <= d - l; k++)   SKL[k][l] = 0.0;
    int uspan = FindSpan(n, p, u, U);
    DersBasisFuns(uspan, u, p, du, U, Nu);
    int vspan = FindSpan(m, q, v, V);
    DersBasisFuns(vspan, v, q, dv, V, Nv);
    double* temp = new double[q + 1];
    for (int k = 0; k <= du; k++)
    {
        for (int s = 0; s <= q; s++)
        {
            temp[s] = 0.0;
            for (int r = 0; r <= p; r++)
                temp[s] = temp[s] + Nu[k][r] * P[uspan - p + r][vspan - q + s];
        }
        int dd = d - k < dv ? (d - k) : dv;
        for (int l = 0; l <= dd; l++)
        {
            SKL[k][l] = 0.0;
            for (int s = 0; s <= q; s++)
                SKL[k][l] = SKL[k][l] + Nv[l][s] * temp[s];
        }
    }
}
void NRUBS::SurfacePoint(int n, int p, double* U, int m, int q, double* V, double** Pw, double u, double v, double S,double* Nu,double* Nv)
{ 
    /*A4.3计算NRUBS（有理曲面）上的点，当然这里只能计算分量*/
    int uspan = FindSpan(n, p, u, U);
    BasisFuns(uspan, u, p, U, Nu);
    int vspan = FindSpan(m, q, v, V);
    BasisFuns(vspan, v, q, V, Nv);
    S = 0.0;
    double* temp = new double[q + 1];
    for (int l = 0; l <= q; l++)
    {
        temp[l] = 0.0;
        for (int k = 0; k <= p; k++)
            temp[l] = temp[l] + Nu[k] * Pw[uspan - p + k][vspan - q + l];
    }
    double Sw = 0.0;
    for (int l = 0; l <= q; l++)
        Sw = Sw + Nv[l] * temp[l];
    S = Sw / w;
}
void NRUBS::RatSurfaceDerivs(double** Aders, double** wders, int d, double** SKL)
{
    /*A4.4计算S(u,v)的偏导矢*/
    for (int k = 0; k <= d; k++)
        for (int l = 0; l <= d - k; l++)
        {
            double v = Aders[k][l];
            for (int j = 1; j <= l; j++)
                v = v - Bin[l][j] * wders[0][j] * SKL[k][l - j];
            for (int i = 1; i <= k; i++)
            {
                v = v - Bin[k][i] * wders[i][0] * SKL[k - i][l];
                double v2 = 0.0;
                for (int j = 1; j <= l; j++)
                    v2 = v2 + Bin[l][j] * wders[i][j] * SKL[k - i][l - j];
                v = v - Bin[k][i] * v2;
            }
            SKL[k][l] = v / wders[0][0];
        }
}
void NRUBS::binCalculate(int n) {
    /*计算组合数，初始化Bin i:下标，j：上标*/
    Bin = new double* [n];
    for (int i = 0; i < n; i++) {
        Bin[i] = new double[n];
        memset(Bin[i], 0, sizeof(Bin[i]));
        Bin[i][0] = 1;
    }
    for (int i = 1; i < n; i++) {
        for (int j = 1; j <= i; j++) {
            Bin[i][j] = Bin[i - 1][j] + Bin[i - 1][j - 1];
        }
    }

}


