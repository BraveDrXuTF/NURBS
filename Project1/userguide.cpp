#include "NRUBS.h"
#include<iostream>
int main() {
	
	/*�㷨A2.1*/
	/*ʹ������3.2��B����������Ϊʾ��*/
	int m = 10 ,p = 3;
	int n = m - p - 1;
	NRUBS mynrubs(10);//���캯����ô�ظ����ã�
	double U[11] = { 0,0,0,0,2 / 5,3 / 5,3 / 5,1,1,1,1 };
	int ids = mynrubs.FindSpan(n, p, 0.75, U);
	printf("0.75�ڵ�%d������", ids);
	/*�㷨A2.2*/
	/*ʹ������2.3(�Լ�2.4)��B����������Ϊʾ��*/
	m = 10, p = 2;
	n = m - p - 1;
	double u = 2.5;
	double U1[11] = { 0,0,0,1,2,3,4,4,5,5,5 };
	double* N = new double[p + 1];
	ids = mynrubs.FindSpan(n, p, u, U1);
	printf("A2.2 %f�ڵ�%d������", u,ids);
	mynrubs.BasisFuns(ids, u, p, U1, N);
	double** der_2_4 = new double*[n + 1];
	for (int o = 0; o < n+1; o++) {
		der_2_4[o] = new double[p + 1];
	}
	mynrubs.DersBasisFuns(ids, u, p, n, U1, der_2_4);

	std::cout <<"p=2ʱ��B������������"<< N[0] <<' '<< ' '<<N[1] <<' '<< N[2]<<'\n';
	for (int k = 0; k < n+1; k++) {
		for (int j = 0; j < p+1; j++) {
			printf("N%d ,%d��%d�׵��� :der_2_4[%d][%d] %f\n", ids - p + j, p, k,k,j,der_2_4[k][j]);
		}
	}
	double U_Bsurface[] = {0,0,0,2/5,3/5,1,1,1};
	double V_Bsurface[] = { 0,0,0,1 / 5,1 / 2,4 / 5,1,1,1 };
	double u = 1 / 5;
	double v = 3 / 5;
	mynrubs.SurfaceDerivsAlg1()

}