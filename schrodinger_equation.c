#include <stdio.h>
#include <math.h>
#include <complex.h> 

#define N_steps 100000
#define L 1000
#define dt 1
#define dx 4.0
#define w 0.002

double V(int x_rel)
{
	return - pow(w,2) * pow(x_rel - 500,2) / 2.0; //SHO

	//return 0; 
}

double complex b(int i)
{
	return 0.5 * dt * I * (1.0 / pow(dx,2.0) + V(i)) + 1.0;
}

double u_0(int x_rel)
{
	//return sqrt(2 / L) * sin(5*x_rel * M_PI / L); //eigenstate of infinite square well

	//return (1.0/ sqrt(5 * M_PI)) * exp(I * 0.5 * x_rel) * exp(-pow(x_rel - 500, 2) / (2 * 25)); //gaussian packet

	return  (w / M_PI) * exp(-w * pow(x_rel - 500, 2) / 2.0);
}

void CN(double complex *u, double complex *u_aux, double complex *u_next, double complex a, int bi)
{
	//bi = 1 (bounded); bi = 0 (periodic)

	int i;

	for(i = bi; i <= L - bi; i++) u_aux[i] = conj(a) * (u[(i-1+L)%L] + u[(i+1+L)%L]) + conj(b(i)) * u[i];

	//thomas algorithm

	double complex c_new[L+1], d_new[L+1];

	c_new[bi] = a / b(bi);
	for(i = 1 + bi; i <= L - bi; i++) c_new[i] = a / (b(i) - c_new[i-1] * a);

	d_new[bi] = u_aux[bi] / b(bi);
	for(i = 1 + bi; i <= L - bi; i++) d_new[i] = (u_aux[i] - d_new[i-1] * a) / (b(i) - c_new[i-1] * a);

	if (bi == 1) u_next[0] = u_next[L] = 0;
	u_next[L-bi] = d_new[L-bi];
	for(i = L-1-bi; i >= bi; i --) u_next[i] = d_new[i] - c_new[i] * u_next[i+1];
	
	//u = u_next

	for(i = 0; i <= L; i++) u[i] = u_next[i];
} 

int main(void)
{
	int i, j, n = 0;

	double complex u[L+1], u_next[L+1], u_aux[L+1], a = - 0.25 * I * dt / pow(dx,2.0);

	//Initial Contition
	for(i = 0; i <= L; i ++) u[i] = u_0(i);

	while(n < N_steps)
	{
		CN(u, u_aux, u_next, a, 0);
		printf("set title 'Time = %d'\nplot \'-' w lp pt 7 ps 0.1", n);
		for(i = 0; i <= L; i ++) printf("%d\t%.10lf\n",i,pow(creal(u[i]),2) + pow(cimag(u[i]),2));
		//for(i = 0; i <= L; i ++) printf("%d\t%lf\n",i,cimag(u[i]));
		//for(i = 0; i <= L; i ++) printf("%d\t%lf\n",i,creal(u[i]));
		printf("e\npause 0.1\n"); 
		n ++;
	}

	return 0;
	
}	
