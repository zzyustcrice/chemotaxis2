#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>// OpenMP编程需要包含的头文件
#include<string.h>
#include<vector>
using namespace std;
#define DEBUG 1
#define CHEMOTAXIS 1
#define NCPU 8
#define N 300000
#define PI 3.14159
#define grid 500 // grid number
#define length 2 // grid length 2um
typedef struct Particle {
	double x, y, h; // recording particle position and heading direction
	double x0, y0, h0;
	double y1, y2, y3;
	int state;
	int tc;
}Particles;
typedef struct PeriodReve {
	int i, t, t0;
	double  x0, y0, x1, y1;
}PeriodRev;
typedef struct node {
	int val;
	struct node *next;
}node_t;
int Rev = 0;
int main() {
	int  i, j, l, i1, i2, j1, j2;
	double r = 0;
	double periodlength[400] = { 0 };
	int number[400] = { 0 };
	int particlestep(FILE *p8, vector<vector<PeriodRev> > &revdata, int t, int number[400], double periodlength[400], double ff, double f[50], int i, double lumda, double v0, double timelength, node_t Grid[grid][grid], Particles temp[N], Particles par[N], double d, double ea, double b, int n[grid][grid], double M[2 * length*grid][2 * length*grid], double R[2 * length*grid][2 * length*grid], double A[length*grid / 5][length*grid / 5], double R1[length*grid / 5][length*grid / 5], double coe[length*grid / 5][length*grid / 5], double tavg);
	char s[50];
	static double M[2 * length*grid][2 * length*grid] = {}, R[2 * length*grid][2 * length*grid] = {};//recording memory
	static double A[length*grid / 5][length*grid / 5] = {}, R1[length*grid / 5][length*grid / 5] = {}, coe[length*grid / 5][length*grid / 5] = {};//recording chemotax
	double timelength = 0.05, v0 = 5;
	double f[50] = { 0 }, ff = 0;
	double noisetime = 1;
	int k = 40;
	int tempint;
	double tempdouble;
	double q = 2;
	int t = 0, tstep = 18000; //time step
	double totnoise = 0;
	static int n[grid][grid] = { 0 };//record visit frequency
	double d = 75, a = 2 * timelength, ea, b = 0.25*timelength;//decay exponent
	double lumda = d * timelength / 25;//lumda=Ddt/dl^2
	double *coep[length*grid / 5] = {};
	double tavg = 3.3;
	static Particles par[N] = {};
	static Particles temp[N] = {};
	static node_t Grid[grid][grid] = {};
	int matrixInversion(double **a, int n);
	double myrandom();
	void push(node_t * head, int val);
	double initializeAngle();
	vector<vector<double> > rundurations(48, vector<double>(2, 0));
	double position = 0;
	int sizeofchefield = length * grid / 5;
	//7/20/2016{
	vector<vector<PeriodRev> > revdata(NCPU);

	//}
	node_t *ListDelete(node_t *currP, int value);
	//7/20/2016{
	FILE *p1, *p2, *p3, *p4, *p5, *p6, *p0, *p7, *p8, *p9;
	p0 = fopen("disaway.txt", "w");
	p1 = fopen("distoward.txt", "w");
	p6 = fopen("numaway.txt", "w");
	p3 = fopen("numtoward.txt", "w");
	p7 = fopen("signal.txt", "w");
	p2 = fopen("memory-x.txt", "w");
	p4 = fopen("periodlength.txt", "w");
	p8 = fopen("RevData.txt", "w");
	p9 = fopen("runduration.txt", "rb");
	//}
	for (i = 0; i < 48; i++) {
		fscanf(p9, "%lf	%lf", &rundurations[i][0], &rundurations[i][1]);
	}
	// initialize particles coordinates
	srand(time(NULL));
	ea = exp(-a);
	memset(M, 0, sizeof(M));
	//initialize linked list

	for (i = 0; i<grid; i++)
		for (j = 0; j<grid; j++)
		{
			Grid[i][j].next = NULL;
			Grid[i][j].val = -1;
		}

	for (i = 0; i<50; i++)
	{
		f[i] = pow(i, q) / (pow(i, q) + pow(k, q));
	}
	memset(A, 0, sizeof(A));
	//initialize particles' status
#pragma omp parallel for num_threads(NCPU) private(r)
	for (i = 0; i<N; i++)
	{
		//		  fscanf(p0,"%lf, %lf, %lf\r\n",&par[i].x,&par[i].y,&par[i].h);
		r = myrandom(); if (r == 1)r = 0;
		par[i].x = grid * length*r;
		r = myrandom(); if (r == 1)r = 0;
		par[i].y = grid * length*r;
		par[i].h = initializeAngle();
		r = myrandom();
		par[i].tc = tavg / timelength * r;
		par[i].state = 0;
		par[i].y1 = par[i].y2 = par[i].y3 = 0;

		push(&Grid[(int)(par[i].x / length)][(int)(par[i].y / length)], i);
	}
	for (i = 0; i<sizeofchefield; i++)
	{
		coep[i] = coe[i];
		i1 = i - 1; i2 = i + 1;
		if (i1<0)i1 += sizeofchefield;
		if (i2 >= sizeofchefield)i2 -= sizeofchefield;
		coe[i][i] = 1 + 2 * lumda;
		coe[i][i1] = -lumda;
		coe[i][i2] = -lumda;
	}
	matrixInversion(coep, sizeofchefield);
	for (t = 1; t < tstep; t++)
	{
		//		printf("%d\n", t);
		//		while (rundurations[position][0] < (t - 6000)*timelength / 60.0 && position<47)position++;
		//		tavg = rundurations[position][1];
		memset(R1, 0, sizeof(R1));
#pragma omp parallel for num_threads(NCPU) private(j)
		for (i = 0; i < 2 * length*grid; i++)
		{
			for (j = 0; j < 2 * length*grid; j++)
			{
				M[i][j] = M[i][j] * ea;// memory decay
			}
		}

		//noise in direction
		totnoise += timelength / noisetime;
		if (totnoise>1) {
			totnoise--;
#pragma omp parallel for num_threads(NCPU) private(r)
			for (i = 0; i<N; i++) {
				r = 2 * myrandom();
				r = (r - 1)*sqrt(PI*PI / 225 * noisetime);
				par[i].h += r;
				while (par[i].h >= PI) { par[i].h -= 2 * PI; }
				while (par[i].h<-PI) { par[i].h += 2 * PI; }

			}
		}
#pragma omp parallel for num_threads(NCPU)
		for (i = 0; i < N; i++)
		{
			temp[i] = par[i];// record the particles position
			n[(int)(par[i].x / length)][(int)(par[i].y / length)]++;
		}
#pragma omp parallel for num_threads(NCPU) private(r)
		for (i = 0; i < N; i++)
		{
			r = myrandom()*0.4 - 0.2 + 1;
			particlestep(p8, revdata, t, number, periodlength, ff, f, i, lumda, v0*r, timelength, Grid, temp, par, d, ea, b, n, M, R, A, R1, coe, tavg);

		}
		if (CHEMOTAXIS) {

			for (i = 0; i < sizeofchefield; i++)
				for (j = 0; j < sizeofchefield; j++)
				{
					A[i][j] = A[i][j] + R1[i][j]; R1[i][j] = A[i][j];
				}
			//diffuse
			for (int diffuloop = 0; diffuloop < 4; diffuloop++) {
				if (diffuloop % 2 == 0) {

					for (i = 0; i < sizeofchefield; i++) {
						i1 = i - 1; i2 = i + 1;
						if (i1 < 0)i1 += sizeofchefield;
						if (i2 >= sizeofchefield)i2 = 0;
						for (j = 0; j < sizeofchefield; j++)
							R1[i][j] = (1 - 2 * lumda - b)*A[i][j] + lumda * A[i1][j] + lumda * A[i2][j];

					}

					for (j = 0; j < sizeofchefield; j++)
						for (i = 0; i < sizeofchefield; i++) {
							A[i][j] = 0;
							for (l = 0; l < sizeofchefield; l++)
								A[i][j] += R1[i][l] * coe[l][j];
						}
				}
				else {
					for (j = 0; j < sizeofchefield; j++)
					{
						j1 = j - 1; j2 = j + 1;
						if (j1 < 0)j1 += sizeofchefield;
						if (j2 >= sizeofchefield)j2 = 0;
						for (i = 0; i < sizeofchefield; i++)
							R1[i][j] = (1 - 2 * lumda - b)*A[i][j] + lumda * A[i][j1] + lumda * A[i][j2];
					}
					for (j = 0; j < sizeofchefield; j++)
						for (i = 0; i < sizeofchefield; i++) {
							A[i][j] = 0;
							for (l = 0; l < sizeofchefield; l++)
								A[i][j] += R1[l][j] * coe[i][l];
						}
				}
			}
		}
		//update linked list
		for (i = 0; i<N; i++) {
			ListDelete(&Grid[(int)(temp[i].x / length)][(int)(temp[i].y / length)], i);
			push(&Grid[(int)(par[i].x / length)][(int)(par[i].y / length)], i);
		}


		//          update memory
#pragma omp parallel for num_threads(NCPU) private(j)
		for (i = 0; i<2 * length*grid; i++)
			for (j = 0; j<2 * length*grid; j++)
			{
				M[i][j] = M[i][j] + R[i][j]; R[i][j] = 0;
			}

		//out put reversal data

		for (tempint = 0; tempint < NCPU; tempint++) {
			for (i = 0; i<revdata[tempint].size(); i++)
				fprintf(p8, "%d, %d, %lf, %lf, %lf, %lf\n", revdata[tempint][i].t0, revdata[tempint][i].t, revdata[tempint][i].x0, revdata[tempint][i].y0, revdata[tempint][i].x1, revdata[tempint][i].y1);
			revdata[tempint].clear();
		}




		if (t % 120 == 0)
		{
			sprintf(s, "particles%d.txt", t / 120);
			p5 = fopen(s, "w");
			//err handle
			if (p5 == NULL)printf("t=%d, fopen err", t);
			else {
				for (i = 0; i < N; i++)
				{
					fprintf(p5, "%lf, %lf, %lf, %lf, %d, %d\r\n", par[i].x, par[i].y, par[i].y2, par[i].h, par[i].tc, par[i].state);

				}
			}
			fclose(p5);

		}
	}//printf("%d\n",rotation);

	for (i = 0; i<grid; i++)
	{
		for (j = 0; j<grid; j++)
		{
			fprintf(p4, "%d ", n[j][i]);

		}
		fprintf(p4, "\r\n");

	}
	for (i = 0; i<2 * length*grid; i++)
	{
		for (j = 0; j<2 * length*grid; j++) {
			fprintf(p2, "%lf ", M[j][i]);
		}fprintf(p2, "\r\n");
	}
	fclose(p2);
	fclose(p4);

	for (i = 0; i< sizeofchefield; i++)
	{
		for (j = 0; j< sizeofchefield; j++) {
			fprintf(p7, "%lf ", A[j][i]);
		}fprintf(p7, "\r\n");
	}
	fclose(p6);
	fclose(p0);
	fclose(p1);
	fclose(p3);
	fclose(p7);
	fclose(p8);



}
int particlestep(FILE *p8, vector<vector<PeriodRev> > &revdata, int t, int number[400], double periodlength[400], double ff, double f[50], int i, double lumda, double v0, double timelength, node_t Grid[grid][grid], Particles temp[N], Particles par[N], double d, double ea, double b, int n[grid][grid], double M[2 * length*grid][2 * length*grid], double R[2 * length*grid][2 * length*grid], double A[length*grid / 5][length*grid / 5], double R1[length*grid / 5][length*grid / 5], double coe[length*grid / 5][length*grid / 5], double tavg) {

	int  j = 0, l = 0, i2 = 0, j2 = 0, tempint = 0, sum = 0;
	int n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0, n6 = 0, n7 = 0;
	double tempdouble = 0;
	double xxx = 0;
	double ps = 0;//possibility of staying
	double x1 = 0, y1 = 0;
	double v = 0;
	double r = 0;
	double myrandom(), nearbycount(int i, double x, double y, double radius, Particles par[N], node_t Grid[grid][grid]);
	double angle(int i, double hd[5], double beta, double dt, node_t Grid[grid][grid], Particles temp[N], Particles par[N]);
	double slimedeposit(double dt, double x1, double y1, double x, double y, double R[2 * length*grid][2 * length*grid]);
	double ReversingProbability(double m, double kapa, int tc, double dt);
	double beta = 0, c = 0.0;
	double hd[5] = { 0 };//record memory in each sector ahead
	double dis = 0;
	double theta = 0;
	double slimeradius = 3;
	double kapa = 0, m = 0, tdvia = 0.5;
	double tadpt = 1;
	double pchem = 1.5;
	double chemsignal = 0, tempdouble2 = 0;
	double initializeAngle();
	PeriodRev singlerev;
	sum = nearbycount(i, par[i].x, par[i].y, 5, temp, Grid);


	par[i].tc++;
	tempdouble = A[(int)(par[i].x) / 5][(int)(par[i].y) / 5];
	if (tempdouble<0)tempdouble = 0;
	chemsignal = pow(tempdouble, 1) / (pow(45, 1) + pow(tempdouble, 1));
	tempdouble2 = (chemsignal - par[i].y2);
	par[i].y2 += (tempdouble2) / tadpt * timelength;
	if (par[i].state == 0)
	{
		tempdouble = tavg * exp(11.3* chemsignal - 12.4*par[i].y2);
		m = tempdouble * tempdouble / tdvia / tdvia;
		kapa = tempdouble / tdvia / tdvia;
		r = myrandom();
		ff = ReversingProbability(m, kapa, par[i].tc, timelength);
		if (r < ff) {
			if (t > 6000) {
				singlerev.i = i;
				singlerev.t0 = t - par[i].tc;
				singlerev.t = t;
				singlerev.x0 = par[i].x0;
				singlerev.y0 = par[i].y0;
				singlerev.x1 = par[i].x;
				singlerev.y1 = par[i].y;
				revdata[i / (N / NCPU)].push_back(singlerev);
			}



			//		number[t / 30]++;
			//		periodlength[t / 30] += par[i].tc*timelength;
			par[i].tc = 0;
			par[i].x0 = par[i].x;
			par[i].y0 = par[i].y;
			tempdouble = A[(int)(par[i].x) / 5][(int)(par[i].y) / 5];
			if (tempdouble < 0)tempdouble = 0;
			tempdouble = 0.25 *(pow(15, 2) + 2.5*pow(tempdouble, 2)) / (pow(15, 2) + pow(tempdouble, 2));
			r = myrandom();
			if (r < tempdouble)par[i].state = 1;
			else {
				if (par[i].h>0)par[i].h -= PI;//reversal
				else par[i].h += PI;
			}

		}
		if (sum >= 50)ps = 1;
		else ps = f[sum];
		v = v0 * (1 - 0.25* ps);
		dis = v * timelength;
	}
	if (par[i].state == 1) {
		tempdouble = A[(int)(par[i].x) / 5][(int)(par[i].y) / 5];
		if (tempdouble < 0)tempdouble = 0;
		tempdouble = 2.25 *(pow(24, 2) + 3.7*pow(tempdouble, 2)) / (pow(24, 2) + pow(tempdouble, 2));
		if (par[i].tc>tempdouble / timelength)
		{
			par[i].x0 = par[i].x;
			par[i].y0 = par[i].y;
			par[i].tc = 0;
			par[i].state = 0;
			par[i].h = initializeAngle();
		}
		if (sum >= 50)ps = 1;
		else ps = f[sum];
		v = v0 * (1 - 0.25* ps);
		dis = v / 4 * timelength;
	}



	//calculate slime value in each sector
	for (j = floor(2 * par[i].x - 2 * slimeradius - 1); j <= floor(2 * par[i].x + 2 * slimeradius + 1); j++)
	{
		i2 = j;
		if (j < 0)i2 = j + 2 * length*grid;
		else if (j>2 * length*grid - 1)i2 = j - 2 * length*grid;
		for (l = floor(2 * par[i].y - 2 * slimeradius - 1); l <= floor(2 * par[i].y + 2 * slimeradius + 1); l++)
		{
			j2 = l;
			if (l < 0)j2 = l + 2 * length*grid;
			else if (l>2 * length*grid - 1)j2 = l - 2 * length*grid;
			if (pow(j + 0.5 - 2 * par[i].x, 2) + pow(l + 0.5 - 2 * par[i].y, 2) < 4 * slimeradius*slimeradius&&M[i2][j2]>exp(-3)) {
				theta = atan(((double)(l)*0.5 + 0.25 - par[i].y) / ((double)(j)*0.5 + 0.25 - par[i].x));
				if ((double)(j)*0.5 + 0.25 - par[i].x < 0)theta += PI;
				theta -= par[i].h;
				while (theta >  PI)theta -= 2 * PI;
				while (theta <= -PI)theta += 2 * PI;

				if (theta > -PI / 2 && theta <= -0.3*PI)hd[0] += M[i2][j2];
				else if (theta > -0.3*PI && theta <= -0.1*PI)hd[1] += M[i2][j2];
				else if (theta > -0.1*PI && theta <= 0.1*PI)hd[2] += M[i2][j2];
				else if (theta > 0.1*PI &&theta <= 0.3*PI)hd[3] += M[i2][j2];
				else if (theta > 0.3*PI && theta <= 0.5*PI)hd[4] += M[i2][j2];
			}
		}
	}


	r = angle(i, hd, beta, timelength, Grid, temp, par);
	r = r + par[i].h;
	x1 = par[i].x + dis * cos(r);
	y1 = par[i].y + dis * sin(r);
	//record deposited memory
	slimedeposit(timelength, x1, y1, par[i].x, par[i].y, R);
	//update particle status

	par[i].h = r;
	while (par[i].h >= PI) { par[i].h -= 2 * PI; }
	while (par[i].h<-PI) { par[i].h += 2 * PI; }
	//apply periodic boundary condition
	if (x1<0)x1 += grid * length;
	else if (x1 >= grid * length)x1 -= grid * length;
	if (y1<0)y1 += grid * length;
	else if (y1 >= grid * length)y1 -= grid * length;
	if (CHEMOTAXIS&&t > 6000)
		R1[(int)(x1) / 5][(int)(y1) / 5] += pchem * timelength*(0.8*pow(A[(int)(x1) / 5][(int)(y1) / 5], 3) / (pow(A[(int)(x1) / 5][(int)(y1) / 5], 3) + pow(10, 3)) + 0.2);
	par[i].x = x1;
	par[i].y = y1;
	n[(int)(par[i].x / length)][(int)(par[i].y / length)]++;
	return(0);
}
double ReversingProbability(double m, double kapa, int tc, double timelength) {
	double tempdouble = 0;
	double *ff;
	int i;
	ff = new double[tc];
	for (i = 0; i < tc; i++) {
		ff[i] = pow(kapa, m)*pow(i*timelength, m - 1)*exp(-kapa * i*timelength) / tgamma(m)*timelength;
		tempdouble += ff[i];

	}
	tempdouble = ff[tc - 1] / (1 - tempdouble + ff[tc - 1]);
	free(ff);
	return(tempdouble);

}

double slimedeposit(double dt, double x1, double y1, double x, double y, double R[2 * length*grid][2 * length*grid]) {
	double y2, x2;
	int j, j1, j2, j3;
	int radius = 0;
	double psli = 20;
	if (fabs(x1 - x)>fabs(y1 - y)) {
		if (fabs(x1 - x) < 0.5) {
			for (j = (int)(2 * x) - 2 * radius; j <= (int)(2 * x) + 2 * radius; j++)
				for (j1 = (int)(2 * y) - 2 * radius; j1 <= (int)(2 * y) + 2 * radius; j1++)
				{
					j2 = j, j3 = j1;
					if (j2 < 0)j2 += 2 * length*grid;
					else if (j2 >= 2 * length*grid)j2 -= 2 * length*grid;
					if (j3 < 0)j3 += 2 * length*grid;
					else if (j3 >= 2 * length*grid)j3 -= 2 * length*grid;
					R[j2][j3] += psli * dt;
				}return 0;
		}
		else {
			if (x1 > x) {
				for (j1 = floor(2 * x); j1 <= floor(2 * x1); j1++) {
					y2 = y + (j1*0.5 + 0.25 - x)*(y1 - y) / (x1 - x);
					j = j1;
					for (j2 = (int)(2 * y2) - 2 * radius; j2 <= (int)(2 * y2) + 2 * radius; j2++) {
						j3 = j2;
						//now apply periodic boundary
						if (j < 0)j += 2 * length*grid;
						else if (j >= 2 * length*grid)j -= 2 * length*grid;
						if (j3<0)j3 += 2 * grid*length;
						else if (j3>grid*length)j3 -= 2 * grid*length;
						R[j][j3] += psli * dt;
					}
				}
			}
			else {
				for (j1 = floor(2 * x); j1 >= floor(2 * x1); j1--) {
					y2 = y + (j1*0.5 + 0.25 - x)*(y1 - y) / (x1 - x);
					j = j1;
					for (j2 = (int)(2 * y2) - 2 * radius; j2 <= (int)(2 * y2) + 2 * radius; j2++) {
						j3 = j2;
						//now apply periodic boundary
						if (j < 0)j += 2 * length*grid;
						else if (j >= 2 * length*grid)j -= 2 * length*grid;
						if (j3<0)j3 += 2 * grid*length;
						else if (j3>grid*length)j3 -= 2 * grid*length;
						R[j][j3] += psli * dt;
					}
				}
			}
		}
	}
	else {
		if (fabs(y1 - y)<0.5) {
			for (j = (int)(2 * x) - 2 * radius; j <= (int)(2 * x) + 2 * radius; j++)
				for (j1 = (int)(2 * y) - 2 * radius; j1 <= (int)(2 * y) + 2 * radius; j1++)
				{
					j2 = j, j3 = j1;
					if (j2<0)j2 += 2 * length*grid;
					else if (j2 >= 2 * length*grid)j2 -= 2 * length*grid;
					if (j3<0)j3 += 2 * length*grid;
					else if (j3 >= 2 * length*grid)j3 -= 2 * length*grid;
					R[j2][j3] += psli * dt;
				}return 0;
		}
		else {
			if (y1 > y) {
				for (j1 = floor(2 * y); j1 <= floor(2 * y1); j1++) {
					x2 = x + (j1*0.5 + 0.25 - y)*(x1 - x) / (y1 - y);
					j = j1;
					for (j2 = (int)(2 * x2) - 2 * radius; j2 <= (int)(2 * x2) + 2 * radius; j2++) {
						j3 = j2;
						//now apply periodic boundary
						if (j < 0)j += 2 * length*grid;
						else if (j >= 2 * length*grid)j -= 2 * length*grid;
						if (j3<0)j3 += 2 * grid*length;
						else if (j3>grid*length)j3 -= 2 * grid*length;
						R[j][j3] += psli * dt;
					}
				}
			}
			else {
				for (j1 = floor(2 * y); j1 >= floor(2 * y1); j1--) {
					x2 = x + (j1*0.5 + 0.25 - y)*(x1 - x) / (y1 - y);
					j = j1;
					for (j2 = (int)(2 * x2) - 2 * radius; j2 <= (int)(2 * x2) + 2 * radius; j2++) {
						j3 = j2;
						//now apply periodic boundary
						if (j < 0)j += 2 * length*grid;
						else if (j >= 2 * length*grid)j -= 2 * length*grid;
						if (j3<0)j3 += 2 * grid*length;
						else if (j3>grid*length)j3 -= 2 * grid*length;
						R[j][j3] += psli * dt;
					}
				}
			}
		}
	}
}
double myrandom() {
	return(rand() / (double)RAND_MAX);

}
int matrixInversion(double **a, int n)
{
	int *is = new int[n];
	int *js = new int[n];
	int i, j, k;
	double d, p;
	for (k = 0; k < n; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				p = fabs(a[i][j]);
				if (p>d) { d = p; is[k] = i; js[k] = j; }
			}
		if (0.0 == d)
		{
			delete[]is; delete[]js; printf("err**not inv\n");
			return(0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				p = a[k][j];
				a[k][j] = a[is[k]][j];
				a[is[k]][j] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				p = a[i][k];
				a[i][k] = a[i][js[k]];
				a[i][js[k]] = p;
			}
		a[k][k] = 1.0 / a[k][k];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				a[k][j] *= a[k][k];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						a[i][j] -= a[i][k] * a[k][j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				a[i][k] = -a[i][k] * a[k][k];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				p = a[k][j];
				a[k][j] = a[js[k]][j];
				a[js[k]][j] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				p = a[i][k];
				a[i][k] = a[i][is[k]];
				a[i][is[k]] = p;
			}
	}
	delete[] is; delete[] js;
	return(1);
}
double nearbycount(int i, double x, double y, double radius, Particles par[N], node_t Grid[grid][grid]) {
	int i0, i1, j1, j, k, tag = 0, sum = 0;
	double theta = 0;
	node_t * current;
	k = radius / length + 1;
	for (i0 = (int)(floor(x / length - k)); i0 <= x / length + k; i0++)
		for (j = (int)(floor(y / length - k)); j <= y / length + k; j++) {
			x = par[i].x; y = par[i].y;
			//periodic boundary condition
			i1 = i0; j1 = j;
			if (i0<0) { i1 = i0 + grid; x += grid * length; }
			else if (i0>grid - 1) { i1 = i0 - grid; x -= grid * length; }
			if (j<0) { j1 = j + grid; y += grid * length; }
			else if (j>grid - 1) { j1 = j - grid; y -= grid * length; }
			current = &Grid[i1][j1];
			while (current->next != NULL) {
				current = current->next;
				tag = current->val;
				if (pow((x - par[tag].x), 2) + pow((y - par[tag].y), 2) <= radius * radius&&i != tag) {

					//					theta = atan((y - par[tag].y) / (x - par[tag].x));
					//					if ((x - par[tag].x) > 0)theta += PI;
					//					theta -= par[i].h;
					//					while (theta > 1.5*PI)theta -= 2 * PI;
					//					while (theta <= -0.5*PI)theta += 2 * PI;
					//					if (theta <= 0.5*PI)
					//					{
					sum++; //alignment += fabs(sin(par[tag].h - par[i].h));
						   //					}
				}
			}
		}

	return(sum);
}


double angle(int i, double hd[5], double beta, double dt, node_t Grid[grid][grid], Particles temp[N], Particles par[N]) {
	node_t * current;
	double F(double x); double Fad(double x);
	double eali = 0.7, esli = 1.0, eden = 18, r = 0;
	double x, y, x1, y1, f, radius = 5;
	double thetalign = 0, thetag = 0, thetaslime = 0, dis = 0;
	double hdmax = 0;
	double dtr = 0.75, kad = 0, fad;
	int i0, i1, j1, j, k, tag, sum = 0;
	int count = 0;
	//slime direction
	for (j = 0; j<5; j++)
		if (hdmax<hd[j])hdmax = hd[j];
	hdmax = hdmax * 0.8;
	if (hd[2] >= hdmax) { thetaslime = 0; }
	else if (hd[1]>hd[3]) {
		if (hd[1]>hdmax) { thetaslime = -0.2*PI; }
		else if (hd[0]>hd[4]) {
			if (hd[0]>hdmax) { thetaslime = -0.4*PI; }
		}
		else if (hd[4]>hdmax) { thetaslime = 0.4*PI; }
	}
	else if (hd[3]>hdmax) { thetaslime = 0.2*PI; }
	else if (hd[0]>hd[4]) {
		if (hd[0]>hdmax) { thetaslime = -0.4*PI; }
	}
	else if (hd[4]>hdmax) { thetaslime = 0.4*PI; }

	//alignment
	x = temp[i].x; y = temp[i].y;

	k = radius / length + 1;
	for (i0 = (int)(floor(x / length - k)); i0 <= x / length + k; i0++)
		for (j = (int)(floor(y / length - k)); j <= y / length + k; j++) {
			x = temp[i].x; y = temp[i].y;
			//periodic boundary condition
			i1 = i0; j1 = j;
			if (i0<0) { i1 = i0 + grid; x += grid * length; }
			else if (i0>grid - 1) { i1 = i0 - grid; x -= grid * length; }
			if (j<0) { j1 = j + grid; y += grid * length; }
			else if (j>grid - 1) { j1 = j - grid; y -= grid * length; }
			current = &Grid[i1][j1];
			while (current->next != NULL) {
				current = current->next;
				tag = current->val;
				if (pow((x - temp[tag].x), 2) + pow((y - temp[tag].y), 2) < radius*radius) {
					count++;
					thetalign += sin(2 * (temp[i].h - temp[tag].h));


				}
			}
		}


	//density gradient
	x = temp[i].x; y = temp[i].y;
	k = 2; f = 0, fad = 0;
	for (i0 = (int)(floor(x / length - k)); i0 <= x / length + k; i0++)
		for (j = (int)(floor(y / length - k)); j <= y / length + k; j++) {
			x = temp[i].x; y = temp[i].y;
			//periodic boundary condition
			i1 = i0; j1 = j;
			if (i0<0) { i1 = i0 + grid; x += grid * length; }
			else if (i0>grid - 1) { i1 = i0 - grid; x -= grid * length; }
			if (j<0) { j1 = j + grid; y += grid * length; }
			else if (j>grid - 1) { j1 = j - grid; y -= grid * length; }
			current = &Grid[i1][j1];
			while (current->next != NULL) {
				current = current->next;
				tag = current->val;
				if (pow((x - temp[tag].x), 2) + pow((y - temp[tag].y), 2) < 9 && tag != i) {
					x1 = temp[tag].x - x; y1 = temp[tag].y - y;
					dis = x1 * cos(temp[i].h + PI / 2) + y1 * sin(temp[i].h + PI / 2);
					if (fabs(dis) <= 0.5) {
						f += F(dis);
					}
					else if (fabs(dis) <= dtr) {
						fad += Fad(dis) / (dtr - 1);
					}
				}
			}
		}

	par[i].x += eden * f*dt*cos(par[i].h + PI / 2) + kad * fad*dt*cos(par[i].h + PI / 2);
	par[i].y += eden * f*dt*sin(par[i].h + PI / 2) + kad * fad*dt*sin(par[i].h + PI / 2);
	if (count == 0)thetalign = 0;
	else thetalign = thetalign * eali / count;
	r = dt * (0 - thetalign + esli * sin(2 * thetaslime));
	return(r);
}
double F(double x) {
	if (x >= 0)return(x - 0.5);
	else return(-F(-x));
}
double Fad(double x) {
	if (x >= 0)return(x - 0.5);
	else return(-F(-x));
}

double normal(double t) {
	double myrandom();
	double r, r1, r2;
	double sigma = (PI*PI) / 225 / 3 * t*0.88;
	double sqrtsigma = sqrt(sigma);
	r1 = myrandom()*((0.003 + 6 / sqrt(2 * PI)));
	if (r1 < 0.0015) {
		r = -0.5*PI + r1 / (0.003 / (PI - 6 * sqrtsigma));

	}
	else if (r1 < 0.003) {
		r1 -= 0.0015;
		r = 3 * sqrtsigma + r1 / (0.003 / (PI - 6 * sqrtsigma));

	}
	else {

		r1 -= 0.003;
		r = -3 * sqrtsigma + r1 * sqrt(2 * PI)*sqrtsigma;
		r2 = myrandom();
		if (r2>exp(-r * r / (2 * sigma))) {
			r = normal(t);
		}
	}
	return(r);
}
void push(node_t * head, int val) {
	node_t * current = head;
	while (current->next != NULL) {
		current = current->next;
	}

	/* now we can add a new variable */
	current->next = (node_t *)malloc(sizeof(node_t));
	current->next->val = val;
	current->next->next = NULL;
}
node_t *ListDelete(node_t *currP, int value)
{
	/* See if we are at end of list. */
	if (currP == NULL) {
		printf("error\n");
		return NULL;
	}
	/*
	* Check to see if current node is one
	* to be deleted.
	*/
	if (currP->val == value) {
		node_t *tempNextP;

		/* Save the next pointer in the node. */
		tempNextP = currP->next;

		/* Deallocate the node. */
		free(currP);

		/*
		* Return the NEW pointer to where we
		* were called from.  I.e., the pointer
		* the previous call will use to "skip
		* over" the removed node.
		*/
		return tempNextP;
	}

	/*
	* Check the rest of the list, fixing the next
	* pointer in case the next node is the one
	* removed.
	*/
	currP->next = ListDelete(currP->next, value);


	/*
	* Return the pointer to where we were called
	* from.  Since we did not remove this node it
	* will be the same.
	*/
	return currP;
}
double initializeAngle() {
	double myrandom();
	double r = myrandom();
	return(PI*(2 * r - 1));
}
