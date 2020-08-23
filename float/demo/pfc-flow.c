// Author: Petri Hirvonen, petenez@gmail.com, August 2020

// See my GitHub repository github.com/petenez/pfc and the paper Real-Time Fluid Dynamics for Games by Jos Stam https://pdfs.semanticscholar.org/847f/819a4ea14bd789aca8bc88e85e906cfc657c.pdf for further details of the PFC and CFD solvers.

// compile: gcc pfc-flow.c -fopenmp -lfftw3_omp -lfftw3 -lm -O3 -Wall -o pfc-flow
// run: ./pfc-flow

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define pi 3.14159265358979323846264338327

struct S {

	int t;			// iteration
	int tend;		// last iteration
	
	int Nx;			// system width
	int Ny;			// system height
	
	int Tprint;		// print interval
	int Twrite;		// write interval
	int pprint;		// print precision
	int pwrite;		// write precision
	
	double dx;		// horizontal discretization (PFC*)
	double dy;		// vertical discretization (*)
	double dt;		// time step (*)
	
	double alpha;	// quadratic term coefficient (*)
	double sigma;	// smoothing spread (*)
	
	double* A;		// operator for linear part (*)
	double* B;		// operator for nonlinear part (*)
	double* n;		// density (*)
	double* m;		// smoothed density (*)
	
	double mumin;	// minimum viscosity
	double mumax;	// maximum viscosity
	
	double* vx;		// horizontal velocity
	double* vy;		// vertical velocity
	double* mu;		// viscosity
	double* p;		// arrays for intermediate results
	double* q;
	double* r;
	
	fftw_plan n_N;	// forward FFT from n to n
	fftw_plan m_M;	// forward FFT from m to m
	fftw_plan N_n;	// backward FFT from n to n
	fftw_plan M_m;	// backward FFT from m to m

};

// write PFC density
void write_n(struct S* s) {

	// physical array width (padding due to real-data FFT) (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;

	char str[BUFSIZ];
	sprintf(str, "test-%d.n", s->t);
	FILE* file = fopen(str, "w");
	
	fprintf(file, "%d %d %.*lf %.*lf\n", s->Nx, s->Ny, s->pwrite, s->dx, s->pwrite, s->dy);
	
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			fprintf(file, "%.*e\n", s->pwrite, s->n[ij]);
		}
	}
	
	fclose(file);
	
}

// write smoothed PFC density
void write_m(struct S* s) {

	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;

	char str[BUFSIZ];
	sprintf(str, "test-%d.m", s->t);
	FILE* file = fopen(str, "w");
	
	fprintf(file, "%d %d %.*lf %.*lf\n", s->Nx, s->Ny, s->pwrite, s->dx, s->pwrite, s->dy);
	
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			fprintf(file, "%.*e\n", s->pwrite, s->m[ij]);
		}
	}
	
	fclose(file);
	
}

// write flow velocity
void write_v(struct S* s) {

	int Nxp = s->Nx + 2;

	char str[BUFSIZ];
	sprintf(str, "test-%d.v", s->t);
	FILE* file = fopen(str, "w");
	
	fprintf(file, "%d %d %.*lf %.*lf\n", s->Nx, s->Ny, s->pwrite, s->dx, s->pwrite, s->dy);
	
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			fprintf(file, "%.*e %.*e\n", s->pwrite, s->vx[ij], s->pwrite, s->vy[ij]);
		}
	}
	
	fclose(file);
	
}

// allocate arrays
void allocate_data(struct S* s) {

	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;
	int NxhNy = Nxh*s->Ny;
	int NxpNy = Nxp*s->Ny;

	s->A = fftw_malloc(NxhNy*sizeof(double));
	s->B = fftw_malloc(NxhNy*sizeof(double));
	s->n = fftw_malloc(NxpNy*sizeof(double));
	s->m = fftw_malloc(NxpNy*sizeof(double));
	
	Nxp = s->Nx + 2;
	int Nyp = s->Ny + 2;
	int NxpNyp = Nxp*Nyp;
	
	s->vx = fftw_malloc(NxpNyp*sizeof(double));
	s->vy = fftw_malloc(NxpNyp*sizeof(double));
	s->mu = fftw_malloc(NxpNyp*sizeof(double));
	s->p = fftw_malloc(NxpNyp*sizeof(double));
	s->q = fftw_malloc(NxpNyp*sizeof(double));
	s->r = fftw_malloc(NxpNyp*sizeof(double));

}

// set up FFT plans
void plan_FFTs(struct S* s) {
	
	s->n_N = fftw_plan_dft_r2c_2d(s->Ny, s->Nx, s->n, (fftw_complex*)s->n, FFTW_MEASURE);
	s->m_M = fftw_plan_dft_r2c_2d(s->Ny, s->Nx, s->m, (fftw_complex*)s->m, FFTW_MEASURE);
	
	s->N_n = fftw_plan_dft_c2r_2d(s->Ny, s->Nx, (fftw_complex*)s->n, s->n, FFTW_MEASURE);
	s->M_m = fftw_plan_dft_c2r_2d(s->Ny, s->Nx, (fftw_complex*)s->m, s->m, FFTW_MEASURE);
	
}

// get average of PFC array
double ave_PFC(struct S* s, double* a) {

	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;
	
	int T = omp_get_max_threads();
	double a_[T];
	
	for(int t = 0; t < T; ++t) {
		a_[t] = 0.0;
	}
	
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			a_[omp_get_thread_num()] += a[ij];
		}
	}
	
	double sum = 0.0;
	
	for(int t = 0; t < T; ++t) {
		sum += a_[t];
	}
	
	return sum/s->Nx/s->Ny;

}

// set average of PFC array
void set_ave_PFC(struct S* s, double* a, double ave_) {

	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;
	
	double ave0 = ave_PFC(s, a);
	
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			a[ij] += ave_ - ave0;
		}
	}
	
}

// randomize PFC array
void randomize_PFC(struct S* s, double* a, double nave, double namp) {
	
	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;
	
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			a[ij] = nave + namp*(2.0*rand()/RAND_MAX - 1.0);
		}
	}
	
}

// initialize PFC operators
// see update_AB() in pfc-relax.c at github.com/petenez/pfc
void init_operators(struct S* s) {

	int Nxh = s->Nx/2 + 1;
	double _NxNy = 1.0/s->Nx/s->Ny;
	
	double dkx = 2.0*pi/s->Nx/s->dx;
	double dky = 2.0*pi/s->Ny/s->dy;

	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxh*j;
		double ky2;
		if(j < s->Ny/2) ky2 = j*dky;
		else ky2 = (j - s->Ny)*dky;
		ky2 *= ky2;
		for(int i = 0; i < Nxh; ++i) {
			int ij = ij0 + i;
			double kx = i*dkx;
			double k2 = ky2 + kx*kx;
			double k2_1 = k2 - 1.0;
			double a = s->alpha + k2_1*k2_1;
			double ex = exp(-k2*a*s->dt);
			s->A[ij] = ex*_NxNy;
			if(a == 0.0) s->B[ij] = -k2*s->dt*_NxNy;
			else s->B[ij] = (ex - 1.0)/a*_NxNy;
		}
	}

}

// update PFC density and smoothed density
// see step() in pfc-relax.c at github.com/petenez/pfc
void update_PFC(struct S* s) {

	// (*)
	int Nxh = s->Nx/2 + 1;
	int Nxp = 2*Nxh;
	
	if(s->n[0] != s->n[0]) {	// if NaN
		printf("Error: NaN value detected.\n");
		exit(1);
	}
	
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			s->m[ij] = s->n[ij]*s->n[ij]*s->n[ij];
		}
	}
	
	fftw_execute(s->n_N);
	fftw_execute(s->m_M);
	
	double dkx = 2.0*pi/s->Nx/s->dx;
	double dky = 2.0*pi/s->Ny/s->dy;
	double a = -0.5/s->sigma/s->sigma;
	
	fftw_complex* N = (fftw_complex*)s->n;
	fftw_complex* M = (fftw_complex*)s->m;
	
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxh*j;
		double ky2;
		if(j < s->Ny/2) ky2 = j*dky;
		else ky2 = (j - s->Ny)*dky;
		ky2 *= ky2;
		for(int i = 0; i < Nxh; ++i) {
			int ij = ij0 + i;
			double kx = i*dkx;
			double k2 = ky2 + kx*kx;
			
			// Gaussian kernel
			double G = exp(a*k2);
			
			N[ij] = s->A[ij]*N[ij] + s->B[ij]*M[ij];
			
			// smoothed density from convolution
			M[ij] = G*N[ij];
		}
	}
	
	fftw_execute(s->N_n);
	fftw_execute(s->M_m);
	
}

// get norm of an array
double abs1(struct S* s, double* a) {

	// physical array width (padding due to periodic boundaries) (**)
	int Nxp = s->Nx + 2;
	
	int T = omp_get_max_threads();
	double a_[T];
	
	for(int t = 0; t < T; ++t) {
		a_[t] = 0.0;
	}
	
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			a_[omp_get_thread_num()] += a[ij]*a[ij];
		}
	}
	
	double result = 0.0;
	
	for(int t = 0; t < T; ++t) {
		result += a_[t];
	}
	
	return sqrt(result);
	
}

// get norm of difference of two arrays
double abs2(struct S* s, double* a, double* b) {

	// (**)
	int Nxp = s->Nx + 2;
	
	int T = omp_get_max_threads();
	double ab_[T];
	
	for(int t = 0; t < T; ++t) {
		ab_[t] = 0.0;
	}
	
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			double diff = a[ij] - b[ij];
			ab_[omp_get_thread_num()] += diff*diff;
		}
	}
	
	double result = 0.0;
	
	for(int t = 0; t < T; ++t) {
		result += ab_[t];
	}
	
	return sqrt(result);
	
}

// update periodic boundaries
void update_boundaries(struct S* s, double* a) {

	// (**)
	int Nxp = s->Nx + 2;
	
	int NxpNy = Nxp*s->Ny;
	int NxpNy_ = Nxp*(s->Ny + 1);
	
	// horizontal boundaries
	#pragma omp parallel for
	for(int i = 1; i <= s->Nx; ++i) {
		a[i] = a[NxpNy + i];
		a[NxpNy_ + i] = a[Nxp + i];
	}
	
	// vertical boundaries
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int Nxpj = Nxp*j;
		a[Nxpj] = a[Nxpj + s->Nx];
		a[Nxpj + s->Nx + 1] = a[Nxpj + 1];
	}
	
	// corners
	a[0] = a[NxpNy + s->Nx + 1];
	a[Nxp - 1] = a[NxpNy + 1];
	a[NxpNy_] = a[Nxp + s->Nx];
	a[NxpNy_ + s->Nx + 1] = a[Nxp + 1];

}

// get average of array
double ave(struct S* s, double* a) {

	// (**)
	int Nxp = s->Nx + 2;
	
	int T = omp_get_max_threads();
	double a_[T];
	
	for(int t = 0; t < T; ++t) {
		a_[t] = 0.0;
	}
	
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			a_[omp_get_thread_num()] += a[ij];
		}
	}
	
	double sum = 0.0;
	
	for(int t = 0; t < T; ++t) {
		sum += a_[t];
	}
	
	return sum/s->Nx/s->Ny;

}

// set average of array
void set_ave(struct S* s, double* a, double ave_) {

	// (**)
	int Nxp = s->Nx + 2;
	
	double ave0 = ave(s, a);
	
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			a[ij] += ave_ - ave0;
		}
	}
	
	update_boundaries(s, a);
	
}

// randomize array
void randomize(struct S* s, double* a, double nave, double namp) {

	// (**)
	int Nxp = s->Nx + 2;
	
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			a[ij] = nave + namp*(2.0*rand()/RAND_MAX - 1.0);
		}
	}
	
	update_boundaries(s, a);
	
}

// diffuse array
// apply diffusion to array a using Jacobi method
void diffuse(struct S* s, double** a) {

	// (**)
	int Nxp = s->Nx + 2;
	
	double epsilon = 0.01;
	
	// initialize p with a
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			s->p[ij] = (*a)[ij];
		}
	}
	
	update_boundaries(s, s->p);
	
	// iterate Jacobi method until converged
	for(int k = 0; k == 0 || abs2(s, s->p, s->q)/abs1(s, s->p) > epsilon; ++k) {
	
		#pragma omp parallel for
		for(int j = 1; j <= s->Ny; ++j) {
			int ij0 = Nxp*j;
			for(int i = 1; i <= s->Nx; ++i) {
				int ij = ij0 + i;
				s->q[ij] = ((*a)[ij] + s->mu[ij]*(s->p[ij + 1] + s->p[ij - 1] + s->p[ij + Nxp] + s->p[ij - Nxp]))/(1.0 + 4.0*s->mu[ij]);
			}
		}
		
		update_boundaries(s, s->q);
		
		// swap p and q
		double* ptr = s->p;
		s->p = s->q;
		s->q = ptr;
	
	}
	
	// swap a and p
	double* tmp = *a;
	*a = s->p;
	s->p = tmp;

}

// project flow
// apply projection to flow using Jacobi method
void project(struct S* s) {

	// (**)
	int Nxp = s->Nx + 2;
	
	double epsilon = 0.01;
	
	// initialize p and q
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			s->p[ij] = 0.0;
			s->q[ij] = -0.5*(s->vx[ij + 1] - s->vx[ij - 1] + s->vy[ij + Nxp] - s->vy[ij - Nxp]);
		}
	}
	
	update_boundaries(s, s->p);
	update_boundaries(s, s->q);
	
	// iterate Jacobi method until converged
	for(int k = 0; k == 0 || abs2(s, s->p, s->r)/abs1(s, s->p) > epsilon; ++k) {
	
		#pragma omp parallel for
		for(int j = 1; j <= s->Ny; ++j) {
			int ij0 = Nxp*j;
			for(int i = 1; i <= s->Nx; ++i) {
				int ij = ij0 + i;
				s->r[ij] = 0.25*(s->q[ij] + s->p[ij + 1] + s->p[ij - 1] + s->p[ij + Nxp] + s->p[ij - Nxp]);
				
			}
		}
	
		update_boundaries(s, s->r);
		
		// swap p and r
		double* ptr = s->p;
		s->p = s->r;
		s->r = ptr;
		
	}
	
	// update velocity
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			s->vx[ij] -= 0.5*(s->p[ij + 1] - s->p[ij - 1]);
			s->vy[ij] -= 0.5*(s->p[ij + Nxp] - s->p[ij - Nxp]);
		}
	}
	
	update_boundaries(s, s->vx);
	update_boundaries(s, s->vy);
	
}

// advect flow
void advect(struct S* s) {

	// (**)
	int Nxp = s->Nx + 2;
	
	// get flow at each grid point by backtracking and bilinear interpolation 
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			
			// origin of flow
			double x = i - s->vx[ij];
			double y = j - s->vy[ij];
			
			// periodic boundaries
			while(x < 1) x += s->Nx;
			while(x >= s->Nx + 1) x -= s->Nx;
			while(y < 1) y += s->Ny;
			while(y >= s->Ny + 1) y -= s->Ny;
			
			// surrounding grid points' indices
			int i0 = (int)x;
			int j0 = (int)y;
			int i1 = i0 + 1;
			int j1 = j0 + 1;
			
			// interpolate horizontal velocity
			double vx00 = s->vx[Nxp*j0 + i0];
			double vx10 = s->vx[Nxp*j0 + i1];
			double vx01 = s->vx[Nxp*j1 + i0];
			double vx11 = s->vx[Nxp*j1 + i1];
			double vx0 = vx00 + (vx10 - vx00)*(x - i0);
			double vx1 = vx01 + (vx11 - vx01)*(x - i0);
			s->p[ij] = vx0 + (vx1 - vx0)*(y - j0);
			
			// interpolate vertical velocity
			double vy00 = s->vy[Nxp*j0 + i0];
			double vy10 = s->vy[Nxp*j0 + i1];
			double vy01 = s->vy[Nxp*j1 + i0];
			double vy11 = s->vy[Nxp*j1 + i1];
			double vy0 = vy00 + (vy10 - vy00)*(x - i0);
			double vy1 = vy01 + (vy11 - vy01)*(x - i0);
			s->q[ij] = vy0 + (vy1 - vy0)*(y - j0);
		}
	}
	
	update_boundaries(s, s->p);
	update_boundaries(s, s->q);
	
	// swap vx and p, and vy and q
	double* tmp = s->vx;
	s->vx = s->p;
	s->p = tmp;
	tmp = s->vy;
	s->vy = s->q;
	s->q = tmp;

}

// apply gravity
// add upward velocity proportional to smoothed PFC density ("buoyancy")
void gravity(struct S* s) {

	// (*, **)
	int Nxh_PFC = s->Nx/2 + 1;
	int Nxp_PFC = 2*Nxh_PFC;
	int Nxp_flow = s->Nx + 2;
	
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0_PFC = Nxp_PFC*j;
		int ij0_flow = Nxp_flow*(j + 1) + 1;
		for(int i = 0; i < s->Nx; ++i) {
			s->vy[ij0_flow + i] += 0.00005*s->m[ij0_PFC + i];
		}
	}
	
	update_boundaries(s, s->vy);
	
}

// update viscosity
// map viscosity linearly to smoothed PFC density
void update_mu(struct S* s) {
	
	// (*)
	int Nxh_PFC = s->Nx/2 + 1;
	int Nxp_PFC = 2*Nxh_PFC;
	
	double min = 1.0e300;
	double max = -1.0e300;
	
	// determine extrema of smoothed PFC density
	for(int j = 0; j < s->Ny; ++j) {
		int ij0 = Nxp_PFC*j;
		for(int i = 0; i < s->Nx; ++i) {
			int ij = ij0 + i;
			if(s->m[ij] < min) min = s->m[ij];
			if(s->m[ij] > max) max = s->m[ij];
		}
	}
	
	// (**)
	int Nxp_flow = s->Nx + 2;
	
	// constant smoothed density
	if(min == max) {
		#pragma omp parallel for
		for(int j = 0; j < s->Ny; ++j) {
			int ij0_flow = Nxp_flow*(j + 1) + 1;
			for(int i = 0; i < s->Nx; ++i) {
				s->mu[ij0_flow + i] = s->mumin;
			}
		}
	}
	// nonconstant smoothed density
	else {
		double _range = 1.0/(max - min);
		#pragma omp parallel for
		for(int j = 0; j < s->Ny; ++j) {
			int ij0_PFC = Nxp_PFC*j;
			int ij0_flow = Nxp_flow*(j + 1) + 1;
			for(int i = 0; i < s->Nx; ++i) {
				double frac = (s->m[ij0_PFC + i] - min)*_range;
				s->mu[ij0_flow + i] = s->mumin + (s->mumax - s->mumin)*0.5*(1.0 + tanh(100.0*(frac - 0.2)));
			}
		}
	}
	
}

// advect PFC density
// advect PFC density with flow field
void advect_PFC(struct S* s) {

	// (*, **)
	int Nxh_PFC = s->Nx/2 + 1;
	int Nxp_PFC = 2*Nxh_PFC;
	int Nxp_flow = s->Nx + 2;
	
	// initialize q
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0_PFC = Nxp_PFC*j;
		int ij0_flow = Nxp_flow*(j + 1) + 1;
		for(int i = 0; i < s->Nx; ++i) {
			s->q[ij0_flow + i] = s->n[ij0_PFC + i];
		}
	}
	
	update_boundaries(s, s->q);
	
	// see advect()
	#pragma omp parallel for
	for(int j = 1; j <= s->Ny; ++j) {
		int ij0 = Nxp_flow*j;
		for(int i = 1; i <= s->Nx; ++i) {
			int ij = ij0 + i;
			double x = i - s->vx[ij];
			double y = j - s->vy[ij];
			while(x < 1) x += s->Nx;
			while(x >= s->Nx + 1) x -= s->Nx;
			while(y < 1) y += s->Ny;
			while(y >= s->Ny + 1) y -= s->Ny;
			int i0 = (int)x;
			int j0 = (int)y;
			int i1 = i0 + 1;
			int j1 = j0 + 1;
			double a00 = s->q[Nxp_flow*j0 + i0];
			double a10 = s->q[Nxp_flow*j0 + i1];
			double a01 = s->q[Nxp_flow*j1 + i0];
			double a11 = s->q[Nxp_flow*j1 + i1];
			double a0 = a00 + (a10 - a00)*(x - i0);
			double a1 = a01 + (a11 - a01)*(x - i0);
			s->p[ij] = a0 + (a1 - a0)*(y - j0);
		}
	}
	
	// update PFC density
	#pragma omp parallel for
	for(int j = 0; j < s->Ny; ++j) {
		int ij0_PFC = Nxp_PFC*j;
		int ij0_flow = Nxp_flow*(j + 1) + 1;
		for(int i = 0; i < s->Nx; ++i) {
			s->n[ij0_PFC + i] = s->p[ij0_flow + i];
		}
	}

}

// evolve system
void evolve(struct S* s) {
	
	// iterate
	for(s->t = 0; s->t <= s->tend; ++(s->t)) {
		
		// print progress
		if(!(s->t%s->Tprint)) {
			printf("%d\n", s->t);
		}
		
		// write output
		if(!(s->t%s->Twrite)) {
			write_n(s);
			write_m(s);
			write_v(s);
		}
		
		if(s->t < s->tend) {
		
			// multiple PFC density updates per each flow velocity update
			// (more rigid and nicer looking crystals)
			for(int i = 0; i < 10; ++i) {
				update_PFC(s);
			}
			
			// initial average flow velocity and PFC density
			double vx0 = ave(s, s->vx);
			double vy0 = ave(s, s->vy);
			double nave0 = ave_PFC(s, s->n);
			
			// apply gravity and update viscosity
			gravity(s);
			update_mu(s);
			
			// update flow field
			diffuse(s, &s->vx);
			diffuse(s, &s->vy);
			project(s);
			advect(s);
			project(s);
			
			// eliminate velocity drift
			set_ave(s, s->vx, vx0);
			set_ave(s, s->vy, vy0);
			
			// advect PFC density and eliminate drift
			advect_PFC(s);
			set_ave_PFC(s, s->n, nave0);
		}
	}
	
}

// free arrays and destroy FFT plans
void clean(struct S* s) {

	fftw_free(s->A);
	fftw_free(s->B);
	fftw_free(s->n);
	fftw_free(s->m);
	
	fftw_free(s->vx);
	fftw_free(s->vy);
	fftw_free(s->mu);
	fftw_free(s->p);
	fftw_free(s->q);
	fftw_free(s->r);
	
	fftw_destroy_plan(s->n_N);
	fftw_destroy_plan(s->m_M);
	fftw_destroy_plan(s->N_n);
	fftw_destroy_plan(s->M_m);
	
}

int main() {
	
	struct S s;
	
	s.tend = 1000000;		// last iteration
	
	s.Nx = 512;				// system width
	s.Ny = 512;				// system height
	
	s.dx = 0.75;			// system discretization (PFC*)
	s.dy = 0.75;
	s.dt = 0.5;				// time step (*)
	
	s.alpha = -0.6;			// quadratic term coefficient (*)
	s.sigma = 0.15;			// smoothing spread (*)
	
	s.mumax = 1.0e9;		// maximum viscosity
	s.mumin = 1.0e-3;		// minimum viscosity
	
	s.Tprint = 100;			// print interval
	s.Twrite = 100;			// write inteval
	s.pprint = 12;			// print precision
	s.pwrite = 6;			// write precision
	
	double nave0 = -0.5;	// initial average PFC density
	
	// initialize random number generator
	srand(time(NULL));

	// initialize OpenMP threads
	if(!fftw_init_threads()) {
		printf("Error: Initializing threads failed.\n");
		exit(1);
	}
	
	// use maximum number of threads for FFTs
	fftw_plan_with_nthreads(omp_get_max_threads());
	
	// allocate data arrays, set up FFT plans and initialize PFC operators
	allocate_data(&s);
	plan_FFTs(&s);
	init_operators(&s);
	
	// initialize PFC data
	randomize_PFC(&s, s.n, nave0, 1.0);
	set_ave_PFC(&s, s.n, nave0);
	randomize_PFC(&s, s.m, nave0, 0.0);
	
	// initialize flow data
	randomize(&s, s.vx, 0.0, 0.0);
	randomize(&s, s.vy, 0.0, 0.0);
	randomize(&s, s.mu, s.mumin, 0.0);

	// evolve system	
	evolve(&s);
	
	// clean up
	clean(&s);

	return 0;
	
}
