#include <stdio.h>
#include <math.h>

void dgefa(int lda, int n, double a[lda][n], int ipvt[], int info); 

void dgesl(int lda, int n, double a[lda][n], int ipvt[], double b[], int job);

void daxpy(int n, double da, double dx[], int incx, double dy[], int incy);

void dscal(int n, double da, double dx[], int incx);

void dmxpy(int n1, double y[], int n2, int ldm, double x[], double m[ldm][n2]);

double ddot(int n, double dx[], int incx, double dy[], int incy);

void tranpose(int n, double a[n][n]);

int idamax(int n, double dx[], int incx);

void aprint(char* sre, int n, double a[n]);

int main() {
	int i, i1, i2, j, job, k, m, n1; 
	const int n = 9;
	const int ndp = 400;
	int ipvt[n], info;
	double d1, dc[n], dc1[n], dv[n], dv1[n];
	double c[n], r[n], x[n][n], v[n], v1[n], t1, t2, t3, t4;

	job = 0;

	// set data - fortran line 17
	c[0] = 5.0;
	c[1] = 2304.0;
	c[2] = 118101.0;
	c[3] = 1838336.0;
	c[4] = 14855109.0;
	c[5] = 79514880.0;
	c[6] = 321537749.0;
	c[7] = 1062287616.0;
	c[8] = 3014530821.0;

	for (i1 = 0; i1 < n; i1++) {
		dc[i1] = c[i1];
	}

	aprint("input data =", 9, dc);

	// calculate polynomial regression coefficients
	for (i2 = 1; i2 < n+1; i2++) {
		for (i1 = 1; i1 < n+1; i1++) {
			t1 = 0.0;

			for (k = 0; k < n; k++) {
				t2 = k;
				if (i1 == 1 && i2 == 1) {
					t1 = t1 + 1.0;
				} else {
					t1 = t1 + pow(t2, i1+i2-2);
				}
			}

			x[i1-1][i2-1] = t1;
		}
	}

	for (i1 = 1; i1 < n+1; i1++) {
		t1 = 0.0;

		for (k = 0; k < n; k++) {
			t2 = k;
			if (i1 == 1) { 
				t1 = t1 + c[k];
			} else {
				t1 = t1 + pow(t2, i1-1)*c[k];
			}
		}

		v[i1-1] = t1;
	}

	// fortran line 70-71
	dgefa(n, n, x, ipvt, info);
	dgesl(n, n, x, ipvt, v, job);

	t1 = 0.0;

	for (i1 = 1; i1 < n+1; i1++) {
		v1[i1-1] = rint(v[i1-1]);
		t1 = fmax(t1, fabs(v[i1-1]-v1[i1-1]));
		dv[i1-1] = v[i1-1];
		dv1[i1-1] = v1[i1-1];
	}

	printf("max deviation from integer value = %g \n", (double) t1);
	aprint("raw regression coefficients =", 9, dv);
	aprint("rounded regression coefficients =", 9, dv1);

	// fortran line 86
	d1 = 0.0;

	for (k = 0; k < n; k++) {
		t1 = 0.0;
		t2 = (double) k;

		for (i1 = 1; i1 < n+1; i1++) {
			if (k == 0 && i1 == 1) {
				t1 = t1 + dv1[i1-1];
			} else {
				t1 = t1 + dv1[i1-1] * pow(t2, i1 - 1);
			}
		}

		dc1[k] = t1;
		d1 = d1 + fabs(dc1[k] - dc[k]);
	}

	// fortran line 104
	aprint("regenerated data based on rounded regression coefficients =", 9, dc1);
	printf("total error from original data = %g \n", d1);
}

// fortran line 111
void dgefa(int lda, int n, double a[lda][n], int ipvt[], int info) {
	double t;
	int j, k, kp1, l, nm1;

	info = 0;
	nm1 = n - 1;
	if (nm1 < 1) goto l70;

	for (k = 1; k < nm1+1; k++) {
		kp1 = k + 1;
		 // find l = pivot index

		// fortran line 176
		l = idamax(n-k+1, &a[k-1][k-1], 1) + k - 1;
		ipvt[k-1] = l;

		// zero pivot implies this column already triangularized
		if (a[l-1][k-1] == 0.0) goto l40;

		// interchange if neccessary
		if (l == k) goto l10;
		t = a[k-1][l-1];
		a[k-1][l-1] = a[k-1][k-1];
		a[k-1][k-1] = t;

l10:
		// compute multipliers
		t = -1.0/a[k-1][k-1];
		// DEBUG
		dscal(n-k,t, &a[k-1][k],1);

		// row elimination with column indexing
		for (j = kp1; j < n+1; j++) {
			t = a[j-1][l-1];
			if (l == k) goto l20;
			a[j-1][l-1] = a[j-1][k-1];
			a[j-1][k-1] = t;
l20:
			daxpy(n-k,t,&a[k-1][k],1,&a[j-1][k],1);
		} // l30
		goto l50;
l40:
		info = k;
l50:
		continue;
	} // l60

l70:
	ipvt[n-1] = n;
	if (a[n-1][n-1] == 0.0) info = n;
	return;
}

void dgesl(int lda, int n, double a[lda][n], int ipvt[], double b[], int job) {
	double t;
	int k, kb, l, nm1;

	nm1 = n - 1;
	if (job != 0) goto l50;

	// job = 0, solve a * x = b
	// first solve l*y = b

	if (nm1 < 1) goto l30;
	for (k = 1; k < nm1+1; k++) {
		l = ipvt[k-1];
		t = b[l-1];
		if (l == k) goto l10;
		b[l-1] = b[k-1];
		b[k-1] = t;
l10:
		daxpy(n-k, t, &a[k-1][k], 1, &b[k], 1);
	} // l20
l30:
	// now solve u*x = y
	
	for (kb = 1; kb < n + 1; kb++) {
		k = n + 1 - kb;
		b[k-1] = b[k-1]/a[k-1][k-1];
		t = -b[k-1];
		daxpy(k-1,t,&a[k-1][0],1,&b[0],1);
	} // l40
	goto l100;

l50: 

	// job = nonzero, solve trans(a) * x = b
	// first solve trans(u)*y = b

	for (k = 1; k < n+1; k++) {
		t = ddot(k-1, &a[k-1][0], 1, &b[0], 1);
		b[k-1] = (b[k-1]-t)/a[k-1][k-1];
	} // l60

	// now solve trans(l)*x = y
	
	if (nm1 < 1) goto l90;
	// fortran line 322
	for (kb = 1; kb < nm1+1; kb++) {
		k = n - kb;
		b[k-1] = b[k-1] + ddot(n-k,&a[k-1][k],1,&b[k],1);
		l = ipvt[k-1];
		if (l == k) goto l70;
		t = b[l-1];
		b[l-1] = b[k-1];
		b[k-1] = t;
l70: 
		continue;
	} // l80

l90:
l100:
	return;
}

// fortran line 336
void daxpy(int n, double da, double dx[], int incx, double dy[], int incy){
	// constant times a vector plus a vector
	// jack dongarra, linkpack, 3/11/78
	
	int i, ix, iy, m, mp1;

	if (n < 0) return;
	if (da == 0.0) return;
	if (incx == 1 && incy == 1) goto l20;

	// code for unequal increments or equal increments
	// not equal to 1
	
	ix = 1;
	iy = 1;
	if (incx < 0) ix = (-n+1)*incx + 1;
	if (incy < 0) iy = (-n+1)*incy + 1;
	for (i = 1; i < n + 1; i++) {
		dy[iy-1] = dy[iy-1] + da*dx[ix-1];
		ix = ix + incx;
		iy = iy + incy;
	} // l10
	return;

	// code for both increments equal to 1
l20:
	for (i = 1; i < n + 1; i ++) {
		dy[i-1] = dy[i-1] + da*dx[i-1];
	} // l30

	return;
}

double ddot(int n, double dx[], int incx, double dy[], int incy) {
	double res, dtemp;
	int i, ix, iy, m, mp1;

	res = 0.0;
	dtemp = 0.0;

	if (n < 0) return res;
	if (incx == 1 && incy == 1) goto l20;

	// code for unequal increments or equal increments
	// not equal to 1
	
	ix = 1;
	iy = 1;
	if (incx < 0) ix = (-n+1)*incx + 1;
	if (incy < 0) iy = (-n+1)*incy + 1;
	for (i = 1; i < n+1; i++) {
		dtemp = dtemp + dx[ix-1]*dy[iy-1];
		ix = ix + incx;
		iy = iy + incy;
	} // l10
	res = dtemp;
	return res;

	// code for both increments equal to 1
	
l20:
	for (i = 1; i < n + 1; i++) {
		dtemp = dtemp + dx[i-1]*dy[i-1];
	} // l30
	res = dtemp;
	return res;
}

// fortran line 414
void dscal(int n, double da, double dx[], int incx) {
	// scales a vector by a constant
	int i, m, mp1, nincx;

	if (n < 0) return;
	if (incx == 1) goto l20;

	// code for increment not equal to 1
	nincx = n*incx;
	for (i = 1; i < nincx + 1; i += incx) {
		dx[i-1] = da*dx[i-1];
	} // l10
	return;

	// code for increment equal to 1
l20:
	for (i = 1; i < n + 1; i ++) {
		dx[i-1] = da*dx[i-1];
	} // l30

	return;
}

int idamax(int n, double dx[], int incx){
	// finds the index of element having max.
	// absolute value
	// jack dongara, linkpack, 3/11/78
	double dmax;
	int res;
	int i, ix;
	// aprint("idamax dx =", n, dx);

	res = 0;
	if (n < 1) return res;
	res = 1;
	if (n == 1) return res;
	if (incx == 1) goto l20;

	// code for increment not equal to 1
	ix = 1;
	dmax = abs(dx[1-1]);
	ix = ix + incx;
	for (i = 2; i < n+1; i++) {
		if (abs(dx[ix-1]) < dmax) goto l5;
		res = i;
		dmax = abs(dx[ix-1]);
l5:
		ix = ix + incx;
	} // l10
	return res;

	// code for increment equal to 1
l20:
	dmax = fabs(dx[1-1]);
	for (i = 2; i < n+1; i++) {
		if (fabs(dx[i-1]) < dmax) goto l30;
		res = i;
		dmax = fabs(dx[i-1]);

l30:
		continue;
	} // l30

	return res;
}

void dmxpy(int n1, double y[], int n2, int ldm, double x[], double m[ldm][n2]) {
	double s;
	int j, k;

	for (k = 1; k < n1 + 1; k++) {
		s = 0.0;

		for (j = 1; j < n2+1; j++) {
			s = s + m[k][j]*x[j];
		} // l100

		y[k] = y[k] + s;
	} // l200

	return;
}

void aprint(char* str, int n, double a[n]) {
	printf("%s\n", str);
	int i;
	for (i = 0; i < n; i+=4) {
		printf("\t%g", a[i]);
		if (i+1 < n) printf("\t%g", a[i+1]);
		if (i+2 < n) printf("\t%g", a[i+2]);
		if (i+3 < n) printf("\t%g", a[i+3]);
		printf("\n");
	}
}

void tranpose(int n, double a[n][n]) {
	int i;
	int j;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			double save = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = save;
		}
	}
}
