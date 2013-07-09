
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* (c) 1988-1992 Numerical Recipes Software */
/* With some modifications by Michael Froh */

void spline(float x[], float y[], int n, float yp1, float ypn, float *y2) {
		int i,k;
		float p, qn, sig, un, *u;

		u = (float *)malloc(sizeof(float)*n-1);
		if (yp1 > 0.99e30)
				y2[0]=u[0]=0.0;
		else {
				y2[0] = -0.5;
				u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0]) - yp1);
		}
		for (i=1; i< n-1; i++) {
				sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
				p = sig*y2[i-1]+2.0;
				y2[i] = (sig - 1.0)/p;
				u[i] = (y[i+1] - y[i])/(x[i+1]-x[i]) 
						- (y[i]-y[i-1])/(x[i]-x[i-1]);
				u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
		}
		if (ypn > 0.99e30)
				qn = un = 0;
		else {
				qn = 0.5;
				un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
		}
		y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
		for (k=n-1;k>0;k--)
				y2[k-1]=y2[k-1]*y2[k]+u[k-1];
		free(u);
}

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y) {
		int klo, khi, k;
		float h,b,a;

		klo = 0;
		khi = n-1;
		while (khi-klo > 1) {
				k = (khi+klo+1) >> 1;
				if (xa[k] > x)
						khi = k;
				else
						klo = k;
		}
		h = xa[khi]-xa[klo];
		assert(h != 0.0);
		if (h < 0.0 ) {
				printf("h is negative, x is %f\n",x);
		}
		a = (xa[khi] - x)/h;
		b = (x - xa[klo])/h;
		*y = a * ya[klo] + b*ya[khi]+((a*a*a-a)*y2a[klo] 
						+ (b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
