#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi 6.28318530717959
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
//////////////////////////////////////////////////////////////////////////////////

	double b=529.610525;     //kT=0.00094409//
	double wb=0.006835371;
	double n=0.29;
	double r=0.2*wb;

///////////////////////////////////////////////////
double funj(double w)
	{
	double j;

	j=n*w*wb*wb*wb*wb/((wb*wb-w*w)*(wb*wb-w*w)+4*w*w*r*r);

	return j;
	}

//////////////////////////////////////////////////
double fungr(double t)
	{
	double w,gr=0;

	for(w=0.0001;w<=0.11;w=w+0.0001)
		{
		double valj=funj(w)*(1-cos(w*t))*0.0001/(w*w*tanh(b*w));
		if(isnan(valj))
			{
			printf("Got a w=%lf\n",w);
			exit(0);
			}

		gr=gr+valj;
		}
	return gr;
	}

/////////////////////////////////////////////////
double fungi(double t)
	{
	double w,gi=0;

	for(w=0.0001;w<=0.11;w=w+0.0001)
	{
		gi=gi+funj(w)*(sin(w*t)-w*t)*0.0001/(w*w);
	}

	return gi;
	}
////////////////////////////////////////////////
void swap(double *a,double *b)
{
	double tempr;

	tempr=*a;
	*a=*b;
	*b=tempr;
}
/////////////////////////////////////////////////
int main(void)
{
	int i,m,j;
	long int k=32768;
	double x,w,t=-k/2;
	double data[2*k];


	for(i=0;i<2*k;i+=2)
	{
		x=i/2;

			data[i]=exp(-fungr(t+x))*cos(fungi(t+x));

			data[i+1]=exp(-fungr(t+x))*sin(fungi(t+x));


	}


	for(j=0;j<k;j++)
			{

				swap(&data[j],&data[j+k]);
			}

	four1 (data-1, k, 1);

	for(j=0;j<k;j++)
		{

			swap(&data[j],&data[j+k]);
		}



	for(m=-300;m<220;m+=2)

	{
		w=27.21*Pi*(m/2)/k;


		printf("%lf  %lf\n",w,data[m+k]);

	}


	return 0;
}
