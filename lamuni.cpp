#include <iostream>
#include <cmath>
#include "mex.h"
// #include "mexAdapter.hpp"
using namespace std;

// subfunctions
double stumpffs(double z)
{
	double s = 0;
	if (z > 0)
	{
		s = (sqrt(z) - sin(sqrt(z))) / pow(sqrt(z), 3);
	}

	else if (z < 0)
	{
		s = (sinh(sqrt(-z)) - sqrt(-z)) / pow(sqrt(-z), 3);
	}

	else
	{
		s = 1 / 6;
	}

	return s;
}

double stumpffc(double z)
{
	double c = 0;
	if (z > 0)
	{
		c = (1 - cos(sqrt(z))) / z;
	}

	else if (z < 0)
	{
		c = (cosh(sqrt(-z)) - 1) / (-z);
	}

	else
	{
		c = 1 / 2;
	}

	return c;
}

double yfun(double z, double r1, double r2, double s, double c, double A)
{
	double yf = r1 + r2 + A * (z * s - 1) / sqrt(c);
	return yf;
}

double ffun(double yf, double c, double s, double A, double mu, double dt)
{
	double ff = (pow((yf / c), 1.5) * s) + A * sqrt(yf) - sqrt(mu) * dt;
	return ff;
}

double dFdz(double z, double yf0, double yf, double c, double s, double A)
{
	double dF = 0;

	if (z == 0)
	{
		dF = sqrt(2) / 40 * pow(yf0,1.5) + A / 8 * (sqrt(yf0)) + A * sqrt(1 / 2 / yf0);
	}

	else
	{
		dF = pow((yf / c), 1.5) * (1 / 2 / z * (c - 3 * s / 2 / c) + 3 * pow(s,2) / 4 / c) + A / 8 * (3 * s / c * sqrt(yf) + A * sqrt(c / yf));
	}
	return dF;
}

// main start
	void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
	*prhs[]) {

	int R1[3], R2[3], i;
	char comp[3] = {'x', 'y', 'z'};
	double mu, dt;

	// Requesting R1 and R2 components to fill arrays
	for (i = 0; i < 3; i++)
	{
		cout << "Enter value for component " << comp[i] << " (R1): ";
		cin >> R1[i];
	}

	for (i = 0; i < 3; i++)
	{
		cout << "Enter value for component " << comp[i] << " (R2): ";
		cin >> R2[i];
	}
	
	// Requesting body mu and dt
	cout << "Enter value for mu: ";
	cin >> mu;
	cout << "Enter value for dt: ";
	cin >> dt;
	
	double r1, r2;
	double sum;

	// Finding magnitudes of R1 and R2
	sum = (pow(R1[0],2)) + (pow(R1[1],2)) + (pow(R1[2],2));
	r1 = sqrt(sum);
	cout << "\nr1 is " << r1 << "\n";

	sum = (pow(R2[0],2)) + (pow(R2[1],2)) + (pow(R2[2],2));
	r2 = sqrt(sum);
	cout << "r2 is " << r2 << "\n\n";

	// Doing cross product to find orthagonal vector
	double orvec;
	orvec = (R1[0] * R2[1]) - (R1[1] * R2[0]);
	
	// Doing dot product to use in phi formula
	double dotp = 0;
	for (int i = 0; i < 3; i++)
	{
		dotp = dotp + R1[i] * R2[i];
	}

	double phi = acos((dotp) / (r1 * r2));
	double pi = 2 * acos(0.0); // creating pi

	// Prograde orbit assumption
	if (orvec <= 0)
	{
		phi = 2 * pi - phi;
	}

	double A = sin(phi) * sqrt((r1 * r2) / (1 - cos(phi)));

	double z;
	z = -100; // initial guess
	double s, c, yf, yf0, ff, dF;

	// Calling subfunctions
	s = stumpffs(z);
	c = stumpffc(z);
	yf = yfun(z, r1, r2, s, c, A);
	yf0 = yfun(0, r1, r2, s, c, A);
	ff = ffun(yf, c, s, A, mu, dt);
	dF = dFdz(z, yf0, yf, c, s, A);

	// Finding where ff = 0 
	while (ff < 0)
	{
		z = z + 0.1;
		s = stumpffs(z);
		c = stumpffc(z);
		yf = yfun(z, r1, r2, s, c, A);
		ff = ffun(yf, c, s, A, mu, dt);
		dF = dFdz(z, yf0, yf, c, s, A);
	}

	// tolerance and maximum amount of iterations
	double tol = 1e-8;
	int nMax = 5000;

	// Finding z fixed point
	double ratio = 1.0;
	int n = 0;

	while ((abs(ratio) > tol) && (n <= nMax))
	{
		n = n + 1;
		ratio = ff / dF;
	    z = z - ratio;
		s = stumpffs(z);
		c = stumpffc(z);
		yf = yfun(z, r1, r2, s, c, A);
		yf0 = yfun(0, r1, r2, s, c, A);
		ff = ffun(yf, c, s, A, mu, dt);
		dF = dFdz(z, yf0, yf, c, s, A);
	}

	if (n >= nMax)
	{
		printf("Iteration limit reached.\nRatio = %lf\nz = %lf\n", ratio, z);
	}

	// Lagrange coefficients
	double f = 1 - yf / r1;
	double g = A * sqrt(yf / mu);
	double gdot = 1 - yf / r2;


	// Calculating and printing V1Plus and V2Minus
    double V1Plus[3];
    double V2Minus[3];

    for (i = 0; i <= 2; i++)
    {
    	V1Plus[i] = 1 / g * (R2[i] - f * R1[i]);
    	V2Minus[i] = 1 / g * (gdot * R2[i] - R1[i]);
    }

    printf("V1Plus: \n");
    for (i = 0; i <= 2; i++)
        {
        	printf("%lf \t",V1Plus[i]);
        }
    printf("\n\n");

    printf("V2Minus: \n");
    for (i = 0; i <= 2; i++)
            {
            	printf("%lf \t",V2Minus[i]);
            }

	// return 0;
	}

