/*
**  Author: Guillermo Alicea
**  31.01.17
**  UCF COT4500 - Professor Glinos
**  Assignment 1: Polynomial Root Finder using Newton's Method and Secant Method
*/

#include <stdio.h>
#include <stdlib.h>

double computeFunction(double *polynomial, double x, double pN_2, double pN_1, int degree, int mode);

double power(double x, int power);

double absolute(double x);

void doNewton(double p0, double tol, int max, double* polynomial, int degree);

void doSecant(double p0, double p1, double tol, int max, double* polynomial, int degree);

int main(int argc, char **argv)
{
	int     i = 0;
  int     j = 0;
  int     max = 0;
  int     lim = 0;
  double  p0 = 0;
  double  p1 = 0;
  double  tol = 0;
  double *polynomial = NULL;

  //read in command-line arguments, ensuring they are there and of correct count
	if (argc > 6)
	{
      p0 = atof(argv[1]);
      p1 = atof(argv[2]);
      tol = atof(argv[3]);
      max = atoi(argv[4]);
      lim = atoi(argv[5]);
      polynomial = malloc(sizeof(double) * (lim + 1));

		for (i = 6, j = 0; i < argc; i++, j++)
	     polynomial[j] = atof(argv[i]);

    if (j != (lim + 1))
    {
        printf("Incorrect number of command-line arguments.\n"
               "Try again. <p0> <p1> <tol> <max> <lim> <{n_lim,..,n_0}>\n\n");
        exit(1);
    }
	}
  else
  {
    printf("Failed to enter command-line arguments.\n"
           "Try again. <p0> <p1> <tol> <max> <lim> <{n_lim,..,n_0}>\n\n");
    exit(1);
  }

  printf("\nPolynomial Root Finder by Guillermo Alicea\n\n"
        "Input Parameters:\n\n"
        "\tp0 = %.1f\n"
        "\tp1 = %.1f\n"
        "\ttol = %.1e\n"
        "\tmax = %d\n\n", p0, p1, tol, max);

  printf("Polynomial is of order: %d\n\n", lim);

  printf("Terms of polynomial: ");

  if (lim >= 0)
    for (i = 0, j = lim; i <= lim; i++, j--)
      printf("%.1f*x^%d%s", polynomial[i], (j), (i == lim) ? "\n\n": " + ");

  printf("Newton's Method:\n\n");
  doNewton(p0, tol, max, polynomial, lim);

  printf("Secant Method:\n\n");
  doSecant(p0, p1, tol, max, polynomial, lim);

  free(polynomial);

  return 0;
}

//mode == 0: compute secant f' approximation
//mode == 1: compute f(x)
//mode == 2(rather, anything other than 0 or 1): compute f'(x)
double computeFunction(double *polynomial, double x, double pN_2, double pN_1, int degree, int mode)
{
  int    i = 0;
  int    temp = degree;
  double result = 0;

  if (!mode)
  {
    return ((computeFunction(polynomial, pN_1, 0, 0, degree, 1) - computeFunction(polynomial, pN_2, 0, 0, degree, 1)) /
            (pN_1 - pN_2));
  }

  for (i = 0; i <= degree; i++)
  {
    result  += polynomial[i] * ((mode == 1) ? power(x, temp--) : (temp-- * power(x, temp)));
  }

  return result;
}

//Assumes y will always be positive or 0 (constraint of the problem)
double power(double x, int y)
{
  if (y <= 0) return 1;
  return (x * power(x, y - 1));
}

double absolute(double x)
{
  return (x > 0) ? x : -x;
}

//straight from the pseudo code on the lecture slide
void doNewton(double p0, double tol, int max, double* polynomial, int degree)
{
  int    i = 1;
  int    success = 0;
  int    dbz = 0; //divide by zero
  double pN_1 = p0; //p_(n-1)
  double pN = 0; //p_(n)
  double result = 0;

  while (i <= max)
  {
    if ((result = computeFunction(polynomial, pN_1, 0, 0, degree, 2)))
    {
      pN = pN_1 - (computeFunction(polynomial, pN_1, 0, 0, degree, 1) / result);

      if (absolute(pN - pN_1) < tol)
      {
        success = 1;
        printf("\tp%d = %.17lg\n\n", i, pN);
        break;
      }

      printf("\tp%d = %.17lg\n", i, pN);

      i++;
      pN_1 = pN;
    }
    else
    {
      dbz = 1;
      break;
    }
  }

  if (success == 0 && dbz)
    printf("\tSolution could not be found: Divide by Zero detected\n\n");
  else if (success == 0)
    printf("\tSolution could not be found: Maximum iterations reached\n\n");
  else
    printf("\tSolution found after %d iterations: %.17lg\n\n", i, pN);

}

//Also straight from the pseudo code on the lecture slides
void doSecant(double p0, double p1, double tol, int max, double* polynomial, int degree)
{
  int    i = 2;
  int    success = 0;
  int    dbz = 0; //divide by zero
  double pN_2 = p0; //p_(n-2)
  double pN_1 = p1; //p_(n-1)
  double pN = 0; //p_(n)
  double temp = 0;
  double result = 0;

  while (i <= max)
  {
    if ((result = computeFunction(polynomial, 0, pN_2, pN_1, degree, 0)))
    {
      double temp = (pN_1 - pN_2);

      pN = pN_1 - (computeFunction(polynomial, pN_1, 0, 0, degree, 1) * temp) / (result * temp);

      if (absolute(pN - pN_1) < tol)
      {
        success = 1;
        printf("\tp%d = %.17lg\n\n", i, pN);
        break;
      }

      printf("\tp%d = %.17lg\n", i, pN);

      i++;
      pN_2 = pN_1;
      pN_1 = pN;
    }
    else
    {
      dbz = 1;
      break;
    }
  }

  if (success == 0 && dbz)
    printf("\tSolution could not be found: Divide by Zero detected\n\n");
  else if (success == 0)
    printf("\tSolution could not be found: Maximum iterations reached\n\n");
  else
    printf("\tSolution found after %d iterations: %.17lg\n\n", (i - 1), pN);

}
