/*
**  Author: Guillermo Alicea
**  03.04.17
**  UCF COT4500 - Professor Glinos
**  Assignment 2: LU Factorization
*/

#include <stdio.h>
#include <stdlib.h>

#define WIDTH 7

int widthU = 0, widthL = 0;

double absolute(double num);

int findDigits(double num);

void findU(double **matrix, double **L, int rows);

void updateRow(double *row1, double *row2, double multiple, int rows);

int main(int argc, char **argv)
{
  //filename not found or read priviledges not granted for current user -> error
  FILE *fp;
  if ((fp = fopen(argv[1], "r")) == NULL)
  {
    printf("\nError detected: input file does not exist, or cannot be read. Exiting..\n\n");
    exit(0);
  }

  int i = 0, j = 0, rows = 0;

  //scan number of rows/columns (guaranteed square matrix)
  if ((fscanf(fp, "%d", &rows)) != 1)
  {
      printf("\nError detected: input file has improper format. Exiting..\n\n");
      exit(0);
  }

  double **matrix = (double **)malloc(sizeof(double) * rows);
  double **L = (double **)malloc(sizeof(double) * rows);

  for (i = 0; i < rows; i++)
  {
    matrix[i] = (double *)malloc(sizeof(double) * rows);
    L[i] = (double *)malloc(sizeof(double) * rows);
  }

  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < rows; j++)
    {
      //number of rows scanned and number of elements parsed do not match -> error
      if ((fscanf(fp, "%lf", &matrix[i][j])) != 1)
      {
        printf("\n\nError detected during input file parsing. Exiting..\n\n");
        exit(0);
      }
      widthU = (widthU > absolute(matrix[i][j])) ? widthU : absolute(matrix[i][j]);
    }
  }

  fclose(fp);

  widthU = findDigits(widthU) + WIDTH;

  printf("----------------------------------------\n\n");
  printf("LU Factorization by Guillermo Alicea\n\n");

  printf("Input Matrix:\n\n");
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < rows; j++)
      printf("%*.4lf ", widthU, matrix[i][j]);
    printf("\n");
  }
  widthU = 0;

  findU(matrix, L, rows);

  widthL = findDigits(widthL) + WIDTH;

  printf("\n\"L\" Matrix:\n\n");
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < rows; j++)
      printf("%*.4lf ", widthL, L[i][j]);
    printf("\n");
  }

  widthU = findDigits(widthU) + WIDTH;

  printf("\n\"U\" Matrix:\n\n");
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < rows; j++)
      printf("%*.4lf ", widthU, matrix[i][j]);
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < rows; i++)
  {
    free(matrix[i]);
    free(L[i]);
  }
  free(matrix);
  free(L);

  return 0;
}

double absolute(double num)
{
  return (num > 0) ? num : -num;
}

int findDigits(double num)
{
  int i = 0;
  int mag = (int) abs(num);
  while (mag > 0)
  {
    mag /= 10;
    i++;
  }
  return i;
}

void updateRow(double *row1, double *row2, double multiple, int rows)
{
  int i = 0;

  for (; i < rows; i++)
  {
    row2[i] = row2[i] - multiple * row1[i];
    widthU = (widthU > absolute(row2[i])) ? widthL : absolute(row2[i]);
  }
}

void findU(double **matrix, double **L, int rows)
{
  int i = 0, j = 0;
  double multiple = 0;

  for (i = 0; i < rows; i++)
    for (j = i + 1; j < rows; j++)
    {
      //matrix[i][i] == 0, matrix cannot be factorized since 1) cannot cancel out the jth row unless
      //all of [j][i]'s below are zero and 2) U will not be upper triangular since an upper element is 0 -> error
      if (matrix[i][i] == 0)
      {
        printf("\nError found during factorization: Division by Zero. Exiting..\n\n");
        exit(0);
      }
      multiple = matrix[j][i] / matrix[i][i];
      updateRow(matrix[i], matrix[j], multiple, rows);
      L[j][i] = multiple;
      widthL = (widthL > absolute(L[j][i])) ? widthL : absolute(L[j][i]);
      L[i][i] = 1;
    }

  L[rows - 1][rows - 1] = 1;
}
