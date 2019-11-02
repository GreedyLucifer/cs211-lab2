#include "../include/for_you_to_do.h"
/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
	register int i, t, j, k, maxind, temps;
	register double max;
	register double *tempv = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < (n - 1); i++)
	{
		// pivoting
		maxind = i;
		max = fabs(A[i * n + i]);

		for (t = i + 1; t < n; t++)
		{
			if (fabs(A[t * n + i]) > max)
			{
				maxind = t;
				max = fabs(A[t * n + i]);
			}
		}
		if (max == 0)
		{
			return -1;
		}
		else
		{
			if (maxind != i)
			{
				// save pivoting information
				temps = ipiv[i];
				ipiv[i] = ipiv[maxind];
				ipiv[maxind] = temps;
				// swap rows
				memcpy(tempv, A + i * n, n * sizeof(double));
				memcpy(A + i * n, A + maxind * n, n * sizeof(double));
				memcpy(A + maxind * n, tempv, n * sizeof(double));
			}
		}

		// factorization
		for (j = i + 1; j < n; j++)
		{
			A[j * n + i] = A[j * n + i] / A[i * n + i];
			for (k = i + 1; k < n; k++)
			{
				A[j * n + k] -= A[j  *n + i] * A[i * n + k];
			}
		}
	}
	free(tempv);

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
	register int i, j;
	register double sum;
	register double *y = (double*)malloc(n * sizeof(double));
	if (UPLO == 'L')
	{
		y[0] = B[ipiv[0]];
		for (i = 1; i < n; i++)
		{
			sum = 0;
			for (j = 0; j < i; j++)
			{
				sum += y[j] * A[i*n + j];
			}
			y[i] = B[ipiv[i]] - sum;
		}
	}
	else if (UPLO == 'U')
	{
		y[n - 1] = B[n - 1] / A[(n - 1)*n + n - 1];
		for (i = n - 2; i >= 0; i--)
		{
			sum = 0;
			for (j = i + 1; j < n; j++)
			{
				sum += y[j] * A[i*n + j];
			}
			y[i] = (B[i] - sum) / A[i*n + i];
		}
	}

	memcpy(B, y, sizeof(double) * n);
	free(y);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(const double *A, const double *B, double *C, int n, int i, int j, int k, int b)
{
	/* add your code here */
	/* please just copy from your lab1 function optimal( ... ) */
	int i1, j1, k1;
		/* B x B mini matrix multiplications */
		for (i1 = i; (i1 < i + b) && (i1 < n); i1 += 3)
			for (j1 = j; (j1 < j + b) && (j1 < n); j1 += 3)
			{
				int c0 = i1 * n + j1;
				int c1 = c0 + n;
				int c2 = c1 + n;
				register double c00 = C[c0];
				register double c01 = C[c0 + 1];
				register double c02 = C[c0 + 2];
				register double c10 = C[c1];
				register double c11 = C[c1 + 1];
				register double c12 = C[c1 + 2];
				register double c20 = C[c2];
				register double c21 = C[c2 + 1];
				register double c22 = C[c2 + 2];

				for (k1 = k; (k1 < k + b) && (k1 < n); k1 += 3)
				{
					int a0 = i1 * n + k1;
					int a1 = a0 + n;
					int a2 = a1 + n;
					int b0 = k1 * n + j1;
					int b1 = b0 + n;
					int b2 = b1 + n;
					register double a00 = A[a0];
					register double a10 = A[a1];
					register double a20 = A[a2];
					register double b00 = B[b0]; register double b01 = B[b0 + 1]; register double b02 = B[b0 + 2];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02;

					a00 = A[a0 + 1];
					a10 = A[a1 + 1];
					a20 = A[a2 + 1];
					b00 = B[b1]; b01 = B[b1 + 1]; b02 = B[b1 + 2];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02;

					a00 = A[a0 + 2];
					a10 = A[a1 + 2];
					a20 = A[a2 + 2];
					b00 = B[b2]; b01 = B[b2 + 1]; b02 = B[b2 + 2];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02;

				}
				C[c0] = c00;
				C[c0 + 1] = c01;
				C[c0 + 2] = c02;
				C[c1] = c10;
				C[c1 + 1] = c11;
				C[c1 + 2] = c12;
				C[c2] = c20;
				C[c2 + 1] = c21;
				C[c2 + 2] = c22;
			}
				/*register int c0 = i1 * n + j1;
				register int c1 = c0 + n;
				register int c2 = c1 + n;
				register int c3 = c2 + n;
				register double c00 = C[c0];
				register double c01 = C[c0 + 1];
				register double c02 = C[c0 + 2];
				register double c03 = C[c0 + 3];
				register double c10 = C[c1];
				register double c11 = C[c1 + 1];
				register double c12 = C[c1 + 2];
				register double c13 = C[c1 + 3];
				register double c20 = C[c2];
				register double c21 = C[c2 + 1];
				register double c22 = C[c2 + 2];
				register double c23 = C[c2 + 3];
				register double c30 = C[c3];
				register double c31 = C[c3 + 1];
				register double c32 = C[c3 + 2];
				register double c33 = C[c3 + 3];

				for (k1 = k; (k1 < k + b) && (k1 < n); k1 += 4)
				{
					register int a0 = i1 * n + k1;
					register int a1 = a0 + n;
					register int a2 = a1 + n;
					register int a3 = a2 + n;
					register int b0 = k1 * n + j1;
					register int b1 = b0 + n;
					register int b2 = b1 + n;
					register int b3 = b2 + n;
					register double a00 = A[a0];
					register double a10 = A[a1];
					register double a20 = A[a2];
					register double a30 = A[a3];
					register double b00 = B[b0]; register double b01 = B[b0 + 1]; register double b02 = B[b0 + 2]; register double b03 = B[b0 + 3];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02; c03 -= a00 * b03;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02; c13 -= a10 * b03;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02; c23 -= a20 * b03;
					c30 -= a30 * b00; c31 -= a30 * b01; c32 -= a30 * b02; c33 -= a30 * b03;

					a00 = A[a0 + 1];
					a10 = A[a1 + 1];
					a20 = A[a2 + 1];
					a30 = A[a3 + 1];
					b00 = B[b1]; b01 = B[b1 + 1]; b02 = B[b1 + 2]; b03 = B[b1 + 3];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02; c03 -= a00 * b03;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02; c13 -= a10 * b03;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02; c23 -= a20 * b03;
					c30 -= a30 * b00; c31 -= a30 * b01; c32 -= a30 * b02; c33 -= a30 * b03;

					a00 = A[a0 + 2];
					a10 = A[a1 + 2];
					a20 = A[a2 + 2];
					a30 = A[a3 + 2];
					b00 = B[b2]; b01 = B[b2 + 1]; b02 = B[b2 + 2]; b03 = B[b2 + 3];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02; c03 -= a00 * b03;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02; c13 -= a10 * b03;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02; c23 -= a20 * b03;
					c30 -= a30 * b00; c31 -= a30 * b01; c32 -= a30 * b02; c33 -= a30 * b03;

					a00 = A[a0 + 3];
					a10 = A[a1 + 3];
					a20 = A[a2 + 3];
					a30 = A[a3 + 3];
					b00 = B[b3]; b01 = B[b3 + 1]; b02 = B[b3 + 2]; b03 = B[b3 + 3];

					c00 -= a00 * b00; c01 -= a00 * b01; c02 -= a00 * b02; c03 -= a00 * b03;
					c10 -= a10 * b00; c11 -= a10 * b01; c12 -= a10 * b02; c13 -= a10 * b03;
					c20 -= a20 * b00; c21 -= a20 * b01; c22 -= a20 * b02; c23 -= a20 * b03;
					c30 -= a30 * b00; c31 -= a30 * b01; c32 -= a30 * b02; c33 -= a30 * b03;

				}
				C[c0] = c00;
				C[c0 + 1] = c01;
				C[c0 + 2] = c02;
				C[c0 + 3] = c03;
				C[c1] = c10;
				C[c1 + 1] = c11;
				C[c1 + 2] = c12;
				C[c1 + 3] = c13;
				C[c2] = c20;
				C[c2 + 1] = c21;
				C[c2 + 2] = c22;
				C[c2 + 3] = c23;
				C[c3] = c30;
				C[c3 + 1] = c31;
				C[c3 + 2] = c32;
				C[c3 + 3] = c33;
			}*/

    return;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
	int ib, i, j, k, maxind;
	double max, sum;
	double *tempv = (double*)malloc(sizeof(double) * n);

	for (ib = 0; ib < (n - 1); ib += b)
	{
		for (i = ib; i < ib + b && i < n; i++)
		{
			// pivoting
			maxind = i;
			max = fabs(A[i*n + i]);

			for (j = i + 1; j < n; j++)
			{
				if (fabs(A[j*n + i]) > max)
				{
					maxind = j;
					max = fabs(A[j*n + i]);
				}
			}
			if (max == 0)
			{
				return -1;
			}
			else
			{
				if (maxind != i)
				{
					// save pivoting information
					int temp = ipiv[i];
					ipiv[i] = ipiv[maxind];
					ipiv[maxind] = temp;
					// swap rows
					memcpy(tempv, A + i * n, n * sizeof(double));
					memcpy(A + i * n, A + maxind * n, n * sizeof(double));
					memcpy(A + maxind * n, tempv, n * sizeof(double));
				}
			}

			// factorization
			for (j = i + 1; j < n; j++)
			{
				A[j*n + i] = A[j*n + i] / A[i*n + i];
				for (k = i + 1; k < ib + b && k < n; k++)
				{
					A[j*n + k] -= A[j*n + i] * A[i*n + k];
				}
			}
		}

		// update A(ib:end, end+1:n)
		for (i = ib; i < ib + b && i < n; i++)
		{
			for (j = ib + b; j < n; j++)
			{
				sum = 0;
				for (k = ib; k < i; k++)
				{
					sum += A[i*n + k] * A[k*n + j];
				}
				A[i*n + j] -= sum;
			}
		}

		// update A(end+1:n, end+1:n)
		for (i = ib + b; i < n; i += b)
		{
			for (j = ib + b; j < n; j += b)
			{
				mydgemm(A, A, A, n, i, j, ib, b);
			}
		}
	}
    return 0;
}

