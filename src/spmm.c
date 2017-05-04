#include "spmm.h"
#include <time.h>

struct sparse_matrix_t*
load_sparse_matrix (enum sparse_matrix_file_format_t file_format, 
	const char *const matrix_filename)
{
	struct sparse_matrix_t* A = NULL;
	int errcode = 0;

	if (file_format == MATRIX_MARKET) {
		struct coo_matrix_t* repr = malloc (sizeof (struct coo_matrix_t));
		struct csr_matrix_t* csr_repr = calloc (1, sizeof (struct csr_matrix_t));

		/* Read COO-format matrix from file */
		errcode = read_matrix_market_sparse (matrix_filename, repr, csr_repr);
		if (errcode != 0)
		{
			free (repr);
			free (csr_repr);
			return NULL;
		}
		A = create_sparse_matrix (CSR, csr_repr);
		if ((A == NULL) || (csr_repr == NULL))
		{
			return NULL;
		}
	}

	return A;
}

struct sparse_matrix_t*
load_bcsr_matrix (enum sparse_matrix_file_format_t file_format, 
	const char *const matrix_filename)
{
	struct sparse_matrix_t* A = NULL;
	int errcode = 0;

	if (file_format == MATRIX_MARKET) {
		struct bcsr_matrix_t* repr = malloc (sizeof (struct bcsr_matrix_t));
		errcode = read_bcsr (matrix_filename, repr);
		// printf("Read matrix\n");
		// printf("bn : %d\n", repr->bn);
		// for (int i=0; i<repr->nnzb; i++) {
		// 	printf("Index : %d\n", repr->colind[i]);
		// }
		A = create_sparse_matrix (BCSR, repr);
	}
	return A;
}

struct sparse_matrix_t*
create_sparse_matrix (enum sparse_matrix_storage_format_t format, 
	void* repr)
{
	struct sparse_matrix_t* A = NULL;

	A = malloc (sizeof (struct sparse_matrix_t));
	A->format = format;
	A->repr = repr;
	return A;
}

struct csr_matrix_t*
csr_matrix_matmatmult (struct csr_matrix_t* B, struct csr_matrix_t* A)
{
  /* B is m x p and A is p x n. */
	const int m = B->m;
	const int p = B->n;
	const int n = A->n;
	const enum value_type_t value_type = B->value_type;
	int errcode = 0;
	// printf("Before anything\n");
	if (p != A->m)
	{
		printf("csr_matrix_matmatmult :: Returning NULL\n");
		return NULL;
	}
	else if (B->value_type != A->value_type)
	{
		printf("csr_matrix_matmatmult :: Returning value type NULL\n");
		return NULL;
	}

	double* Cval = NULL;
	int* Cptr = NULL; 
	int* Cind = NULL;
	int nnz = 0;
	double* Bvalues = (double*) B->values;
	double* Avalues = (double*) A->values;
	// printf("Before mul\n");
	errcode = csr_matmatmult_double_real (&Cptr, &Cind, &Cval, &nnz, 1.0, 
		A->rowptr, A->colidx, Avalues,
		B->rowptr, B->colidx, Bvalues,
		m, p, n);
	// printf("After mul\n");
	if (errcode != 0)
	{
		return NULL;
	}

	if (m > 0)
	{
		struct csr_matrix_t* C = create_csr_matrix (m, n, nnz, Cval, Cind, Cptr, 
			UNSYMMETRIC, 0, REAL, 
			LIBRARY_DEALLOCATES, NULL, 
			NO_COPY);
		return C;
	}
}


int
csr_matmatmult_double_real (int** pCptr, int** pCind, double** pCval, 
	int* pCnnz, double alpha, int* Aptr, 
	int* Aind, double* Aval, int* Bptr, 
	int* Bind, double* Bval, 
	const int m, const int p, const int n)
{
  /* Borrowed from the Hypre code hypre_CSRMatrixMultiply */
	int ia, ib, ic, ja, jb, num_nonzeros=0;
	int row_start, counter;
	double a_entry, b_entry;
	int* B_marker;
	int* Cptr;
	int* Cind;
	double* Cval;

	B_marker = calloc (n, sizeof (int));
	Cptr = malloc ((m+1) * sizeof (int));
	for (ic = 1; ic < m; ic++)
    Cptr[ic] = -1; /* flag to detect errors */

		Cptr[0] = 0;

	for (ib = 0; ib < n; ib++)
		B_marker[ib] = -1;

  /* 
   * Count the number of nonzeros in each row of C, and use this
   * to set up the ptr array of C.
   */
  for (ic = 0; ic < m; ic++) /* for each row of C */
	{
      /*// fprintf (3, "\t\tRow %d of C\n", ic);*/
      /* For each row ic of A:
       *   For each entry (ic,ja) in row ic of A:
       *     Look in row ja of B (since entry (ic,ja) of A weights row ja of B):
       *     For each entry (ja,jb) in row ja of B:
       *       If we haven't seen an entry in column jb of B for the current row ic of A
       *         Mark it and increment nnz
       *       EndIf
       *     EndFor
       *   EndFor
       * EndFor
       */
		for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++) 
		{
			ja = Aind[ia];
			if (ja < 0 || ja >= p)
			{
	      // fprintf (0, "*** csr_matmatmult_double_real: Element A(%d,%d) is "
			//"out of range [0,%d] x [0,%d] ***\n", ic, ja, m-1, p-1);
				free (Cptr);
				free (B_marker);
				return -1;
			}
			for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
			{
				jb = Bind[ib];
				if (jb < 0 || jb >= n)
				{
		  // fprintf (0, "*** csr_matmatmult_double_real: Index jb = %d "
			    //"is out of range [0,%d] ***\n", jb, n-1);
					free (Cptr);
					free (B_marker);
					return -1;
				}
				if (B_marker[jb] != ic)
				{
					B_marker[jb] = ic;
					num_nonzeros++;
		  /*// fprintf (4, "\t\t\tIdentified new nonzero (%d,%d)\n", ic, jb);*/
				}
			}
		}
		Cptr[ic+1] = num_nonzeros;            
	}
  // fprintf (3, "\tThere are %d nonzeros in C.\n", num_nonzeros);
  // fprintf (3, "\tAllocating Cval and Cind\n");
	Cval = malloc (num_nonzeros * sizeof (double));
	Cind = malloc (num_nonzeros * sizeof (int));

	for (ic = 0; ic < num_nonzeros; ic++)
	{
      /* 
       * We set Cval[ic] and Cind[ic] to "invalid" values (NaN and -1, 
       * respectively) to make sure that we fill them all in.   Note 
       * that we can't use 0.0/0.0 to generate NaN because when running
       * this function via an FFI under some Lisps (SBCL in particular),
       * 0.0/0.0 raises a floating-point exception instead of being a 
       * silent NaN like we want it to be (mfh 28 July 2006).  That's 
       * why we need the C99 #defined value NAN and the fpclassify()
       * function in <math.h>.  We use the assert to make sure that the 
       * NAN really is a NAN.
       */
#ifdef USE_ISNAN
		Cval[ic] = NAN;
      // assert (fpclassify (Cval[ic]) == FP_NAN);
#else
      /* 
       * We don't have C99 so we can't use NAN; initialize Cval[ic]
       * to some innocuous value instead.  Note that we can't just set
       * it to 0.0/0.0, as this may be evaluated; under plain C this 
       * typically doesn't signal, but under Lisp running the C function,
       * it may signal.
       */
		Cval[ic] = 0.0;
#endif
		Cind[ic] = -1;
	}

	for (ib = 0; ib < n; ib++)
		B_marker[ib] = -1;

  // fprintf (3, "\tComputing nonzero values in C\n");
	counter = 0;
	clock_t start, diff;
	start = clock();
	for (ic = 0; ic < m; ic++)
	{
		row_start = Cptr[ic];
		for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++)
		{
			ja = Aind[ia];
			a_entry = Aval[ia];
			for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
			{
				jb = Bind[ib];
				b_entry = Bval[ib];
				if (B_marker[jb] < row_start)
				{
					B_marker[jb] = counter;
					Cind[B_marker[jb]] = jb;
					Cval[B_marker[jb]] = a_entry * b_entry;
					counter++;
				}
				else
				{
					Cval[B_marker[jb]] += a_entry * b_entry;
				}
			}
		}
	}
	diff = clock() - start;
	printf("\nMULTIPLY Time taken %ld seconds %ld milliseconds\n", diff/1000, diff);
	free (B_marker);
  /* Go back through and scale by alpha */
  // fprintf (3, "\tScaling by alpha\n");
	for (ic = 0; ic < m; ic++)
	{
		int k;
		for (k = Cptr[ic]; k < Cptr[ic+1]; k++)
			Cval[k] = alpha * Cval[k];
	}

	*pCnnz = num_nonzeros;
	*pCval = Cval;
	*pCptr = Cptr;
	*pCind = Cind;
  // fprintf (2, "=== Done with csr_matmatmult_double_real ===\n");
	return 0;
}  

struct bcsr_matrix_t*
bcsr_matrix_matmatmult (struct bcsr_matrix_t* B, struct bcsr_matrix_t* A)
{
  /* B is m x p and A is p x n. */
	const int bm = B->bm;
	const int bp = B->bn;
	const int bn = A->bn;
	const enum value_type_t value_type = B->value_type;
	int errcode = 0;
	if (bp != A->bm)
	{
		printf("bcsr_matrix_matmatmult :: Returning NULL\n");
		return NULL;
	}
	else if (B->value_type != A->value_type)
	{
		printf("bcsr_matrix_matmatmult :: Returning value type NULL\n");
		return NULL;
	}

	double* Cval = NULL;
	int* Cptr = NULL; 
	int* Cind = NULL;
	int nnzb = 0;
	double* Bvalues = (double*) B->values;
	double* Avalues = (double*) A->values;
	errcode = bcsr_matmatmult_double_real (&Cptr, &Cind, &Cval, &nnzb, 1.0, 
		A->rowptr, A->colind, Avalues,
		B->rowptr, B->colind, Bvalues,
		bm, bp, bn, A->r, A->c);
	if (errcode != 0)
	{
		return NULL;
	}

	if (bm > 0)
	{
		struct bcsr_matrix_t* C = create_bcsr_matrix (bm, bn, A->r, A->c, nnzb, Cval, Cind, Cptr,
			     A->symmetry_type, A->symmetric_storage_location, 
			     A->value_type, A->col_oriented_p,
			     LIBRARY_DEALLOCATES, &free, NO_COPY);
		return C;
	}
}

int
bcsr_matmatmult_double_real (int** cRowptr, int** cColind, double** cValues, int* cNnzb, double alpha, 
	int* aRowptr, int* aColind, double* aValues,
	int* bRowptr, int* bColind, double* bValues,
	const int mb, const int pb, const int nb, const int r, const int c)
{
  /* Borrowed from the Hypre code hypre_CSRMatrixMultiply */
	int ia, ib, ic, ja, jb, num_nonzerosb=0;
	int row_start, counter;
	double *a_entry, *b_entry;
	int* B_marker;
	int* Cptr;
	int* Cind;
	double* Cval;

	B_marker = calloc (nb, sizeof (int));
	Cptr = malloc ((mb+1) * sizeof (int));
	for (ic = 1; ic < mb; ic++)
    	Cptr[ic] = -1; /* flag to detect errors */

	Cptr[0] = 0;

	for (ib = 0; ib < nb; ib++)
		B_marker[ib] = -1;

  /* 
   * Count the number of nonzeros in each row of C, and use this
   * to set up the ptr array of C.
   */
  for (ic = 0; ic < mb; ic++) /* for each row of C */
	{
      /*// fprintf (3, "\t\tRow %d of C\n", ic);*/
      /* For each row ic of A:
       *   For each entry (ic,ja) in row ic of A:
       *     Look in row ja of B (since entry (ic,ja) of A weights row ja of B):
       *     For each entry (ja,jb) in row ja of B:
       *       If we haven't seen an entry in column jb of B for the current row ic of A
       *         Mark it and increment nnz
       *       EndIf
       *     EndFor
       *   EndFor
       * EndFor
       */
		for (ia = aRowptr[ic]; ia < aRowptr[ic+1]; ia++) 
		{
			ja = aColind[ia];
			if (ja < 0 || ja >= pb)
			{
				free (Cptr);
				free (B_marker);
				return -1;
			}
			for (ib = bRowptr[ja]; ib < bRowptr[ja+1]; ib++)
			{
				jb = bColind[ib];
				if (jb < 0 || jb >= nb)
				{
					free (Cptr);
					free (B_marker);
					return -1;
				}
				if (B_marker[jb] != ic)
				{
					B_marker[jb] = ic;
					num_nonzerosb++;
				}
			}
		}
		Cptr[ic+1] = num_nonzerosb;            
	}

	Cval = calloc (num_nonzerosb * r * c, sizeof (double));
	Cind = malloc (num_nonzerosb * sizeof (int));

	for (ic = 0; ic < num_nonzerosb; ic++)
	{
		Cind[ic] = -1;
	}
       
//        * We set Cval[ic] and Cind[ic] to "invalid" values (NaN and -1, 
//        * respectively) to make sure that we fill them all in.   Note 
//        * that we can't use 0.0/0.0 to generate NaN because when running
//        * this function via an FFI under some Lisps (SBCL in particular),
//        * 0.0/0.0 raises a floating-point exception instead of being a 
//        * silent NaN like we want it to be (mfh 28 July 2006).  That's 
//        * why we need the C99 #defined value NAN and the fpclassify()
//        * function in <math.h>.  We use the assert to make sure that the 
//        * NAN really is a NAN.
       
// #ifdef USE_ISNAN
// 		Cval[ic] = NAN;
//       // assert (fpclassify (Cval[ic]) == FP_NAN);
// #else
//       /* 
//        * We don't have C99 so we can't use NAN; initialize Cval[ic]
//        * to some innocuous value instead.  Note that we can't just set
//        * it to 0.0/0.0, as this may be evaluated; under plain C this 
//        * typically doesn't signal, but under Lisp running the C function,
//        * it may signal.
//        */
// 		Cval[ic] = 0.0;
// #endif
		

	for (ib = 0; ib < nb; ib++)
		B_marker[ib] = -1;


	clock_t start, diff;
	start = clock();
	counter = 0;
	for (ic = 0; ic < mb; ic++)
	{
		row_start = Cptr[ic];
		for (ia = aRowptr[ic]; ia < aRowptr[ic+1]; ia++)
		{
			ja = aColind[ia];
			a_entry = &aValues[ia*r*c];
			for (ib = bRowptr[ja]; ib < bRowptr[ja+1]; ib++)
			{
				jb = bColind[ib];
				b_entry = &bValues[ib*r*c];
				if (B_marker[jb] < row_start)
				{
					B_marker[jb] = counter;
					Cind[B_marker[jb]] = jb;
					//Cval[B_marker[jb]] = a_entry * b_entry;
					for (int rr=0; rr<c; rr++) {
						for (int pp=0; pp<r; pp++) {
							for (int qq=0; qq<c; qq++) {
								Cval[B_marker[jb]*r*c + pp * r + qq] += a_entry[pp * r + rr] * b_entry[rr * r + qq];
								// if (counter == 0) {
								// 	printf("Value : %lg\t", Cval[B_marker[jb] + pp * r + qq]);
								// 	printf("Index : %d\n", B_marker[jb] + pp * r + qq);
								// 	printf("A value : %lg\n", a_entry[pp * r + rr]);
								// 	printf("B value : %lg\n", b_entry[rr * r + qq]);
								// 	//printf("Bmarker : %d\n", B_marker[jb]);
								// }
							}
						}
					}
					counter++;
				}
				else
				{
					for (int rr=0; rr<c; rr++) {
						for (int pp=0; pp<r; pp++) {
							for (int qq=0; qq<c; qq++) {
								Cval[B_marker[jb]*r*c + pp * r + qq] += a_entry[pp * r + rr] * b_entry[rr * r + qq];
							}
						}
					}
				}
			}
		}
	}
	diff = clock() - start;
	printf("\nMULTIPLY Time taken %ld seconds %ld milliseconds\n", diff/1000, diff);
	free (B_marker);

	*cNnzb = num_nonzerosb;
	*cValues = Cval;
	*cRowptr = Cptr;
	*cColind = Cind;
	return 0;
}


int
print_csr_matrix_in_matrix_market_format (FILE* out, const struct csr_matrix_t* A)
{
	int start, end;
	int i, k;
	char value_type_label[20];
	char symmetry_type_label[20];

	if (A->value_type == REAL)
		strncpy (value_type_label, "real", 19);

	if (A->symmetry_type == UNSYMMETRIC)
		strncpy (symmetry_type_label, "general", 19);
	else if (A->symmetry_type == SYMMETRIC)
		strncpy (symmetry_type_label, "symmetric", 19);
	else if (A->symmetry_type == SKEW_SYMMETRIC)
		strncpy (symmetry_type_label, "skew-symmetric", 19);
	else
		return -1;

	fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", 
		value_type_label, symmetry_type_label);

	fprintf (out, "%d %d %d\n", A->m, A->n, A->nnz);

	if (A->value_type == REAL)
	{
		const double* const values = (const double* const) (A->values);

		for (i = 0; i < A->m; i++)
		{
			start = A->rowptr[i];
			end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
			for (k = start; k < end; k++)
				fprintf (out, "%d %d %.13e\n", i+1, A->colidx[k] + 1, values[k]);
		}
	}
	return 0;
}

int 
save_bcsr_matrix_in_matrix_market_format (const char* const filename, 
					  struct bcsr_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  // bebop_log (2, "=== save_bcsr_matrix_in_matrix_market_format ===\n");
  out = fopen (filename, "w");
  if (out == NULL)
    {
      // bebop_log (0, "*** save_bcsr_matrix_in_matrix_market_format: failed "
	       // "to open output file \"%s\" ***\n", filename);
      return -1;
    }

  errcode = print_bcsr_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      // bebop_log (0, "*** save_bcsr_matrix_in_matrix_market_format: failed "
	       // "to close output file \"%s\" ***\n", filename);
      return -1;
    }
  // bebop_log (2, "=== Done with save_bcsr_matrix_in_matrix_market_format ===\n");
  return errcode;
}

/*======================================================================*/
int
print_bcsr_matrix_in_matrix_market_format (FILE* out, const struct bcsr_matrix_t* A)
{
  int start, end;
  int i, j, k;
  const int col_oriented_p = A->col_oriented_p;
  const int nnzb           = A->nnzb;
  const int bm             = A->bm;
  const int bn             = A->bn;
  const int r              = A->r;
  const int c              = A->c;
  const int block_size     = r * c;
  const int nnz            = nnzb * block_size;
  const int OUTPUT_INDEX_BASE = 1;
  char symmetry_type_label[20];
  char value_type_label[20];

  // bebop_log (2, "=== print_bcsr_matrix_in_matrix_market_format ===\n");

  if (A->symmetry_type == UNSYMMETRIC)
    strncpy (symmetry_type_label, "general", 19);
  else if (A->symmetry_type == SYMMETRIC)
    strncpy (symmetry_type_label, "symmetric", 19);
  else if (A->symmetry_type == SKEW_SYMMETRIC)
    strncpy (symmetry_type_label, "skew-symmetric", 19);
  else 
    {
      // bebop_log (0, "*** print_bcsr_matrix_in_matrix_market_format: "
	       // "Invalid symmetry type %d of the given bcsr_matrix_t sparse "
	       // "matrix! ***\n", A->symmetry_type);
      // bebop_log (2, "=== Done with print_bcsr_matrix_in_matrix_market_format ===\n");
      return -1;
    }

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else
    {
      // bebop_log (0, "*** print_bcsr_matrix_in_matrix_market_format: Uns"
	       // "upported value type! ***\n");
      // bebop_log (2, "=== Done with print_bcsr_matrix_in_matr"
		// "ix_market_format ===\n");
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", value_type_label, symmetry_type_label);
  fprintf (out, "%d %d %d\n", bm * r, bn * c, nnz);

  if (block_size == 1)
    {
      for (i = 0; i < bm; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  for (j = start; j < end; j++)
	    {
	      if (A->value_type == REAL)
		{
		  double* values = (double*) (A->values);
		  fprintf (out, "%d %d %.13e\n", 
			   OUTPUT_INDEX_BASE + i*r,
			   OUTPUT_INDEX_BASE + A->colind[j],
			   values[j * block_size]);
		}
	    }
	}
    }
  else if (col_oriented_p)
    {
      for (i = 0; i < A->bm; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  for (j = start; j < end; j++)
	    {
	      for (k = 0; k < block_size; k++) /* k goes down columns */
		{
		  if (A->value_type == REAL)
		    {
		      double* values = (double*) (A->values);
		      fprintf (out, "%d %d %.13e\n", 
			       OUTPUT_INDEX_BASE + i*r + (k % r), 
			       OUTPUT_INDEX_BASE + A->colind[j] + (k / r), 
			       values[j*block_size + k]);
		    }
		}
	    }
	}
    }
  else /* Nonzero blocks in the values array are row-oriented */
    {
      for (i = 0; i < A->bm; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  for (j = start; j < end; j++)
	    {
	      for (k = 0; k < r*c; k++)  /* k goes across rows */
		{
		  if (A->value_type == REAL)
		    {
		      double* values = (double*) (A->values);
		      fprintf (out, "%d %d %.13e\n", 
			       OUTPUT_INDEX_BASE + i*r + (k / c), 
			       OUTPUT_INDEX_BASE + A->colind[j] + (k % c), 
			       values[j*block_size + k]);
		    }
		}
	    }
	}
    }
  return 0;
}
