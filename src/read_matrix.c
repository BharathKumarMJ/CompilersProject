#include "matrix_formats.h"
#include "read_matrix.h"
#include "mmio.h"
#include "merge_sort.h"

int
read_matrix_market_sparse (const char* filename, struct coo_matrix_t* A, struct csr_matrix_t* B)
{
	int ret_code;
	char _matcode[4];
	FILE *f;
	int M, N, nnz;   
	int *II; 
	int *JJ;
	void *val;
	int i;
	enum value_type_t value_type;
  	/* Set to nonzero if the MatrixMarket file specifies that symmetric storage is used. */
	int symmetric_p = 0;
	enum symmetry_type_t symmetry_type = UNSYMMETRIC;
	enum symmetric_storage_location_t symmetric_storage_location = -1;

	char* matcode = (char*) _matcode;

	f = fopen (filename, "r");
	if (f == NULL)
	{
		return MM_COULD_NOT_READ_FILE;
	}

	if (mm_read_banner (f, &matcode) != 0)
	{
		return MM_NO_HEADER;
	}

	if (! mm_is_matrix (matcode))
	{
		return MM_UNSUPPORTED_TYPE;
	}
	else if (! mm_is_sparse (matcode))
	{
		return MM_UNSUPPORTED_TYPE;
	}

	if (mm_is_real (matcode))
		value_type = REAL;
	else
	{
		return MM_UNSUPPORTED_TYPE;
	}

	if (mm_is_general (matcode))
		symmetry_type = UNSYMMETRIC;
	else 
	{
		if (mm_is_symmetric (matcode))
		{
			symmetry_type = SYMMETRIC;
			symmetric_p = 1;
		}
		else if (mm_is_skew (matcode))
		{
			symmetry_type = SKEW_SYMMETRIC;
			symmetric_p = 1;
		}
		else 
		{
			return MM_UNSUPPORTED_TYPE;
		}
	}

  /* 
   * Find out size of sparse matrix.  If the matrix is symmetric, then only 
   * the nonzeros on one side of the matrix are shown.  If we wanted to know
   * how many nonzeros there really are, we would need to examine all of them
   * -- in particular, we would need to know how many nonzeros are along the 
   * diagonal (the diagonal may not be full).  But we don't care, since we
   * are only storing the symmetric part.  What we do need to do is figure 
   * out on what side the nonzeros are stored, whether in the upper triangle
   * or in the lower triangle.
   */
	ret_code = mm_read_mtx_crd_size (f, &M, &N, &nnz);
	if (ret_code != 0)
	{
		return -1;
	}
	if (M < 0 || N < 0 || nnz < 0)
	{
		return -2;
	}

  /* 
   * Reserve memory for matrix.  
   */
	II = (int *) malloc (nnz * sizeof(int));
	JJ = (int *) malloc (nnz * sizeof(int));
	if (value_type == REAL)
		val = (double *) malloc (nnz * sizeof(double));
	if (value_type == REAL)
	{
		for (i = 0; i < nnz; i++)
		{
			double x;
			double* __val = (double*) (val);

			fscanf (f, "%d %d %lg\n", &II[i], &JJ[i], &x);
			__val[i] = x;
	 		II[i]--;  /* adjust from 1-based to 0-based */
			JJ[i]--;
		}
	}
	if (symmetric_p)
	{
      /* Figure out whether the nonzeros are stored in the upper triangle or 
       * the lower triangle.  We assume correctness here.  TODO:  For 
       * robustness, we should have a flag that, if set, requires that all 
       * the stored nonzeros be checked to make sure that they are in the 
       * same triangle. */
		for (i = 0; i < nnz; i++)
		{
			if (II[i] < JJ[i])
			{
				symmetric_storage_location = UPPER_TRIANGLE;
				break;
			}
			else if (II[i] > JJ[i])
			{
				symmetric_storage_location = LOWER_TRIANGLE;
				break;
			}
		}

      /* If the matrix is diagonal (i.e. if i reaches nnz), then we arbitrarily
       * decide that the nonzeros are in the upper triangle (which they are,
       * though they are also in the lower triangle). */
		if (i == nnz)
			symmetric_storage_location = UPPER_TRIANGLE;
	}

	if (f != stdin) 
	{
		if (fclose (f))
		{
			free (II);
			free (JJ);
			free (val);
			A->II = NULL;
			A->JJ = NULL;
			A->val = NULL;
			A->m = 0;
			A->n = 0;
			A->nnz = 0;
			return -3;
		}
	}

	A->II   = II;
	A->JJ   = JJ;
	A->val = val;
	A->m   = M;
	A->n   = N;
	A->nnz = nnz;
	A->symmetry_type = symmetry_type;
	A->symmetric_storage_location = symmetric_storage_location;
	A->value_type = value_type;
	A->index_base = 0;
	A->ownership = LIBRARY_DEALLOCATES;
	A->deallocator = &free;

  // Convert COO to CSR
	coo_to_csr (A, B);
	return 0;
}

int
read_bcsr (const char* filename, struct bcsr_matrix_t* A)
{
	int ret_code;
	FILE *f;
	int bm, bn, r, c, nnzb, num_indices, num_ptrs, data_size;   
	int *colind; 
	int *rowptr;
	void *values;
	int MAX_LINE_SIZE = 1024;
	enum value_type_t value_type;
  	/* Set to nonzero if the MatrixMarket file specifies that symmetric storage is used. */
	int symmetric_p = 0;
	enum symmetry_type_t symmetry_type = UNSYMMETRIC;
	enum symmetric_storage_location_t symmetric_storage_location = -1;

	char line[MAX_LINE_SIZE];

	//printf("Inside read bcsr\n");

	f = fopen (filename, "r");

    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // First segment bn
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &bn);
    }
    //printf("Seg 1 done\n");

    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Second segment bm
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &bm);
    }
    //printf("Seg 2 done\n");

    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Third segment r
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &r);
    }

	//printf("Seg 3 done\n");
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Fourth segment c
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &c);
    }

    //printf("Seg 4 done\n");
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Fifth segment nnzb
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &nnzb);
    }

    //printf("Seg 5 done\n");
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Sixth segment numIndices
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &num_indices);
    }

    //printf("Seg 6 done\n");
    // Allocate num_indices array 
    colind = (int *) malloc(num_indices * sizeof(int));
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Seventh segment Indices
    if (line[0] == '%') {
    	for (int i=0; i<num_indices; i++) {
    		fgets(line, MAX_LINE_SIZE, f);
    		sscanf(line, "%d", &colind[i]);
    	}
    }

    //printf("Seg 7 done\n");
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Eigth segment numPtrs
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &num_ptrs);
    }

    //printf("Seg 8 done\n");
    // Allocate rowptr array 
    rowptr = (int *) malloc(num_ptrs * sizeof(int));
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Ninth segment Indices
    if (line[0] == '%') {
    	for (int i=0; i<num_ptrs; i++) {
    		fgets(line, MAX_LINE_SIZE, f);
    		sscanf(line, "%d", &rowptr[i]);
    	}
    }

    //printf("Seg 9 done\n");
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // Tenth segment datasize
    if (line[0] == '%') {
    	fgets(line, MAX_LINE_SIZE, f);
    	sscanf(line, "%d", &data_size);
    }

    //printf("Seg 10 done\n");
    // Allocate data array 
    values = (double *) malloc(data_size * sizeof(double));
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    // 11 segment data
    if (line[0] == '%') {
    	for (int i=0; i<data_size; i++) {
    		double x;
    		double* __val = (double*) (values);
    		fgets(line, MAX_LINE_SIZE, f);
    		sscanf(line, "%lg", &x);
    		__val[i] = x;
    	}
    }

    //printf("Seg 11 done\n");
    char end[MAX_LINE_SIZE];
    if (fgets(line,MAX_LINE_SIZE,f) == NULL) 
            return MM_PREMATURE_EOF;
    sscanf(line, "%s", end);
    //printf("End is %s\n", end);
        
	A->bn   = bn;
	A->bm   = bm;
	A->r = r;
	A->c   = c;
	A->nnzb   = nnzb;
	A->values = values;
	A->colind = colind;
	A->rowptr = rowptr;
	A->symmetry_type = symmetry_type;
	A->symmetric_storage_location = symmetric_storage_location;
	A->value_type = value_type;
	A->ownership = LIBRARY_DEALLOCATES;
	A->deallocator = &free;
	A->col_oriented_p = 0;

	fclose(f);
	return 0;
}


struct csr_matrix_t*
coo_to_csr (struct coo_matrix_t* A, struct csr_matrix_t* B)
{
	int m   = A->m;
	int n   = A->n;
	int nnz = A->nnz;
	enum symmetry_type_t symmetry_type = A->symmetry_type;
	enum symmetric_storage_location_t symmetric_storage_location = A->symmetric_storage_location;
	enum value_type_t value_type = A->value_type;
	int    *rowptr = NULL;
	int    *colidx = NULL;
	void   *__values = NULL;
	void   *__coord_array = NULL;
	int i, j, currow = 0;
	int index_base = A->index_base;

	if (A->value_type == REAL)
		__values = calloc (nnz, sizeof (double));
	else 
		exit (EXIT_FAILURE);

	rowptr = calloc ((m+1), sizeof (int));
	colidx = calloc (nnz,   sizeof (int));

	if (nnz == 0) 
	{
      /* calloc fills in rowptr with zeros, which is all 
         we need if the matrix is empty */
		init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, UNSYMMETRIC, 
			UPPER_TRIANGLE, value_type, LIBRARY_DEALLOCATES,
			&free, NO_COPY);
		return B; 
	}

  /* Intermediate conversion to coordinate array, so we can sort the entries */
	coo_matrix_to_coord_elem_array (&__coord_array, &nnz, A);

	sort_coord_elem_array_for_csr_conversion (__coord_array, nnz, value_type);

	if (value_type == REAL)
	{
		struct coord_array_t
		{
			int r;
			int c;
			double val;
		};
		struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
		double* values = (double*) __values;
      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * row, then within each row by initial column, now the first row with 
       * nonzeros is coord_array[0].r.
       */
		currow = coord_array[0].r - index_base;
		for (i = 0; i <= currow; i++)
		{
	  /* 
	   * Until we get to first row with an entry in it, all the rowptrs 
	   * before then are zero.  The rowptr for that first column is also 
	   * zero. 
	   */
			rowptr[i] = 0;  
		}
		for (i = 0; i < nnz; i++)
		{
			if (coord_array[i].r - index_base > currow)
			{
	      /* 
	       * We may jump more than one row at a time, so set the rowptr 
	       * entries for the empty rows in between.
	       */
				for (j = currow+1; j <= coord_array[i].r - index_base; j++) 
				{
					if (j - index_base < 0 || j - index_base > m)
					{
						free (values);
						free (rowptr);
						free (colidx);
						free (B);
						return NULL;
					}
					rowptr[j] = i;
				}

				currow = coord_array[i].r - index_base;
			}

			values[i] = coord_array[i].val;   
			colidx[i] = coord_array[i].c - index_base;
		}

      /* Set the last entries in rowptr appropriately */
		for (j = currow+1; j <= m; j++)
			rowptr[j] = nnz;

		init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, symmetry_type, 
			symmetric_storage_location, value_type, 
			LIBRARY_DEALLOCATES, &free, NO_COPY);
	}
	free (__coord_array);
	return B;
}

void
coo_matrix_to_coord_elem_array (void** p_coord_array, 
	int *p_length, 
	const struct coo_matrix_t* A)
{
	const int nnz = A->nnz;
	int k;

	*p_length = nnz;

	if (A->value_type == REAL)
	{
		struct coord_elem_t
		{
			int r;
			int c;
			double val;
		};
		struct coord_elem_t* coord_array = malloc (nnz * sizeof (struct coord_elem_t));
		double* values = (double*) (A->val);

		for (k = 0; k < nnz; k++)
		{
			coord_array[k].r = A->II[k];
			coord_array[k].c = A->JJ[k];
			coord_array[k].val = values[k];
		}
		*p_coord_array = coord_array;
	}
}

/*========================================================================*/
void
sort_coord_elem_array_for_csr_conversion (void* coord_array, 
	const int length,
	enum value_type_t value_type)
{ 
	if (value_type == REAL)
	{
		struct coord_elem_t
		{
			int r;
			int c;
			double val;
		};
		merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
			compare_coord_elem_by_col_real);
		merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
			compare_coord_elem_by_row_real);
	}
}

/**
 * Compares two coordinate-array elements by their row indices, and returns 
 * values like strcmp does.  Comparison function for sorting.
 */
int
compare_coord_elem_by_row_real (const void* a, const void* b)
{
	struct coord_elem_t
	{
		int r;
		int c;
		double val;
	};

	const struct coord_elem_t* x = ((struct coord_elem_t*) a);
	const struct coord_elem_t* y = ((struct coord_elem_t*) b);

	if (x->r < y->r) 
		return -1;
	else if (x->r > y->r)
		return +1;

  /* else */
	return 0;
}

/**
 * Compares two coordinate-array elements by their column indices, and 
 * returns values like strcmp does.  Comparison function for sorting.
 */
int
compare_coord_elem_by_col_real (const void* a, const void* b)
{
	struct coord_elem_t
	{
		int r;
		int c;
		double val;
	};

	const struct coord_elem_t* x = ((struct coord_elem_t*) a);
	const struct coord_elem_t* y = ((struct coord_elem_t*) b);

	if (x->c < y->c) 
		return -1;
	else if (x->c > y->c)
		return +1;

  /* else */
	return 0;
}

/*======================================================================*/
struct csr_matrix_t*
create_csr_matrix (const int m, const int n, const int nnz, 
	void* values, int* colidx, int* rowptr,
	enum symmetry_type_t symmetry_type,
	enum symmetric_storage_location_t symmetric_storage_location,
	enum value_type_t value_type,
	enum ownership_mode_t ownership,
	void (*deallocator) (void*),
	enum copy_mode_t copy_mode)
{
	struct csr_matrix_t *A;

	A = malloc (sizeof (struct csr_matrix_t));
	init_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		symmetric_storage_location, value_type, ownership, 
		deallocator, copy_mode);
	return A;
}

/*=====================================================================*/
void
init_csr_matrix (struct csr_matrix_t* A, 
	const int m, const int n, const int nnz, 
	void* values, int* colidx, int* rowptr,
	enum symmetry_type_t symmetry_type,
	enum symmetric_storage_location_t symmetric_storage_location,
	enum value_type_t value_type,
	enum ownership_mode_t ownership,
	void (*deallocator) (void*),
	enum copy_mode_t copy_mode)
{
	pack_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		symmetric_storage_location, value_type, ownership,
		deallocator, copy_mode);
}

/*=====================================================================*/
void
pack_csr_matrix (struct csr_matrix_t* A, 
	const int m, const int n, const int nnz, 
	void* values, int* colidx, int* rowptr,
	enum symmetry_type_t symmetry_type,
	enum symmetric_storage_location_t symmetric_storage_location,
	enum value_type_t value_type,
	enum ownership_mode_t ownership,
	void (*deallocator) (void*),
	enum copy_mode_t copy_mode)
{
	A->m = m;
	A->n = n;
	A->nnz = nnz;
	if (copy_mode == NO_COPY)
	{
		A->values = values;
		A->colidx = colidx;
		A->rowptr = rowptr;
	}
  else /* copy mode: this matrix gets (deep) copies of the user's input arrays */
	{
      /* In copy mode, the library is responsible for deallocating the arrays */
		if (ownership != LIBRARY_DEALLOCATES)
		{
			A->ownership = LIBRARY_DEALLOCATES;
		}
		if (value_type == REAL) 
		{
			A->values = malloc (nnz * sizeof (double));
			memcpy (A->values, values, nnz * sizeof (double));
		}
 
		A->colidx = malloc (nnz * sizeof (int));
		memcpy (A->colidx, colidx, nnz * sizeof (int));

		A->rowptr = malloc ((n+1) * sizeof (int));
		memcpy (A->rowptr, rowptr, (n+1) * sizeof (int));
	}
	A->symmetry_type = symmetry_type;
	A->symmetric_storage_location = symmetric_storage_location;
	A->value_type = value_type;
	A->ownership = ownership;
	if (deallocator == NULL)
		A->deallocator = &free;
	else
		A->deallocator = deallocator;
}


/*======================================================================*/
struct bcsr_matrix_t*
create_bcsr_matrix (const int bm, const int bn, const int r, const int c, 
		    const int nnzb, void* values, int* colind, int* rowptr,
		    const enum symmetry_type_t symmetry_type,
		    const enum symmetric_storage_location_t symmetric_storage_location, 
		    const enum value_type_t value_type,
		    const int col_oriented_p,
		    enum ownership_mode_t ownership,
		    void (*deallocator) (void*),
		    enum copy_mode_t copy_mode)
{
  struct bcsr_matrix_t* A = malloc (sizeof (struct bcsr_matrix_t));
  init_bcsr_matrix (A, bm, bn, r, c, nnzb, values, colind, rowptr, 
		    symmetry_type, symmetric_storage_location, value_type,
		    col_oriented_p, ownership, deallocator, copy_mode);
  return A;
}

/*======================================================================*/
void
init_bcsr_matrix (struct bcsr_matrix_t* A, const int bm, const int bn, 
		  const int r, const int c, const int nnzb, void* values, 
		  int* colind, int* rowptr, const enum symmetry_type_t symmetry_type,
		  const enum symmetric_storage_location_t symmetric_storage_location,
		  const enum value_type_t value_type, const int col_oriented_p,
		  enum ownership_mode_t ownership,
		  void (*deallocator) (void*),
		  enum copy_mode_t copy_mode)
{
  A->bm = bm;
  A->bn = bn;
  A->r = r;
  A->c = c;
  A->nnzb = nnzb;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  A->col_oriented_p = col_oriented_p;
  A->ownership = ownership;
  if (deallocator == NULL)
    A->deallocator = &free;
  else
    A->deallocator = deallocator;
 
  if (copy_mode == NO_COPY)
    {
      A->values = values;
      A->colind = colind;
      A->rowptr = rowptr;
    }
  else 
    {
      const int nnz = nnzb * r * c;
      /* 
       * Copy mode means that the library takes responsibility
       * for deallocation, using the standard deallocator. 
       */
      A->ownership = LIBRARY_DEALLOCATES;
      A->deallocator = &free;
      if (value_type == REAL)
	{
	  A->values = malloc (nnz * sizeof (double));
	  memcpy (A->values, values, nnz * sizeof (double));
	}
      A->rowptr = malloc ((bm+1) * sizeof (int)); 
      memcpy (A->rowptr, rowptr, (bm+1) * sizeof (int));
      A->colind = malloc (nnzb * sizeof (int)); 
      memcpy (A->colind, colind, nnzb * sizeof (int));
    }
}
