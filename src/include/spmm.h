#include "matrix_formats.h"
#include "read_matrix.h"
#include <string.h>
#include <time.h>


struct sparse_matrix_t*
load_sparse_matrix (enum sparse_matrix_file_format_t file_format, 
		    const char *const matrix_filename);

struct sparse_matrix_t*
create_sparse_matrix (enum sparse_matrix_storage_format_t format, 
		      void* repr);

struct csr_matrix_t*
csr_matrix_matmatmult (struct csr_matrix_t* B, struct csr_matrix_t* A);

int
csr_matmatmult_double_real (int** pCptr, int** pCind, double** pCval, 
			    int* pCnnz, double alpha, int* Aptr, 
			    int* Aind, double* Aval, int* Bptr, 
			    int* Bind, double* Bval, 
			    const int m, const int p, const int n);
int
print_csr_matrix_in_matrix_market_format (FILE* out, const struct csr_matrix_t* A);


struct sparse_matrix_t*
load_bcsr_matrix (enum sparse_matrix_file_format_t file_format, 
	const char *const matrix_filename);

struct bcsr_matrix_t*
bcsr_matrix_matmatmult (struct bcsr_matrix_t* B, struct bcsr_matrix_t* A);

int
bcsr_matmatmult_double_real (int** cRowptr, int** cColind, double** cValues, int* cNnzb, double alpha, 
	int* aRowptr, int* aColind, double* aValues,
	int* bRowptr, int* bColind, double* bValues,
	const int mb, const int pb, const int nb, const int r, const int c);

int 
save_bcsr_matrix_in_matrix_market_format (const char* const filename, 
					  struct bcsr_matrix_t* A);

int
print_bcsr_matrix_in_matrix_market_format (FILE* out, const struct bcsr_matrix_t* A);