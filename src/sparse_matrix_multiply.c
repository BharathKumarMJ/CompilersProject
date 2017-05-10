#include "mmio.h"
#include "spmm.h"
#include <time.h>
#include <unistd.h>

extern clock_t csr_mul_time;
extern clock_t bcsr_mul_time;

int main(int argc, char **argv) {

	int flags = 0, opt, options = 1;

	while ((opt = getopt(argc, argv, "o:")) != -1) {
	   switch (opt) {
	   case 'o':
	       flags = 1;
	       options = 2;
	       break;
	   default: /* '?' */
	       fprintf(stderr, "Usage: %s [-o] <csr_filename.mtx> <bcsr_filename_from_script>\n",
	               argv[0]);
	       return 0;
	   }
	}

	if (argc < 3) {
		fprintf(stderr, "Usage: %s <csr_filename.mtx> <bcsr_filename_from_script>\n",
	               argv[0]);
		return 0;
	}

	clock_t start, diff;
	clock_t msec;
	int numIters = 3;
	int i;
	
	// Process the CSR matrix and perform SpMM
	struct sparse_matrix_t* A = load_sparse_matrix(MATRIX_MARKET, argv[options]);
	start = clock();
	struct csr_matrix_t* C;
	for(i = 0; i < numIters; i++)
		C = csr_matrix_matmatmult ((struct csr_matrix_t*) (A->repr), (struct csr_matrix_t*) (A->repr));
	diff = clock() - start;
	msec = (diff * 1000 / CLOCKS_PER_SEC)/numIters;
	printf("\nCSR Time taken %ld seconds %ld milliseconds\n", msec/1000, msec);

	// Process the BCSR matrix and perform blocked SpMM
	struct sparse_matrix_t* A1 = load_bcsr_matrix(MATRIX_MARKET, argv[options+1]);
	start = clock();
	struct bcsr_matrix_t* C1;
	for(i = 0; i < numIters; i++)
		C1 = bcsr_matrix_matmatmult ((struct bcsr_matrix_t*) (A1->repr), (struct bcsr_matrix_t*) (A1->repr));
	diff = clock() - start;
	msec = (diff * 1000 / CLOCKS_PER_SEC)/numIters;
	printf("\nBCSR Time taken %ld seconds %ld milliseconds\n", msec/1000, msec);

	if (flags) {
		FILE *fp = fopen("output.mtx", "w");
		print_csr_matrix_in_matrix_market_format(fp, C);
		fclose(fp);	
		save_bcsr_matrix_in_matrix_market_format("output_bcsr.mtx", C1);
	}
	return 0;
}
