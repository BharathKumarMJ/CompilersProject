#include "mmio.h"
#include "spmm.h"
#include <time.h>

extern clock_t csr_mul_time;
extern clock_t bcsr_mul_time;

int main(int argc, char **argv) {

	clock_t start, diff;
	clock_t msec;
	int numIters = 3;
	int i;
	printf("Num Iters: %d\n", numIters);
	struct sparse_matrix_t* A = load_sparse_matrix(MATRIX_MARKET, argv[1]);
	printf("loaded matrix\n");
	start = clock();
	struct csr_matrix_t* C;
	for(i = 0; i < numIters; i++)
		C = csr_matrix_matmatmult ((struct csr_matrix_t*) (A->repr), (struct csr_matrix_t*) (A->repr));
	diff = clock() - start;
	msec = (diff * 1000 / CLOCKS_PER_SEC)/numIters;
	printf("\nCSR MULTIPLY TIME: %ld seconds %ld milliseconds\n", (csr_mul_time/1000)/numIters, (csr_mul_time)/numIters);
	printf("\nCSR Time taken %ld seconds %ld milliseconds\n", msec/1000, msec);
	// FILE *fp = fopen("output.mtx", "w");
	//print_csr_matrix_in_matrix_market_format(fp, C);
	// fclose(fp);

	struct sparse_matrix_t* A1 = load_bcsr_matrix(MATRIX_MARKET, argv[2]);
	start = clock();
	struct bcsr_matrix_t* C1;
	for(i = 0; i < numIters; i++)
		C1 = bcsr_matrix_matmatmult ((struct bcsr_matrix_t*) (A1->repr), (struct bcsr_matrix_t*) (A1->repr));
	diff = clock() - start;
	msec = (diff * 1000 / CLOCKS_PER_SEC)/numIters;
	printf("\nBCSR MULTIPLY TIME: %ld seconds %ld milliseconds\n", (bcsr_mul_time/1000)/numIters, (bcsr_mul_time)/numIters);
	printf("\nBCSR Time taken %ld seconds %ld milliseconds\n", msec/1000, msec);
	//save_bcsr_matrix_in_matrix_market_format("output_bcsr.mtx", C1);
	return 0;
}
