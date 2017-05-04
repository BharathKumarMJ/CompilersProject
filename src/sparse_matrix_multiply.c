#include "mmio.h"
#include "spmm.h"
#include <time.h>

int main(int argc, char **argv) {

	clock_t start, diff;
	int msec;

	struct sparse_matrix_t* A = load_sparse_matrix(MATRIX_MARKET, argv[1]);
	start = clock();
	struct csr_matrix_t* C = csr_matrix_matmatmult ((struct csr_matrix_t*) (A->repr), (struct csr_matrix_t*) (A->repr));
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("\nTime taken %ld seconds %ld milliseconds\n", diff/1000, diff);
	FILE *fp = fopen("output.mtx", "w");
	//print_csr_matrix_in_matrix_market_format(fp, C);
	fclose(fp);


	struct sparse_matrix_t* A1 = load_bcsr_matrix(MATRIX_MARKET, argv[2]);
	start = clock();
	struct bcsr_matrix_t* C1 = bcsr_matrix_matmatmult ((struct bcsr_matrix_t*) (A1->repr), (struct bcsr_matrix_t*) (A1->repr));
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("\nTime taken %ld seconds %ld milliseconds\n", diff/1000, diff);
	//save_bcsr_matrix_in_matrix_market_format("output_bcsr.mtx", C1);
	return 0;
}
