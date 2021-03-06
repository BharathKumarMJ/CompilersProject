#ifndef _read_mm_h
#define _read_mm_h
/**
 * @file read_mm.h
 * @author Mark Hoemmen
 * @since 06/10/04 09:48:37 PDT
 * @date Time-stamp: <2008-07-16 10:56:55 mhoemmen>
 * 
 * Functions for loading and converting MatrixMarket - format sparse matrix 
 * files.
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

struct coord_elem_t; /* forward declaration */
struct coo_matrix_t; /* forward declaration */
struct csr_matrix_t; /* forward declaration */


/**
 * Converts a coordinate (COO) format sparse matrix to coordinate array 
 * (array-of-structs) format.  Allocates memory for the array of structs.
 * Caller is responsible for freeing that memory when done.
 *
 * @param p_coord_array [OUT]  Pointer to sparse matrix in ``coordinate 
 *                             array'' format 
 * @param p_length      [OUT]  Pointer to number of entries in the above
 *                             coordinate array
 * @param A             [IN]   The COO-format sparse matrix
 */
void
coo_matrix_to_coord_elem_array (void** p_coord_array, 
				int *p_length, 
				const struct coo_matrix_t* A);

struct coo_matrix_t*
csr_to_coo (struct csr_matrix_t* A);

/**
 * Reads a real (non-complex) sparse matrix (general, without any sort of 
 * implicit symmetry) from a Matrix Market (v. 2.0) file, and stores it in 
 * coordinate (COO) format.  For details on the Matrix Market format, 
 * see <a href="http://math.nist.gov/MatrixMarket"> http://math.nist.gov/MatrixMarket </a>.
 *
 * @param filename [IN]  Name of the MatrixMarket file
 * @param A       [OUT]  The coordinate format matrix
 *
 * @return Error code (zero if no errors)
 *
 * @note Matrix Market files are always 1-based, i.e. the index of the 
 *       first element of a matrix is (1,1), not (0,0) as in C.  We adjust 
 *       the offsets accordingly for C-style storage, when reading from and 
 *       writing to these files.
 */
int
read_matrix_market_sparse (const char* filename, struct coo_matrix_t* A, struct csr_matrix_t* B);

/**
 * Reads a real (non-complex) dense matrix from a Matrix Market (v. 2.0) file,
 * allocates space for it, and stores it in column-oriented format in A.
 *
 * @param filename [IN]  Name of the MatrixMarket file
 * @param m       [OUT]  Number of rows in the matrix
 * @param n       [OUT]  Number of columns in the matrix
 * @param A       [OUT]  The dense matrix, in column-oriented format
 *
 * @return Error code (zero if no errors)
 *
 * @note One can also store a vector as a dense matrix in Matrix Market 
 * format, and read it in using this function.
 */
int
read_matrix_market_real_general_dense (const char* filename, 
				       int *m, int *n, double** A);

struct csr_matrix_t*
coo_to_csr (struct coo_matrix_t* A, struct csr_matrix_t* B);

void
coo_matrix_to_coord_elem_array (void** p_coord_array, 
				int *p_length, 
				const struct coo_matrix_t* A);

void
sort_coord_elem_array_for_csr_conversion (void* coord_array, 
					  const int length,
					  enum value_type_t value_type);

int
compare_coord_elem_by_row_real (const void* a, const void* b);

int
compare_coord_elem_by_col_real (const void* a, const void* b);

struct csr_matrix_t*
create_csr_matrix (const int m, const int n, const int nnz, 
		   void* values, int* colidx, int* rowptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type,
		   enum ownership_mode_t ownership,
		   void (*deallocator) (void*),
		   enum copy_mode_t copy_mode);

void
init_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode);

void
pack_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode);


int
read_bcsr (const char* filename, struct bcsr_matrix_t* A);

struct bcsr_matrix_t*
create_bcsr_matrix (const int bm, const int bn, const int r, const int c, 
		    const int nnzb, void* values, int* colind, int* rowptr,
		    const enum symmetry_type_t symmetry_type,
		    const enum symmetric_storage_location_t symmetric_storage_location, 
		    const enum value_type_t value_type,
		    const int col_oriented_p,
		    enum ownership_mode_t ownership,
		    void (*deallocator) (void*),
		    enum copy_mode_t copy_mode);

void
init_bcsr_matrix (struct bcsr_matrix_t* A, const int bm, const int bn, 
		  const int r, const int c, const int nnzb, void* values, 
		  int* colind, int* rowptr, const enum symmetry_type_t symmetry_type,
		  const enum symmetric_storage_location_t symmetric_storage_location,
		  const enum value_type_t value_type, const int col_oriented_p,
		  enum ownership_mode_t ownership,
		  void (*deallocator) (void*),
		  enum copy_mode_t copy_mode);


#endif /* _read_mm_h */
