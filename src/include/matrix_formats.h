#include <stdlib.h>
#include <stdio.h>
enum
index_base_t
{
  ZERO, ONE
};

enum
symmetry_type_t
{
  UNSYMMETRIC, SYMMETRIC, SKEW_SYMMETRIC
};

enum
symmetric_storage_location_t
{
  LOWER_TRIANGLE, UPPER_TRIANGLE
};

enum
value_type_t
{
  REAL, COMPLEX, PATTERN
};

enum ownership_mode_t 
{
  LIBRARY_DEALLOCATES
};

enum copy_mode_t 
{
  NO_COPY
};

/**
 * The supported file formats for sparse matrices.
 */
enum
sparse_matrix_file_format_t
{
  HARWELL_BOEING, MATRIX_MARKET, MATLAB
};

enum sparse_matrix_storage_format_t {
  BCSR, CSR, COO
};


/**
 * Sparse matrix in Block Compressed Sparse Row (BCSR) format.
 */
struct
bcsr_matrix_t
{
  /** 
   * Number of block rows in the matrix (actual number of rows is r*bm) 
   */
  int bm;
  
  /** 
   * Number of block columns in the matrix (actual number of columns is
   * c*bn)
   */
  int bn;

  /**
   * Number of rows in each (dense) subblock
   */
  int r;

  /**
   * Number of columns in each (dense) subblock
   */
  int c;
  
  /** 
   * Number of stored (nonzero) blocks in the matrix.  If the matrix is 
   * stored in a symmetric (or skew, etc.) format, nnz only refers to the 
   * number of stored entries, not the actual number of nonzeros. 
   */
  int nnzb;

  /** Array of stored (nonzero) entries of the matrix */
  void* values;

  /** 
   * Array of column indices of the upper left corners of the nonzero blocks 
   * of the matrix.  These are the actual column indices and not the "block 
   * indices". 
   */
  int* colind;

  /** Array of indices into the colind and values arrays, for each row */
  int* rowptr;

  /**
   * Symmetry type of the matrix.
   */
  enum symmetry_type_t symmetry_type;

  /**
   * If the matrix has a kind of symmetry (or skew-symmetry): Where the actual 
   * elements of the matrix are stored: in the lower triangle or the upper 
   * triangle.
   */
  enum symmetric_storage_location_t symmetric_storage_location;

  /** 
   * Indicates the type of the entries in val:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the "val" array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;

  /** 0 if the nonzero blocks are row-oriented, 1 if column-oriented */
  int col_oriented_p;

  enum ownership_mode_t ownership;
  void (*deallocator) (void*);
};


/**
 * @struct csr_matrix_t
 * @author Mark Hoemmen
 * @since 31 May 2005
 *
 * Sparse matrix in Compressed Sparse Row format.
 */
struct
csr_matrix_t
{
  /** Number of rows in the matrix */
  int m;
  
  /** Number of columns in the matrix */
  int n;
  
  /** 
   * Number of stored (nonzero) entries.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored entries, not the actual number of nonzeros. 
   */
  int nnz;

  /** Array of stored (nonzero) entries of the matrix */
  void* values;

  /** Array of column indices of the stored (nonzero) entries of the matrix */
  int* colidx;

  /** Array of indices into the colidx and values arrays, for each column */
  int* rowptr;

  /**
   * Symmetry type of the matrix.
   */
  enum symmetry_type_t symmetry_type;

  /**
   * If the matrix has a kind of symmetry (or skew-symmetry): Where the actual 
   * elements of the matrix are stored: in the lower triangle or the upper 
   * triangle.
   */
  enum symmetric_storage_location_t symmetric_storage_location;

  /** 
   * Indicates the type of the entries in val:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the "val" array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;

  /**
   * The ownership mode: tells whether this library or the user is
   * responsible for deallocating input arrays.
   */
  enum ownership_mode_t ownership;

  /**
   * The deallocation function to be called on the values, colidx and
   * rowptr arrays, if ownership == LIBRARY_DEALLOCATES.
   */
  void (*deallocator) (void*);
};


/**
 * Wrapper struct for a sparse matrix.  
 */
struct
sparse_matrix_t
{
  /** The internal storage format of the sparse matrix. */
  enum sparse_matrix_storage_format_t format;

  /** 
   * Pointer to the internal representation of the sparse matrix.  For example,
   * If format==CSC, this is a "struct csc_matrix_t*"; if format==COO, this is 
   * a "struct coo_matrix_t*". 
   */
  void* repr;
};

/**
 * A coordinate (COO) format sparse matrix.  This uses the ``struct of
 * arrays'' storage paradigm.  An example of the opposite storage
 * paradigm, ``array of structs,'' would be an array of coord_elem_t
 * (which see).
 */
struct
coo_matrix_t
{
  /** Number of rows */
  int m;
  
  /** Number of columns */
  int n;
  
  /** 
   * Number of stored (nonzero) entries.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored entries, not the actual number of nonzeros. 
   */
  int nnz;
 
  /** For entry e, II[e] is its row (we can't call it I because that is reserved as a macro in C99 for sqrt{-1}) */
  int* II;

  /** For entry e, JJ[e] is its column */
  int* JJ;

  /** For entry e, val[e] is its value.  Type of the entries depends 
   * on the value of value_type. */
  void* val;

  /** Base of indexing (either one or zero) */
  enum index_base_t index_base;

  /**
   * Symmetry type of the matrix.
   */
  enum symmetry_type_t symmetry_type;

  /**
   * If the matrix has a kind of symmetry (or skew-symmetry): Where the actual 
   * elements of the matrix are stored: in the lower triangle or the upper 
   * triangle.
   */
  enum symmetric_storage_location_t symmetric_storage_location;

  /** 
   * Indicates the type of the entries in val:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the "val" array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;

  enum ownership_mode_t ownership;
  void (*deallocator) (void*);
};


/**
 * Initializes the given struct A to be the given COO format sparse matrix.
 *
 * @param A [OUT]  Valid pointer to uninitialized struct.
 * @param m [IN]   Number of rows
 * @param n [IN]   Number of columns
 * @param nnz [IN]  Number of nonzero entries
 * @param II [IN]    Array of row indices
 * @param JJ [IN]    Array of column indices
 * @param val [IN]  Array of nonzero values
 * @param index_base [IN]   Base of the indices (ZERO or ONE).
 * @param symmetry_type [IN]  Symmetry type of the matrix
 * @param symmetric_storage_location [IN]  If the matrix has symmetry, where 
 *  the nonzeros are stored.  If the matrix is to be stored in unsymmetric 
 *  format, the value of this enum is irrelevant.
 * @param value_type [IN]  
 */
void
init_coo_matrix (struct coo_matrix_t *A, int m, int n, int nnz, int *II, 
     int *JJ, void *val, enum index_base_t index_base, 
     enum symmetry_type_t symmetry_type, 
     enum symmetric_storage_location_t symmetric_storage_location,
                 enum value_type_t value_type,
     enum ownership_mode_t ownership,
     void (*deallocator) (void*),
     enum copy_mode_t copy_mode);

struct coo_matrix_t*
create_coo_matrix (int m, int n, int nnz, int *II, 
       int *JJ, void *val, enum index_base_t index_base, 
       enum symmetry_type_t symmetry_type, 
       enum symmetric_storage_location_t symmetric_storage_location,
       enum value_type_t value_type,
       enum ownership_mode_t ownership,
       void (*deallocator) (void*),
       enum copy_mode_t copy_mode);


/**
 * Prints the given COO matrix to the given file stream, in MatrixMarket format.
 *
 * @param out [OUT]   Valid file stream
 * @param A [IN]      Sparse matrix in COO format
 *
 * @return Nonzero if error, else zero.
 */
int
print_coo_matrix_in_matrix_market_format (FILE* out, const struct coo_matrix_t* A);

/**
 * Save the given COO matrix to the given file, in MatrixMarket format.
 *
 * @param filename [IN]   Filename to which to save the given sparse matrix
 * @param A [IN]          Sparse matrix in COO format to save
 *
 * @return Nonzero if error, else zero.
 */
int
save_coo_matrix_in_matrix_market_format (const char* const filename, 
           const struct coo_matrix_t* A);
