/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "dsp_defs.h"
#include "util.h"
#include "Cnames.h"

#ifdef _CRAY
#define c_fortran_dgssv_  C_FORTRAN_DGSSV
#endif

#define HANDLE_SIZE  8

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_r;
    int *perm_c;
} factors_t;

int
c_fortran_dgssv_(int *iopt, int *n, int *nnz, int *nrhs, double *values,
		 int *rowind, int *colptr, double *b, int *ldb,
		 int factors[HANDLE_SIZE], /* a handle containing the pointer
					      to the factored matrices */
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * factors (input/output) int array of size 8
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    char     refact[1], trans[1];
    double   diag_pivot_thresh;
    mem_usage_t   mem_usage;
    factors_t *LUfactors;

    *trans = 'N';

    if ( *iopt == 1 ) { /* LU decomposition */
	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i) --colptr[i];

	dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
			       NC, _D, GE);
	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");

	*refact = 'N';
	diag_pivot_thresh = 1.0;
	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	StatInit(panel_size, relax);

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	permc_spec = 3;
	get_perm_c(permc_spec, &A, perm_c);
	
	sp_preorder(refact, &A, perm_c, etree, &AC);

	dgstrf(refact, &AC, diag_pivot_thresh, 0.0, relax, panel_size, 
	       etree, NULL, 0, perm_r, perm_c, L, U, info);

	if ( *info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
	    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	    dQuerySpace(L, U, panel_size, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		   mem_usage.expansions);
	} else {
	    printf("dgstrf() error returns INFO= %d\n", *info);
	    if ( *info <= *n ) { /* factorization completes */
		dQuerySpace(L, U, panel_size, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		       mem_usage.expansions);
	    }
	}
	
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_r = perm_r;
	LUfactors->perm_c = perm_c;
	factors[0] = (int) LUfactors;

	/* Free un-wanted storage */
	SUPERLU_FREE(etree);
	Destroy_SuperMatrix_Store(&A);
	Destroy_CompCol_Permuted(&AC);
    } else if ( *iopt == 2 ) { /* Triangular solve */
	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) factors[0];
	L = LUfactors->L;
	U = LUfactors->U;
	perm_r = LUfactors->perm_r;
	perm_c = LUfactors->perm_c;

	dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, DN, _D, GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs (trans, L, U, perm_r, perm_c, &B, info);

	Destroy_SuperMatrix_Store(&B);
    } else if ( *iopt == 3 ) { /* Free storage */
	/* Free the LU factors in the factors handle */
	LUfactors = (factors_t*) factors[0];
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
	Destroy_SuperNode_Matrix(LUfactors->L);
	Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
	StatFree();
    } else {
	fprintf(stderr, "Invalid iopt=%d passed to c_bridge_pdgssv()\n");
	exit(-1);
    }
}


