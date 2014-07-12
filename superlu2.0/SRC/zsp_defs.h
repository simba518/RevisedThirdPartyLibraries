

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#ifndef __SUPERLU_zSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_zSP_DEFS

/*
 * File name:		zsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#include <string.h>
#endif

/* Define my integer type int_t */
typedef int int_t; /* default */

#include "Cnames.h"
#include "supermatrix.h"
#include "dcomplex.h"


/* No of marker arrays used in the symbolic factorization,
   each of size n */
#define NO_MARKER     3
#define NUM_TEMPV(m,w,t,b)  ( SUPERLU_MAX(m, (t + b)*w) )

typedef enum {LUSUP, UCOL, LSUB, USUB} MemType;
typedef enum {HEAD, TAIL}              stack_end_t;
typedef enum {SYSTEM, USER}            LU_space_t;

/*
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 */
typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    doublecomplex  *lusup;   /* L supernodes */
    int     *xlusup;
    doublecomplex  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} GlobalLU_t;

typedef struct {
    int panel_size;
    int relax;
    double diag_pivot_thresh;
    double drop_tol;
} factor_param_t;

typedef struct {
    float for_lu;
    float total_needed;
    int   expansions;
} mem_usage_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Driver routines */
extern void
zgssv(SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
	SuperMatrix *, int *);
extern void
zgssvx(char *, char *, char *, SuperMatrix *, factor_param_t *,
       int *, int *, int *, char *, double *, double *,
       SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, 
       SuperMatrix *, double *, double *, double *,
       double *, mem_usage_t *, int *);

/* Supernodal LU factor related */
extern void
zCreate_CompCol_Matrix(SuperMatrix *, int, int, int, doublecomplex *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCreate_Dense_Matrix(SuperMatrix *, int, int, doublecomplex *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, doublecomplex *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_Dense_Matrix(int, int, doublecomplex *, int, doublecomplex *, int);

extern void    Destroy_SuperMatrix_Store(SuperMatrix *);
extern void    Destroy_CompCol_Matrix(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix(SuperMatrix *);
extern void    Destroy_CompCol_Permuted(SuperMatrix *);
extern void    Destroy_Dense_Matrix(SuperMatrix *);
extern void    get_perm_c(int, SuperMatrix *, int *);
extern void    sp_preorder (char*, SuperMatrix*, int*, int*, SuperMatrix*);
extern void    countnz (const int, int *, int *, int *, GlobalLU_t *);
extern void    fixupL (const int, const int *, GlobalLU_t *);

extern void    zallocateA (int, int, doublecomplex **, int **, int **);
extern void    zgstrf (char*, SuperMatrix*, double, double, int, int, int*,
			void *, int, int *, int *, 
                        SuperMatrix *, SuperMatrix *, int *);
extern int     zsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
extern int     zsnode_bmod (const int, const int, const int, doublecomplex *,
                              doublecomplex *, GlobalLU_t *);
extern void    zpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, doublecomplex *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
extern void    zpanel_bmod (const int, const int, const int, const int,
                           doublecomplex *, doublecomplex *, int *, int *,
			   GlobalLU_t *);
extern int     zcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
extern int     zcolumn_bmod (const int, const int, doublecomplex *,
			   doublecomplex *, int *, int *, int, GlobalLU_t *);
extern int     zcopy_to_ucol (int, int, int *, int *, int *,
                              doublecomplex *, GlobalLU_t *);         
extern int     zpivotL (const int, const double, int *, int *, 
                              int *, int *, int *, GlobalLU_t *);
extern void    zpruneL (const int, const int *, const int, const int,
			     const int *, const int *, int *, GlobalLU_t *);
extern void    zreadmt (int *, int *, int *, doublecomplex **, int **, int **);
extern void    zGenXtrue (int, int, doublecomplex *, int);
extern void    zFillRHS (char *, int, doublecomplex *, int, SuperMatrix *,
			SuperMatrix *);
extern void    zgstrs (char *, SuperMatrix *, SuperMatrix *, int *, int *,
			SuperMatrix *, int *);


/* Driver related */

extern void    zgsequ (SuperMatrix *, double *, double *, double *,
			     double *, double *, int *);
extern void    zlaqgs (SuperMatrix *, double *, double *, double,
                             double, double, char *);
extern void    zgscon (char *, SuperMatrix *, SuperMatrix *,
			double, double *, int *);
extern double  zPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    zgsrfs (char *, SuperMatrix *, SuperMatrix *, 
			SuperMatrix *, int *, int *, char *, double *,
			double *, SuperMatrix *, SuperMatrix *, 
			double *, double *, int *);

extern int     sp_ztrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, doublecomplex *, int *);
extern int     sp_zgemv (char *, doublecomplex, SuperMatrix *, doublecomplex *,
			int, doublecomplex, doublecomplex *, int);

extern int     sp_zgemm (char *, char *, int, int, int, doublecomplex,
			SuperMatrix *, doublecomplex *, int, doublecomplex, 
			doublecomplex *, int);

/* Memory-related */
extern int     zLUMemInit (char *, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     GlobalLU_t *, int **, doublecomplex **);
extern void    zSetRWork (int, int, doublecomplex *, doublecomplex **, doublecomplex **);
extern void    zLUWorkFree (int *, doublecomplex *, GlobalLU_t *);
extern int     zLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

extern doublecomplex  *doublecomplexMalloc(int);
extern doublecomplex  *doublecomplexCalloc(int);
extern double  *doubleMalloc(int);
extern double  *doubleCalloc(int);
extern int     zmemory_usage(const int, const int, const int, const int);
extern int     zQuerySpace (SuperMatrix *, SuperMatrix *, int,
				mem_usage_t *);

/* Auxiliary routines */
extern void    zreadhb(int *, int *, int *, doublecomplex **, int **, int **);
extern void    zCompRow_to_CompCol(int, int, int, doublecomplex*, int*, int*,
		                   doublecomplex **, int **, int **);
extern void    zfill (doublecomplex *, int, doublecomplex);
extern void    zinf_norm_error (int, SuperMatrix *, doublecomplex *);
extern void    PrintPerf (SuperMatrix *, SuperMatrix *, mem_usage_t *,
			 doublecomplex, doublecomplex, doublecomplex *, doublecomplex *, char *);

/* Routines for debugging */
extern void    zPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    zPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    zPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    print_lu_col(char *, int, int, int *, GlobalLU_t *);
extern void    check_tempv(int, doublecomplex *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_zSP_DEFS */

