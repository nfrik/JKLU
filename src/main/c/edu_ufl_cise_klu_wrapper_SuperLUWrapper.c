#include "edu_ufl_cise_klu_wrapper_SuperLUWrapper.h"
#include <stdio.h>
#include "slu_ddefs.h"

#include <time.h>

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_dCreate_1CompCol_1Matrix
  (JNIEnv *env, jclass cls, jobject superMat, jint rows, jint cols, jint bum,
	jdoubleArray arra, jintArray size, jintArray iArra, jobject ob, jobject ob2, jobject ob3) {
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_dGenXtrue
  (JNIEnv *env, jclass cls, jint sz, jint sz2, jdoubleArray arra, jint sz3)
{
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_dFillRHS
  (JNIEnv *env, jclass cls, jobject obj, jint sz, jdoubleArray dArra, jint sz2,
	jobject obj1, jobject obj2)
{
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_get_1perm_1c
  (JNIEnv *env, jclass cls, jint sz, jobject obj1, jintArray iArra)
{
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_dgssv
  (JNIEnv *env, jclass cls, jobjectArray objArra, jobject obj1, jintArray iArra,
	jintArray iArra2, jobject obj2, jobject obj3, jobject obj4,
	jobjectArray objArra2, jintArray iArra3)
{
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_dinf_1norm_1error
  (JNIEnv *env, jclass cls, jint sz, jobject obj1, jdoubleArray dArra)
{
	printf("Called %s\n", __func__ );
	return;
}

JNIEXPORT void JNICALL Java_edu_ufl_cise_klu_wrapper_SuperLUWrapper_ccs_1components_1b_1dgssv
  (JNIEnv *env, jclass cld, jint rows, jint cols, jintArray colsPtrs, jintArray rowPointed,
	jdoubleArray values, jdoubleArray b_vector)
{
	SuperMatrix A;
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    superlu_options_t options;


int      *perm_c; /* column permutation vector */
int      *perm_r; /* row permutations from partial pivoting */
int      nrhs, ldx, info, m, n, nnz;
double   *xact, *rhs;
mem_usage_t   mem_usage;
SuperLUStat_t stat;

struct timespec all_start, solve_start, solve_stop, all_stop;

    jsize vals_size, cols_size, row_size, b_size;

    fprintf(stderr, "Called %s\n", __func__ );
    vals_size = (*env)->GetArrayLength(env, values);
    cols_size = (*env)->GetArrayLength(env, colsPtrs);
    row_size = (*env)->GetArrayLength(env, rowPointed);
    b_size = (*env)->GetArrayLength(env, b_vector);

    nnz = vals_size;
    m = rows;
    n = cols;
    nrhs = 1;

    jdouble *values_double = (*env)->GetDoubleArrayElements(env, values, 0);
    jint *rows_int = (*env)->GetIntArrayElements(env, rowPointed, 0);
    jint *cols_int = (*env)->GetIntArrayElements(env, colsPtrs, 0);
    jdouble *b_double = (*env)->GetDoubleArrayElements(env, b_vector, 0);

clock_gettime( CLOCK_MONOTONIC, &all_start);

        dCreate_CompCol_Matrix(&A, rows, cols, vals_size, values_double, rows_int, cols_int, SLU_NC, SLU_D, SLU_GE);
        dCreate_Dense_Matrix(&B, rows, 1,b_double, b_size, SLU_DN, SLU_D, SLU_GE);
        {
            double *sol = (double*) ((DNformat*) B.Store)->nzval;
            for(int i =0; i < n; ++i) {
                printf("[%d] < %f  \t > %f\n", i,b_double[i], sol[i]  );
            }
        }

        set_default_options(&options);
        xact = doubleMalloc(n * nrhs);
      ldx = n;
      dGenXtrue(n, nrhs, xact, ldx);
      dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

      if ( !(perm_c = intMalloc(n)) ) return;
      if ( !(perm_r = intMalloc(m)) ) return;

      /* Initialize the statistics variables. */
      StatInit(&stat);

clock_gettime( CLOCK_MONOTONIC, &solve_start);
      dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

clock_gettime( CLOCK_MONOTONIC, &solve_stop);
printf("info = %d\n", info);
      if ( info == 0 ) {
        // copy
        double *sol = (double*) ((DNformat*) B.Store)->nzval;
        for(int i=0; i < 10; ++i) {
            printf("[%d] -> %f____\n", i, sol[i]);
        }
        for(int i=b_size-10; i < b_size; ++i) {
                    printf("[%d] -> %f____\n", i, sol[i]);
                }
        memcpy(b_double,sol, b_size*sizeof(double));
      }


    (*env)->ReleaseDoubleArrayElements(env, values, values_double, 0);
    (*env)->ReleaseDoubleArrayElements(env, b_vector, b_double, 0);
    (*env)->ReleaseIntArrayElements(env, rowPointed, rows_int, 0);
    (*env)->ReleaseIntArrayElements(env, colsPtrs, cols_int, 0);

clock_gettime( CLOCK_MONOTONIC, &all_stop);

printf("All time = %f\n only solution = %f\n",
	all_stop.tv_sec - all_start.tv_sec + 1e-9 * (all_stop.tv_nsec - all_start.tv_nsec),
	solve_stop.tv_sec - solve_start.tv_sec + 1e-9 * (solve_stop.tv_nsec - solve_start.tv_nsec)
	);
	return;
}
