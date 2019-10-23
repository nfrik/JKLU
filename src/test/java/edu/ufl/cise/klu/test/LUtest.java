package edu.ufl.cise.klu.test;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.utils.CCSMatrixWrap;
import edu.ufl.cise.klu.wrapper.*;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompColMatrix;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import no.uib.cipr.matrix.sparse.SparseVector;
import org.la4j.Vector;
import org.la4j.matrix.DenseMatrix;
import org.la4j.matrix.SparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;
import org.la4j.matrix.sparse.CRSMatrix;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_factor.klu_factor;
import static edu.ufl.cise.klu.tdouble.Dklu_solve.klu_solve;
import static edu.ufl.cise.klu.tdouble.Dklu_version.KLU_SINGULAR;


import java.lang.reflect.Field;
import java.util.Random;

public class LUtest {

    private static final double DELTA = 1e-09 ;

    private static int n = 4000 ;
    private static int [ ] Ap = {0, 2, 5, 9, 10, 12} ;
    private static int [ ] Ai = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
    private static double [ ] Ax = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
    private static double [ ] b = {8., 45., -3., 3., 19.} ;

    /**
     * a simple KLU demo; solution is x = (1,2,3,4,5)
     */

    public static void mtj_ccs_klu_test(){
        long startTime = System.nanoTime();

        FlexCompColMatrix sparseMatrix = new FlexCompColMatrix(n,n);
        CompColMatrix tst = new CompColMatrix(sparseMatrix);

        System.out.println("Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);

        startTime = System.nanoTime();

        b= new double[n];

        Random rand = new Random();
        rand.setSeed(1212);

        for(int i=0;i<sparseMatrix.numColumns();i++){
            for(int m=0;m<rand.nextInt(50);m++){
                int k = rand.nextInt(sparseMatrix.numRows());
                double val = rand.nextDouble();
                sparseMatrix.set(i,k,val);
                sparseMatrix.set(k,i,val);
                sparseMatrix.set(i,i,rand.nextDouble());
            }

            b[i]=rand.nextDouble();
        }
        System.out.println("Matrix initialized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        CompColMatrix ccsm = new CompColMatrix(sparseMatrix);

//        new CompColMatrix()

        Ap= ccsm.getColumnPointers();
        Ai= ccsm.getRowIndices();
        Ax= ccsm.getData();

        System.out.println("Retrieved CCS data time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        KLU_symbolic Symbolic;
        KLU_numeric Numeric;
        KLU_common Common = new KLU_common();

        //Dklu_version.NPRINT = false ;
        //Dklu_internal.NDEBUG = false ;

        klu_defaults (Common);


        Symbolic = klu_analyze (n, Ap, Ai, Common);
        System.out.println("Klu analyzed time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
        System.out.println("Klu factorized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        klu_solve (Symbolic, Numeric, n, 1, b, 0, Common);

        System.out.println("Klu solved time msec: "+ (System.nanoTime()-startTime)/1000000);
        System.out.println("Peak memory:"+ Common.mempeak+" status "+Common.status);

    }

    public static void la4j_ccs_klu_test(){
        long startTime = System.nanoTime();
        SparseMatrix sparseMatrix = CCSMatrix.zero(n, n);
        System.out.println("Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);

        startTime = System.nanoTime();

        b= new double[n];

        Random rand = new Random();
        rand.setSeed(1212);

        for(int i=0;i<sparseMatrix.rows();i++){
            for(int m=0;m<rand.nextInt(50);m++){
                int k = rand.nextInt(sparseMatrix.rows());
                double val = rand.nextDouble();
                sparseMatrix.set(i,k,val);
                sparseMatrix.set(k,i,val);
                sparseMatrix.set(i,i,rand.nextDouble());
            }

            b[i]=rand.nextDouble();
        }
        System.out.println("Matrix initialized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        Ap= CCSMatrixWrap.getColumnPointers(sparseMatrix);
        Ai= CCSMatrixWrap.getRowIndices(sparseMatrix);
        Ax= CCSMatrixWrap.getValues(sparseMatrix);


        System.out.println("Retrieved CCS data time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        KLU_symbolic Symbolic;
        KLU_numeric Numeric;
        KLU_common Common = new KLU_common();

        //Dklu_version.NPRINT = false ;
        //Dklu_internal.NDEBUG = false ;

        klu_defaults (Common);


        Symbolic = klu_analyze (n, Ap, Ai, Common);
        System.out.println("Klu analyzed time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
        System.out.println("Klu factorized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        klu_solve (Symbolic, Numeric, n, 1, b, 0, Common);

        System.out.println("Klu solved time msec: "+ (System.nanoTime()-startTime)/1000000);
        System.out.println("Peak memory:"+ Common.mempeak+" status "+Common.status);
//        for (i = 0 ; i < n ; i++) {
//            System.out.printf("x [%d] = %g\n", i, b [i]) ;
////            assertEquals(i + 1.0, b [i], DELTA) ;
//        }
    }

    private static void superluwr_ccs_klu_test()
    {
        System.out.println(""); System.out.println("starting test with SuperLU native");
        long startTime = System.nanoTime();

        SuperMatrix sparseMatrix = new SuperMatrix();
        Stype_t stype = new Stype_t();
        Dtype_t dtype = new Dtype_t();
        Mtype_t mtype = new Mtype_t();
        double[] nzval = null;
        int[] rowind = null;
        int[] colptr = null;
        SuperLUWrapper.dCreate_CompCol_Matrix(sparseMatrix, n, n, 0, nzval, rowind, colptr, stype, dtype, mtype);

        //SparseMatrix sparseMatrix = CCSMatrix.zero(n, n);
        System.out.println("Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);

        startTime = System.nanoTime();

        b= new double[n];

        Random rand = new Random();
        rand.setSeed(1212);

        for(int i=0;i<sparseMatrix.rows();i++){
            for(int m=0;m<rand.nextInt(50);m++){
                int k = rand.nextInt(sparseMatrix.rows());
                double val = rand.nextDouble();
                sparseMatrix.set(i,k,val);
                sparseMatrix.set(k,i,val);
                sparseMatrix.set(i,i,rand.nextDouble());
            }

            b[i]=rand.nextDouble();
        }
        System.out.println("Matrix initialized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        //Ap= CCSMatrixWrap.getColumnPointers(sparseMatrix);
        //Ai= CCSMatrixWrap.getRowIndices(sparseMatrix);
        //Ax= CCSMatrixWrap.getValues(sparseMatrix);


        System.out.println("Retrieved CCS data time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        KLU_symbolic Symbolic;
        KLU_numeric Numeric;
        KLU_common Common = new KLU_common();

        //Dklu_version.NPRINT = false ;
        //Dklu_internal.NDEBUG = false ;

        klu_defaults (Common);


        Symbolic = klu_analyze (n, Ap, Ai, Common);
        System.out.println("Klu analyzed time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
        System.out.println("Klu factorized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        klu_solve (Symbolic, Numeric, n, 1, b, 0, Common);

        System.out.println("Klu solved time msec: "+ (System.nanoTime()-startTime)/1000000);
        System.out.println("Peak memory:"+ Common.mempeak+" status "+Common.status);
//        for (i = 0 ; i < n ; i++) {
//            System.out.printf("x [%d] = %g\n", i, b [i]) ;
////            assertEquals(i + 1.0, b [i], DELTA) ;
//        }
    }

    private static void short_superLU_wrapped_ccs_test()
    {
        //partial copy ov la4j_ccs_klu

        long startTime = System.nanoTime();
        SparseMatrix sparseMatrix = CCSMatrix.zero(n, n);
        System.out.println("[LA4J]  Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        b= new double[n];
        Random rand = new Random();
        rand.setSeed(1212);

        for(int i=0;i<sparseMatrix.rows();i++){
            for(int m=0;m<rand.nextInt(50);m++){
                int k = rand.nextInt(sparseMatrix.rows());
                double val = rand.nextDouble();
                sparseMatrix.set(i,k,val);
                sparseMatrix.set(k,i,val);
                sparseMatrix.set(i,i,rand.nextDouble());
            }

            b[i]=rand.nextDouble();
        }
        System.out.println("Matrix initialized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        Ap= CCSMatrixWrap.getColumnPointers(sparseMatrix);
        Ai= CCSMatrixWrap.getRowIndices(sparseMatrix);
        Ax= CCSMatrixWrap.getValues(sparseMatrix);

        System.out.println("Retrieved CCS data time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        SuperLUWrapper.ccs_components_b_dgssv(n,n, Ap, Ai, Ax, b);

        System.out.println("Klu solved time msec: "+ (System.nanoTime()-startTime)/1000000);
        //System.out.println("Peak memory:"+ Common.mempeak+" status "+Common.status);
//        for (i = 0 ; i < n ; i++) {
//            System.out.printf("x [%d] = %g\n", i, b [i]) ;
////            assertEquals(i + 1.0, b [i], DELTA) ;
//        }
    }

    public static void main(String[] args) throws NoSuchFieldException, IllegalAccessException {
        mtj_ccs_klu_test();
        la4j_ccs_klu_test();
        //superluwr_ccs_klu_test();
        short_superLU_wrapped_ccs_test();
    }
}
