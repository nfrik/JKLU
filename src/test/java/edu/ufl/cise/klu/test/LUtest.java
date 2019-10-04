package edu.ufl.cise.klu.test;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.utils.CCSMatrixWrap;
import no.uib.cipr.matrix.sparse.CompColMatrix;
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

    public static void main(String[] args) throws NoSuchFieldException, IllegalAccessException {

        SparseMatrix sparseMatrix = CCSMatrix.zero(n, n);
//        SparseMatrix sparseMatrix = CCSMatrix.randomSymmetric(n,0.05,new Random());

        b= new double[n];

        Random rand = new Random();

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

        DenseMatrix dm = sparseMatrix.toDenseMatrix();
//        System.out.println(dm);
        Ap= CCSMatrixWrap.getColumnPointers(sparseMatrix);
        Ai= CCSMatrixWrap.getRowIndices(sparseMatrix);
        Ax= CCSMatrixWrap.getValues(sparseMatrix);

//        sparseMatrix.set(0,0,23f);



        KLU_symbolic Symbolic;
        KLU_numeric Numeric;
        KLU_common Common = new KLU_common();

        //Dklu_version.NPRINT = false ;
        //Dklu_internal.NDEBUG = false ;

        klu_defaults (Common);

        long startTime = System.nanoTime();
        Symbolic = klu_analyze (n, Ap, Ai, Common);
        System.out.println("Execution time msec: "+ (System.nanoTime()-startTime)/1000000);
        Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
        System.out.println("Execution time msec: "+ (System.nanoTime()-startTime)/1000000);
        klu_solve (Symbolic, Numeric, n, 1, b, 0, Common);

        System.out.println("Execution time msec: "+ (System.nanoTime()-startTime)/1000000);
        System.out.println("Peak memory:"+ Common.mempeak+" status "+Common.status);
//        for (i = 0 ; i < n ; i++) {
//            System.out.printf("x [%d] = %g\n", i, b [i]) ;
////            assertEquals(i + 1.0, b [i], DELTA) ;
//        }
    }
}
