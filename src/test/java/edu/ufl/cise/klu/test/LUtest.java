package edu.ufl.cise.klu.test;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.utils.CCSMatrixWrap;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.sparse.CompColMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import org.ejml.data.DMatrixSparseCSC;
import org.la4j.matrix.SparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_factor.klu_factor;
import static edu.ufl.cise.klu.tdouble.Dklu_solve.klu_solve;

import java.util.Iterator;
import java.util.Random;

public class LUtest {

    private static final double DELTA = 1e-09 ;

    private static int n = 10000 ;
    private static int diagm = 1000 ;
    private static int [ ] Ap = {0, 2, 5, 9, 10, 12} ;
    private static int [ ] Ai = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
    private static double [ ] Ax = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
    private static double [ ] b = {8., 45., -3., 3., 19.} ;

    /**
     * a simple KLU demo; solution is x = (1,2,3,4,5)
     */

    public static class SparseUtils{

        public double[] values=null;
        public int[] rows=null;
        public int[] cols=null;
        public CompColMatrix ccm = null;
        public FlexCompColMatrix fccm = null;

        public void get_triplet(CompColMatrix ccm){
            this.ccm=ccm;
            int n = ccm.getRowIndices().length;
            String format = "%d %d %11.4E\n";
            int numRows = ccm.numRows();
            int numCols = ccm.numColumns();
            int nz_length = ccm.getData().length;

            rows = new int[n];
            cols = new int[n];
            values = new double[n];

            int[] col_idx=ccm.getColumnPointers();
            int[] nz_rows=ccm.getRowIndices();
            double[] nz_values=ccm.getData();

            System.out.println(" , rows = " + numRows + " , cols = " + numCols + " , nz_length = " + nz_length);

            int globidx=0;
            for(int col = 0; col < numCols; ++col) {
                int idx0 = col_idx[col];
                int idx1 = col_idx[col + 1];

                for(int i = idx0; i < idx1; ++i) {
                    int row = nz_rows[i];
                    double value = nz_values[i];
                    values[globidx]=value;
                    rows[globidx]=row;
                    cols[globidx]=col;
                    globidx++;
//                    System.out.printf(format, row, col, value);
                }
            }
        }
        public void get_triplet(FlexCompColMatrix fccm){
            this.fccm=fccm;
            this.ccm=new CompColMatrix(fccm);
            get_triplet(fccm);
        }
    }



    public static void mtj_ccs_klu_test(){
        long startTime = System.nanoTime();

        FlexCompColMatrix sparseMatrix = new FlexCompColMatrix(n,n);
//        LinkedSparseMatrix sparseMatrix = new LinkedSparseMatrix(n,n);
//        CompColMatrix tst = new CompColMatrix(sparseMatrix);

        System.out.println("Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);

        startTime = System.nanoTime();

        b= new double[n];

        Random rand = new Random();
        rand.setSeed(1212);

        for(int i=0;i<n;i++){
            int randm=rand.nextInt(diagm);
            for(int m=0;m<randm;m++){
                int k = rand.nextInt(n);
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

        System.out.println("Converted to CCS time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

        SparseUtils su = new SparseUtils();
        su.get_triplet(ccsm);
        System.out.println(su.values.length);
        System.out.println("Extracted triplet time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
//        new CompColMatrix()

//        for(Iterator elem = ccsm.iterator(); elem.hasNext();){
//            Object nxt = (Object) elem.next();
//            System.out.println(nxt);
//        }
//        sparseMatrix.forEach(double a :);

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
            int randm=rand.nextInt(diagm);
            for(int m=0;m<randm;m++){
                int k = rand.nextInt(sparseMatrix.rows());
                double val = rand.nextDouble();
                sparseMatrix.set(i,k,val);
                sparseMatrix.set(k,i,val);
                sparseMatrix.set(i,i,rand.nextDouble());
            }

            b[i]=rand.nextDouble();
        }

//        sparseMatrix.iterator();
        for(Iterator elem = sparseMatrix.iterator(); elem.hasNext();){
            MatrixEntry nxt = (MatrixEntry) elem.next();
            System.out.println(nxt.get());
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


    public static void ejml_ccs_klu_test(){
        long startTime = System.nanoTime();

//        SparseMatrix sparseMatrix = CCSMatrix.zero(n, n);
//        DMatrixSparseTriplet sparseTriplet = new DMatrixSparseTriplet(n,n,10);
        DMatrixSparseCSC sparseTriplet=new DMatrixSparseCSC(n,n,40000);
        System.out.println("Matrix instantiated time msec: "+ (System.nanoTime()-startTime)/1000000);


        startTime = System.nanoTime();

        b= new double[n];

        Random rand = new Random();
        rand.setSeed(1212);

        int count=0;
        for(int i=0;i<n;i++){
            int randm=rand.nextInt(diagm);
            for(int m=0;m<randm;m++){
                int k = rand.nextInt(n);
                double val = rand.nextDouble();
                sparseTriplet.set(i,k,val);
                sparseTriplet.set(k,i,val);
                sparseTriplet.set(i,i,rand.nextDouble());
                count++;
            }

            b[i]=rand.nextDouble();
        }

        System.out.println("Matrix initialized time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();

//        DMatrixSparseCSC sparseCSC = ConvertDMatrixStruct.convert(sparseTriplet,(DMatrixSparseCSC)null);
        DMatrixSparseCSC sparseCSC=sparseTriplet;
//        sparseCSC.printNonZero();
        System.out.println("Matrix converted to CSC time msec: "+ (System.nanoTime()-startTime)/1000000);
        startTime = System.nanoTime();
        Ap= sparseCSC.col_idx;
        Ai= sparseCSC.nz_rows;
        Ax= sparseCSC.nz_values;


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

    public static void main(String[] args) throws NoSuchFieldException, IllegalAccessException {
        mtj_ccs_klu_test();
//        la4j_ccs_klu_test();
//        ejml_ccs_klu_test();

    }
}
