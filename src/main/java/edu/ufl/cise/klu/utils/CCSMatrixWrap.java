package edu.ufl.cise.klu.utils;

import no.uib.cipr.matrix.sparse.SparseVector;
import org.la4j.matrix.SparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;

import java.lang.reflect.Field;

public class CCSMatrixWrap{

//        public CCSMatrixWrap(int r,int c){
////            this(r,c);
//            Field rowIndicesField = null;
//            Field columnPointersField = null;
//            Field valuesField = null;
//            try {
//                rowIndicesField = CCSMatrix.class.getDeclaredField("rowIndices");
//                columnPointersField = CCSMatrix.class.getDeclaredField("columnPointers");
//                valuesField = CCSMatrix.class.getDeclaredField("values");
//            } catch (NoSuchFieldException e) {
//                e.printStackTrace();
//            }
//
//            rowIndicesField.setAccessible(true);
//            columnPointersField.setAccessible(true);
//            valuesField.setAccessible(true);
//        }


//    public SparseMatrix sparseMatrix=null;
////    CCSMatrix.zero(n, n);
//    public CCSMatrixWrap(int n){
//        this.sparseMatrix = CCSMatrix.zero(n,n);
//
//        Field rowIndicesField = null;
//        Field columnPointersField = null;
//        Field valuesField = null;
//        try {
//            rowIndicesField = CCSMatrix.class.getDeclaredField("rowIndices");
//            columnPointersField = CCSMatrix.class.getDeclaredField("columnPointers");
//            valuesField = CCSMatrix.class.getDeclaredField("values");
//        } catch (NoSuchFieldException e) {
//            e.printStackTrace();
//        }
//
//        rowIndicesField.setAccessible(true);
//        columnPointersField.setAccessible(true);
//        valuesField.setAccessible(true);
//    }
//
//
//    public SparseMatrix getSparseMatrix() {
//        return sparseMatrix;
//
//    }
//    public void set(int i, int i1, double v){
//        this.sparseMatrix.set(i,i1,v);
//    }
//
//    public double get(int i, int i1){
//        return this.sparseMatrix.get(i,i1);
//    }
//
////    public void CCSMatrix.randomSymmetric()

    public static int[] getRowIndices(SparseMatrix sparseMatrix){

        Field rowIndicesField = null;
        try {
            rowIndicesField = CCSMatrix.class.getDeclaredField("rowIndices");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert rowIndicesField != null;
        rowIndicesField.setAccessible(true);
        try {
            return (int[]) rowIndicesField.get(sparseMatrix);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
        return null;
    }

    public static void setRowIndices(SparseMatrix sparseMatrix, int[] values){

        Field rowIndicesField = null;
        try {
            rowIndicesField = CCSMatrix.class.getDeclaredField("rowIndices");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert rowIndicesField != null;
        rowIndicesField.setAccessible(true);
        try {
//            return (int[]) rowIndicesField.get(sparseMatrix);
            rowIndicesField.set(sparseMatrix,values);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
//        return null;
    }

    public static int[] getColumnPointers(SparseMatrix sparseMatrix){

        Field columnPointersField = null;

        try {
            columnPointersField = CCSMatrix.class.getDeclaredField("columnPointers");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert columnPointersField != null;
        columnPointersField.setAccessible(true);
        try {
            return (int[]) columnPointersField.get(sparseMatrix);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
        return null;
    }

    public static void setColumnPointers(SparseMatrix sparseMatrix,int[] values){

        Field columnPointersField = null;

        try {
            columnPointersField = CCSMatrix.class.getDeclaredField("columnPointers");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert columnPointersField != null;
        columnPointersField.setAccessible(true);
        try {
//            return (int[]) columnPointersField.get(sparseMatrix);
            columnPointersField.set(sparseMatrix,values);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
//        return null;
    }

    public static double[] getValues(SparseMatrix sparseMatrix){

        Field valuesField = null;
        try {
            valuesField = CCSMatrix.class.getDeclaredField("values");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert valuesField != null;
        valuesField.setAccessible(true);
        try {
            return (double[]) valuesField.get(sparseMatrix);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
        return null;
    }

    public static void setValues(SparseMatrix sparseMatrix,double[] values){

        Field valuesField = null;
        try {
            valuesField = CCSMatrix.class.getDeclaredField("values");
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
        }

        assert valuesField != null;
        valuesField.setAccessible(true);
        try {
//            return (double[]) valuesField.get(sparseMatrix);
              valuesField.set(sparseMatrix,values);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        ;
//        return null;
    }

}
