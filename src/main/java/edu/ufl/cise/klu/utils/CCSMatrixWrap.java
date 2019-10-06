package edu.ufl.cise.klu.utils;

import no.uib.cipr.matrix.sparse.SparseVector;
import org.la4j.matrix.SparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;

import java.lang.reflect.Field;

public class CCSMatrixWrap{

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
