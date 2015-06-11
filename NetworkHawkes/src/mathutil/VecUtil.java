package mathutil;

import java.util.ArrayList;
import java.util.Arrays;

public class VecUtil {
	//  function to normalize vectors
	public static double[] norm(double[] x){
		double[] xnorm = new double[x.length];
		double sumx = 0;
		for(int i = 0; i < x.length; i++) sumx += x[i];
		if(sumx == 0){
			for(int i = 0 ; i < xnorm.length; i++) xnorm[i] = 1;
		}else{
			for(int i = 0; i < xnorm.length; i++) xnorm[i] = x[i] / sumx;
		}
		return(xnorm);
	}

	// function to find max of selected elements in an array.
	public static double array_max(double[] array,  ArrayList<Integer> location){
		double max = Double.MIN_VALUE;
		for(int i : location){
			if(array[i] > max) max = array[i];
		}
		return(max);
	}	

	// function to grab certain column of 2d array
	public static int[] grab2(int[][] matrix, int col){
		int[] out = new int[matrix.length];
		for(int i = 0; i < out.length; i++) out[i] = matrix[i][col];
		return(out);
	}
	//	// function to grab certain column of 2d array
	//	public static int[] grab2(int[][] matrix, int col){
	//				int[] out = new int[matrix.length];
	//				for(int i = 0; i < out.length; i++) out[i] = matrix[i][col];
	//				return(out);
	//	}
	//	
	// function to fill in constant array
	public static double[] fill(double[] vec, double value){
		for(int i = 0; i < vec.length; i++){ vec[i] = value; }
		return(vec);
	}

	// function to sum Adjencency matrix
	public static int sum2d(int[][] matrix){
		int sum = 0;
		for(int i = 0; i < matrix.length; i++){
			for(int j = 0; j < matrix[0].length; j++){
				sum += matrix[i][j];
			}
		}
		return(sum);
	}

	// function to sum Weight matrix
	public static double sum2d(double[][] matrix){
		double sum = 0;
		for(int i = 0; i < matrix.length; i++){
			for(int j = 0; j < matrix[0].length; j++){
				sum += matrix[i][j];
			}
		}
		return(sum);
	}
	// function to select part of the array
	public static double[] selectArray(double[] array, ArrayList<Integer> index){
		double[] out = new double[index.size()];
		for(int i =0; i < index.size(); i++){
			out[i] = array[index.get(i)];
		}
		return(out);
	}
	public static int[] selectArray(int[] array, ArrayList<Integer> index){
		int[] out = new int[index.size()];
		for(int i =0; i < index.size(); i++){
			out[i] = array[index.get(i)];
		}
		return(out);
	}
	// function to calculate sum
	public static double getSum(double[] array){
		double m = 0;
		for(int i = 0; i < array.length; i++){m += array[i];}
		return(m);
	}
	// function to calculate sum
	public static int getSum(int[] array){
		int m = 0;
		for(int i = 0; i < array.length; i++){m += array[i];}
		return(m);
	}	
	// function to calculate mean
	public static double getMean(int[] array){
		double  m = 0;
		for(int i = 0; i < array.length; i++){m += array[i];}
		m = m / (array.length + 0.0);
		return(m);
	}
	// function to calculate mean
	public static double getMean(double[] array){
		double  m = 0;
		for(int i = 0; i < array.length; i++){m += array[i];}
		m = m / (array.length + 0.0);
		return(m);
	}
	// function to calculate mean
	public static double getVar(double[] array){
		double  m = 0;
		double m2 = 0;
		for(int i = 0; i < array.length; i++){m += array[i]; m2 += array[i] * array[i];}
		m = m / (array.length + 0.0);
		m2 = m2 / (array.length + 0.0);
		return(m2 - m*m);
	}
	
	public static String write3d(double[][][] array){
		StringBuilder out = new StringBuilder();
		for(int i = 0; i < array.length; i++){
			for(int j = 0; j < array[0].length; j++){
				out.append(Arrays.toString(array[i][j]).replace("[", "").replace("]", "\n"));				
			}
		}
		return(out.toString());
	}
	public static String write3d(int[][][] array){
		StringBuilder out = new StringBuilder();
		for(int i = 0; i < array.length; i++){
			for(int j = 0; j < array[0].length; j++){
				out.append(Arrays.toString(array[i][j]).replace("[", "").replace("]", "\n"));				
			}
		}
		return(out.toString());
	}
	
	public static void fillselected(double[][][] arraylong, int ind, double[][] array){
		for(int i = 0; i < arraylong[ind].length; i++){
			for(int j = 0; j < arraylong[ind][0].length; j++){
				arraylong[ind][i][j] = array[i][j];
			}
		}
	}
	
	public static void fillselected(int[][][] arraylong, int ind, int[][] array){
		for(int i = 0; i < arraylong[ind].length; i++){
			for(int j = 0; j < arraylong[ind][0].length; j++){
				arraylong[ind][i][j] = array[i][j];
			}
		}
	}

}
