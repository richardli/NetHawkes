package mathutil;

public class MathUtil {
	public static double logSumOfExponentials(double x1, double x2) {
		double max = x1;
		if(x2 > x1) max = x2;
		double sum = 0.0;
		if (x1 != Double.NEGATIVE_INFINITY)
			sum += java.lang.Math.exp(x1 - max);
		if (x2 != Double.NEGATIVE_INFINITY)
			sum += java.lang.Math.exp(x2 - max);
		return max + java.lang.Math.log(sum);
	}
	
	public static int whichmax(double[] x){
		int which = 0;
		double max = x[0];
		for(int i = 0; i < x.length; i++){
			if(x[i] > max){
				max = x[i];
				which = i;
			}
		}
		return(which);
	}
	
	public static int max(int[] x){
		int max = x[0];
		for(int i = 0; i < x.length; i++){
			if(x[i] > max){
				max = x[i];
			}
		}
		return(max);
	}
	
	public static int min(int[] x){
		int min = x[0];
		for(int i = 0; i < x.length; i++){
			if(x[i] < min){
				min = x[i];
			}
		}
		return(min);
	}
}
