package mathutil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class DistUtil {
	public static int discrete_sample(double[] probs, double u){
		int index = probs.length - 1;
		double sum = 0;
		double cumsum = 0;
		for(int i = 0 ; i < probs.length; i++) sum += probs[i];
		for(int i = 0; i < probs.length; i++){
			cumsum += probs[i];
			if(cumsum >= u * sum){
				index = i;
				break;
			}
		}
		return(index);
	}
	public static int discrete_sample(ArrayList<Double> probs, double u){
		int index = probs.size() - 1;
		double sum = 0;
		double cumsum = 0;
		for(int i = 0 ; i < probs.size(); i++) sum += probs.get(i);
		for(int i = 0; i < probs.size(); i++){
			cumsum += probs.get(i);
			if(cumsum >= u * sum){
				index = i;
				break;
			}
		}
		return(index);
	}
	
	public static void main(String[] args) {
		Random rand = new Random();
		double[] probvec = {0.1, 0.01, 0.3, 1.4};
		double[] normprob = new double[4];
		double Z = 0;
		ArrayList<Double> prob = new ArrayList<Double>();
		for(int i = 0; i < probvec.length; i++){
			prob.add(probvec[i]);
			Z += probvec[i];
		}
		for(int i = 0; i < probvec.length; i++){
			normprob[i] = probvec[i] / Z;
		}
		int[] count = new int[4];
		for(int i = 0; i < 10000; i++){
			int a = DistUtil.discrete_sample(prob, rand.nextDouble());
		    count[a]++;
		}
		System.out.println(Arrays.toString(normprob));
		System.out.println(Arrays.toString(count));
	}
	
}
