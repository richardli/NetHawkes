package sampler;

import java.util.ArrayList;

public class Evaluation {
	/*
	 * function to calculate baseline rate ignoring process difference
	 */
	public static double Get_rate_homo(int[][] events, double T){
		double rate = 0;
		for(int i = 0; i < events.length; i++){
			for(int j = 0; j < events[0].length; j++){
				rate += events[i][j];
			}
		}
		return(rate / T / (events.length + 0.0));
	}
	/*
	 * function to calculate baseline log likelihood for multiple chains
	 */
//	public static double Get_loglik_homo_stupid(int[][] events_train, int[][] events, double dt){
//		double loglik = 0.0;
//		double one_rate = 0.0;
//		int Z = events_train.length * events_train[0].length;
//		for(int i = 0; i < events_train.length; i++){
//			for(int j = 0; j < events_train[0].length; j++){
//				one_rate += events_train[i][j];	
//			}
//		}
//		one_rate = one_rate / Z / dt;
//		for(int i = 0; i < events.length; i++){
//			loglik += Get_loglik_homo_single_stupid(events[i], one_rate, dt);
//		}
//		return(loglik);
//	}
//
//	/*
//	 * function to calculate baseline log likelihood
//	 */
//	public static double Get_loglik_homo_single_stupid( int[] events, double rate, double dt){
//		double loglik = 0.0;
//		for(int i = 0; i < events.length; i++){
//			loglik -=  org.apache.commons.math3.special.Gamma.logGamma(events[i] + 1);
//			loglik += events[i] * Math.log(rate);
//			loglik -= rate * dt;
//		}
//		return(loglik);
//	}
	/*
	 * function to calculate baseline log likelihood for multiple chains
	 */
	public static double Get_loglik_homo(int[][] events_train, int[][] events, double dt){
		double loglik = 0.0;
		for(int i = 0; i < events.length; i++){
			loglik += Get_loglik_homo_single(events_train[i], events[i], dt);
		}
		return(loglik);
	}

	/*
	 * function to calculate baseline log likelihood
	 */
	public static double Get_loglik_homo_single(int[] events_train, int[] events, double dt){
		double loglik = 0.0;
		int sum = 0;
		for(int i = 0; i < events_train.length; i++) sum += events_train[i];
		double rate = (sum + 0.0) / (events_train.length + 0.0) / dt;
		for(int i = 0; i < events.length; i++){
			loglik -=  org.apache.commons.math3.special.Gamma.logGamma(events[i] + 1);
			loglik += events[i] * Math.log(rate);
			loglik -= rate * dt;
		}
		return(loglik);
	}
	/*
	 * function to calculate likelihood
	 * Todo: figure out difference between Gamma likelihood for lambda and Poisson for events!
	 * Todo: check the use of dt in Linderman's codes 
	 */
	public static double Get_loglik_poisson_single(double[] lambda, int[] events, double dt){
		double loglik = 0.0;
		for(int i = 0; i < lambda.length; i++){
			loglik -= org.apache.commons.math3.special.Gamma.logGamma(events[i] + 1);
			loglik += events[i] * Math.log(lambda[i]);
			loglik -= lambda[i] * dt;
			if (Double.isNaN(loglik)) {
				System.out.println("Error in bin " + i);
				System.out.println("lambda in that cell: " + lambda[i]);
			}
		}
		return(loglik);
	}
	/*
	 * function to calculate rate
	 */
	public static double[][] getAggRate_with_history(int N, int K, int N_bin, int binmax, double tmax, 
			ArrayList<Event> data, 
			Baseline_model baseline, 
			Network_model network,
			Impulse_model impulse, 
			double dt, 
			ArrayList<Event> history, 
			double T_hist){

		double[][] rate = new double[K][N_bin];

		// initialize to background rate
		for(int i = 0; i < K; i++){
			for(int j = 0; j < N_bin; j++){
				rate[i][j] = baseline.lambda0[i][j];
			}
		}
		// add all historical impulses
		for(int i = 0; i < history.size(); i++){
			// if it is late enough
			if(history.get(i).timebin > T_hist - binmax){
				// for all process j
				for(int j = 0; j < K; j++){
					if(network.A[history.get(i).proc][j] == 0) continue;
					
					// calculate influence
					for(int l = 1; l < binmax; l++){
						// T_hist = 100, binmax = 10, 
						// then timebin ranges 0 - 99, first time bin to consider is:
						// timebin = 91, affect 92, ..., 99 and bin 0 in test
						// so we want to know 9 + 1 - (100 - 91) = 1 -> first cell (cell 0), add impulse 9
						int l_after = l + 1 - (int) (T_hist - history.get(i).timebin);
//						System.out.println(history.get(i).timebin);

						if(l_after > 0){
							double impulse_to_add = impulse.calculateImpulse((l+0.0)*dt, history.get(i).proc, j, tmax);
							rate[j][l_after - 1] += impulse_to_add* network.W[history.get(i).proc][j];
						}
					}
				}
			}
		}
		// add all events' impulses
		// for all event i
		for(int i = 0; i < N; i++){
			// for all process j
			for(int j = 0; j < K; j++){

				// if i-th event could not affect process j
				if(network.A[data.get(i).proc][j] == 0) continue;

				// if it could
				// for all following time points
				for(int l=1; l < binmax; l++){
					// if the bin count plus influence is larger than total bin
					int impacted_bin = data.get(i).timebin + l; 
					if(impacted_bin >= N_bin) break;
					double impulse_to_add = impulse.calculateImpulse((l+0.0)*dt, data.get(i).proc, j, tmax);
					rate[j][impacted_bin] += impulse_to_add* network.W[data.get(i).proc][j];
					//					if(j == 0 & impacted_bin == 2){
					//						System.out.println(i + " " + data.get(i).proc + " " + data.get(i).time);
					//					}
					if(rate[j][impacted_bin] < 0){
						System.out.println("Negative lambda!");
					}
				}
			}
		}			
		return(rate);

	}
	/*
	 * function to calculate rate
	 */
	public static double[][] getAggRate(int N, int K, int N_bin, int binmax, double tmax, 
			ArrayList<Event> data, 
			Baseline_model baseline, 
			Network_model network,
			Impulse_model impulse, 
			double dt){
		double[][] rate = new double[K][N_bin];

		// initialize to background rate
		for(int i = 0; i < K; i++){
			for(int j = 0; j < N_bin; j++){
				rate[i][j] = baseline.lambda0[i][j];
			}
		}
		// add all events' impulses
		// for all event i
		for(int i = 0; i < N; i++){
			// for all process j
			for(int j = 0; j < K; j++){

				// if i-th event could not affect process j
				if(network.A[data.get(i).proc][j] == 0) continue;

				// if it could
				// for all following time points
				for(int l=1; l < binmax; l++){
					// if the bin count plus influence is larger than total bin
					int impacted_bin = data.get(i).timebin + l; 
					if(impacted_bin >= N_bin) break;
					double impulse_to_add = impulse.calculateImpulse((l+0.0)*dt, data.get(i).proc, j, tmax);
					rate[j][impacted_bin] += impulse_to_add* network.W[data.get(i).proc][j];
					//					if(j == 0 & impacted_bin == 2){
					//						System.out.println(i + " " + data.get(i).proc + " " + data.get(i).time);
					//					}
					if(rate[j][impacted_bin] < 0){
						System.out.println("Negative lambda!");
					}
				}
			}
		}			
		return(rate);
	}

	/*
	 *  function to calculate single chain rate
	 */
	public static double[] getAggRate_onechain(int N, int K, int N_bin, int binmax, double tmax, 
			ArrayList<Event> data, 
			Baseline_model baseline, 
			Network_model network,
			Impulse_model impulse, 
			int k, 
			double dt){
		double[] rate = new double [N_bin];

		// initialize to background rate
		for(int j = 0; j < N_bin; j++){
			rate[j] = baseline.lambda0[k][j];
		}
		// add all events' impulses
		// for all event i
		for(int i = 0; i < N; i++){
			// if i-th event could not affect process k
			if(network.A[data.get(i).proc][k] == 0) continue;

			// if it could
			// for all following time points
			for(int l=1; l < binmax; l++){
				// if the bin count plus influence is larger than total bin
				int impacted_bin = data.get(i).timebin + l; 
				if(impacted_bin >= N_bin) break;
				double impulse_to_add = impulse.calculateImpulse((l+0.0)*dt, data.get(i).proc, k, tmax);
				rate[impacted_bin] += impulse_to_add* network.W[data.get(i).proc][k];
				if(rate[impacted_bin] < 0){
					System.out.println("Negative lambda!");
				}
			}
		}
		return(rate);
	}

	/*
	 *  function to calculate single chain rate
	 */
	public static double[] getAggRate_onechain_updatedAdj(int N, int K, int N_bin, int binmax, double tmax, 
			ArrayList<Event> data, 
			Baseline_model baseline, 
			Network_model network,
			Impulse_model impulse, 
			int k, 
			double dt, 
			int[][] A){
		double[] rate = new double [N_bin];

		// initialize to background rate
		for(int j = 0; j < N_bin; j++){
			rate[j] = baseline.lambda0[k][j];
		}
		// add all events' impulses
		// for all event i
		for(int i = 0; i < N; i++){
			// if i-th event could not affect process k
			if(A[data.get(i).proc][k] == 0) continue;

			// if it could
			// for all following time points
			for(int l=1; l < binmax; l++){
				// if the bin count plus influence is larger than total bin
				int impacted_bin = data.get(i).timebin + l; 
				if(impacted_bin >= N_bin) break;
				double impulse_to_add = impulse.calculateImpulse((l+0.0)*dt, data.get(i).proc, k, tmax);
				rate[impacted_bin] += impulse_to_add * network.W[data.get(i).proc][k];
				if(rate[impacted_bin] < 0){
					System.out.println("Negative lambda!");
				}
			}
		}
		return(rate);
	}

}
