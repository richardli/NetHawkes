package sampler;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import mathutil.DistUtil;
import mathutil.MathUtil;
import mathutil.VecUtil;

public class Helper {
	/*
	 * When reading data, the raw time stamps are used instead of bins
	 * start_from_zero: if the process identity starts from 0
	 * Initialize K, T, N, data 
	 */
	public static void read_data(Hawkes_model model, String[] paths, double T, double Tall, int K, boolean start_from_zero) throws IOException{
		model.K = K;
		model.T = T;
		model.Tall = Tall;
		model.network = new Network_model(K);
		model.network.pass_hyper(model.alpha_w, model.beta_w, model.tau0, model.tau1);
		model.impulse = new Impulse_model(K);
		model.impulse.pass_hyper(model.alpha0, model.beta0, model.kap0, model.mu0);

		// read time file
		BufferedReader sc = new BufferedReader(new FileReader(paths[0]));
		String line = sc.readLine();
		String[] field = line.split(",");					
		model.N = field.length;
		model.N_test = 0;
		// index vector for test set, -1 if is train
		int[] test_index = new int[field.length];
		
		// populate data list
		for(int i = 0; i < field.length; i++){
			double time = Double.parseDouble(field[i]);
			if(time <= T){
				model.data.add(new Event(time));		
				test_index[i] = -1;
			}else{
				model.testset.add(new Event(time));
				test_index[i] = model.testset.size() - 1;
				// remove one count from N, add to N_test
				model.N --;
				model.N_test ++;
			}
		}
		sc.close();

		// read process identity file
		sc = new BufferedReader(new FileReader(paths[1]));
		line = sc.readLine();
		field = line.split(",");					
		for(int i = 0; i < field.length; i++){
			// populate train set
			if(test_index[i] == -1){
				if(start_from_zero){
					model.data.get(i).setProc(Integer.parseInt(field[i]));
				}else{
					model.data.get(i).setProc(Integer.parseInt(field[i]) - 1);
				}				
			}else{
				if(start_from_zero){
					model.testset.get(test_index[i]).setProc(Integer.parseInt(field[i]));
				}else{
					model.testset.get(test_index[i]).setProc(Integer.parseInt(field[i]) - 1);
				}
			}
		}
		sc.close();
	}

	/*
	 * Initialize potential parents and matrix events (raw s[] vector is used)
	 */
	public static void binning(Hawkes_model model, int N_bin, double tmax){
		// initiate global parameters
		model.N_bin = N_bin;
		model.dt = model.T / (model.N_bin + 0.0);
		model.N_bin_test = (int) ((model.Tall - model.T) / model.dt);
		model.events_test = new int[model.K][model.N_bin_test];
		
		model.tmax = tmax;
		model.binmax = (int)  (model.tmax / (model.dt + 0.0));
		model.N2 = 0;
		model.events = new int[model.K][N_bin];
		model.baseline = new Baseline_model(model.K, model.T, N_bin);
		model.baseline.pass_hyper(model.alpha_lam, model.beta_lam);
		// message
		System.out.println("binning in process:");
		
		// populate event matrix
		for(int i = 0; i < model.N; i++){
				// update bintime for each data point
			model.data.get(i).binning(model.dt, model.T);

			// update event matrix
			model.events[model.data.get(i).proc][model.data.get(i).timebin]++;
			// update possible parent list
			double t1 = model.data.get(i).time;

			for(int j = 0; j < model.N; j++){
				double t0 = model.data.get(j).time;
				if(t0 < t1 && t0 > t1 - model.tmax){
					model.data.get(i).addPotentialParent(j);
				}					
			}
		}
		
		// binning for test set
		for(int i = 0; i < model.N_test; i++){
			model.testset.get(i).binning(model.dt, model.Tall, model.T);
			model.events_test[model.testset.get(i).proc][model.testset.get(i).timebin]++;
		}
	}

	/*
	 * extract k-k' pair index (k0: parent, k1: event)
	 */
	//		public ArrayList<Integer> getKKP(int k0, int k1){
	//			ArrayList<Integer> list = new ArrayList<Integer>();
	//			for(int i = 0; i < model.N; i++){
	//				if(model.data.get(i).proc == k1 && 
	//				   model.data.get(i).pc == k0) list.add(i);
	//			}
	//			return(list);
	//		}

	/* 
	 * function to sample Z
	 * TODO: not finished
	 */
	public static void Sample_Z(Hawkes_model model, Random rand){
		// reset mean field count to zero
		model.network.resetMF();
		// use an array to track which cell has been calculated
		// code: 0 -> not calculated, 1 -> calculated
		int[][] cell = new int[model.K][model.N_bin];	
		
		for(int i = 0; i < model.N; i++){
			// if no parent possible
			if(model.data.get(i).plist.size() == 0) model.data.get(i).parent = -1;
			// if parent possible
			if(model.data.get(i).plist.size() > 0){
				ArrayList<Integer> choices = new ArrayList<Integer>();
				choices.add(-1);
				choices.addAll(model.data.get(i).plist);
				// populate multi probabilities
				ArrayList<Double> probs = new ArrayList<Double>();
				
				for(int ch : choices){
					if(ch == -1){
						// TODO: check here!!!
						double rate_tmp = model.baseline.lambda0[model.data.get(i).proc][model.data.get(i).timebin];
						probs.add(rate_tmp);
					}else{
						// TODO: check here!!!
						int kk0 = model.data.get(ch).proc;
						int kk1 = model.data.get(i).proc;
						if(model.network.A[kk0][kk1] == 0){
							probs.add(0.0); 
						}else{
							// if connected
							double time_since = model.data.get(i).time - model.data.get(ch).time;
							double impulse_to_add = model.impulse.calculateImpulse(time_since, kk0, kk1, model.tmax);
//							System.out.println(kk1 + "  " + time_since + "   " + model.network.W[kk0][kk1] + " " + impulse_to_add);
							probs.add(model.network.W[kk0][kk1] * impulse_to_add);	
							if(cell[kk1][model.data.get(i).timebin] == 0){
								// note the reverse order
								model.network.MeanField[kk1][kk0 + 1] += model.network.W[kk0][kk1] * impulse_to_add * model.dt;
							}
						}
					}
				}
				// when finished loop, always make sure the cell is marked with 1 to avoid over-counting
				cell[model.data.get(i).proc][model.data.get(i).timebin] = 1;
				
				// sample from them
				int index_choosen = DistUtil.discrete_sample(probs, rand.nextDouble());
				// if the first index (background) is choosen
				if(index_choosen == 0){
					model.data.get(i).updateParent();
				}else{
					model.data.get(i).updateParent(choices.get(index_choosen), 
							model.data.get(choices.get(index_choosen)), 
							model.tmax); 											
				}
			}
		}
		
		// final update of background rate to MF, and normalize to probability
		for(int i = 0; i < model.K; i++){
			for(int j = 0; j < model.N_bin; j++){
				model.network.MeanField[i][0] += model.baseline.lambda0[i][j] * model.dt; 
			}
			model.network.MeanField[i] = VecUtil.norm(model.network.MeanField[i]); 
		}
	}
	
	public static void Sample_Z_MAP(Hawkes_model model, Random rand){
		for(int i = 0; i < model.N; i++){
			// if no parent possible
			if(model.data.get(i).plist.size() == 0) model.data.get(i).parent = -1;
			// if parent possible
			if(model.data.get(i).plist.size() > 0){
				ArrayList<Integer> choices = new ArrayList<Integer>();
				choices.add(-1);
				choices.addAll(model.data.get(i).plist);
				// populate multi probabilities
				double prob0 = 0;
				double[] probs = new double[model.K];
				int[] fake_parent = new int[model.K];
				
				for(int ch : choices){
					if(ch == -1){
						// TODO: check here!!!
						prob0 = (model.baseline.lambda0[model.data.get(i).proc][model.data.get(i).timebin]);
					}else{
						// TODO: check here!!!
						int kk0 = model.data.get(ch).proc;
						int kk1 = model.data.get(i).proc;
						if(model.network.A[kk0][kk1] == 0){
							probs[kk0] += 0; 
						}else{
							// if connected
							double time_since = model.data.get(i).time - model.data.get(ch).time;
							double impulse_to_add = model.impulse.calculateImpulse(time_since, kk0, kk1, model.tmax);
//							System.out.println(kk1 + "  " + time_since + "   " + model.network.W[kk0][kk1] + " " + impulse_to_add);
							probs[kk0] += (model.network.W[kk0][kk1] * impulse_to_add);								
							fake_parent[kk0] = ch;
						}
					}
				}
				// sample from them
				//int index_choosen = DistUtil.discrete_sample(probs, rand.nextDouble());
				double prob1 = 0;
				for(int k = 0; k < model.K; k++){
					prob1 += probs[k];
				}
				
				int index_choosen;
				if(prob0 / (prob0 + prob1) > rand.nextDouble()){
					index_choosen = 0;
				}else{
					index_choosen = MathUtil.whichmax(probs) + 1;
				}
				// if the first index (background) is choosen
				if(index_choosen == 0){
					model.data.get(i).updateParent();
				}else{
					model.data.get(i).updateParent(fake_parent[index_choosen - 1], 
							model.data.get(fake_parent[index_choosen - 1]), 
							model.tmax); 											
				}
			}
		}
	}
	
	
	/*
	 * function to unify all binned events to have the same parent
	 *  note: just a hack to make it like the codes of Linderman
	 */
	public static void unify_Z(ArrayList<Event> data, int K, int N_bin){
		int[][] first_parent = new int[K][N_bin];
		int[][] first_kid = new int[K][N_bin];
		
		for(int i = 0; i < K; i++){
			for(int j = 0; j < N_bin; j++){
				first_parent[i][j] = -2;
			}
		}
		for(int i = 0; i < data.size(); i++){
			if(first_parent[data.get(i).proc][data.get(i).timebin] == -2){
				// if it is the first event in that cell
				first_parent[data.get(i).proc][data.get(i).timebin] = data.get(i).parent;
				first_kid[data.get(i).proc][data.get(i).timebin] = i;
			}else{
				// if it is not the first event
				data.get(i).copyParent(data.get(first_kid[data.get(i).proc][data.get(i).timebin]));
			}
		}
				
	}
	
	/*
	 * function to save output
	 */

	public static void save_output(String header, String runname, 
			int[][][] A_out, double[][][] W_out, double[][][] P_out, double[][][] baseline_out,
			double[] lik_train, double[] lik_test,	double[] bits_per_spike_train,
			double[] bits_per_spike_test, double[] bits_per_second_test) throws IOException{
			
				String outfile1 = header + runname + "A-out.txt";
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile1));
				bw.write(VecUtil.write3d(A_out));
				bw.close();
			
				String outfile2 = header + runname +  "W-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile2));
				bw.write(VecUtil.write3d(W_out));
				bw.close();
				
				String outfile3 = header + runname + "lik-train-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile3));
				bw.write(Arrays.toString(lik_train).replace("[", "").replace("]", ""));
				bw.close();
				
				String outfile4 = header + runname + "lik-test-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile4));
				bw.write(Arrays.toString(lik_test).replace("[", "").replace("]", ""));
				bw.close();
				
				String outfile5 = header + runname + "bps-train-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile5));
				bw.write(Arrays.toString(bits_per_spike_train).replace("[", "").replace("]", ""));
				bw.close();
				
				String outfile6 = header + runname + "bps-test-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile6));
				bw.write(Arrays.toString(bits_per_spike_test).replace("[", "").replace("]", ""));
				bw.close();
				
				String outfile7 = header + runname + "bps-second-test-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile7));
				bw.write(Arrays.toString(bits_per_second_test).replace("[", "").replace("]", ""));
				bw.close();
				
				String outfile8 = header + runname + "P-out.txt";
				 bw = new BufferedWriter(new FileWriter(outfile8));
				bw.write(VecUtil.write3d(P_out));
				bw.close();
			
				String outfile9 = header + runname +  "baseline-out.txt";
				bw = new BufferedWriter(new FileWriter(outfile9));
				bw.write(VecUtil.write3d(baseline_out));
				bw.close();
				
	}
	
	public static double[] getMetric(int K, double[][] currentRate,
			double[][] predictedRate, int[][] events, int[][] events_test, 
			double dt, 
			double lik_homo_train,
			double lik_homo_test,
			double T, double Tall){
		double[] metric = new double[5];

		for(int kk = 0; kk < K; kk++){
			metric[0] += Evaluation.Get_loglik_poisson_single(currentRate[kk], 
					events[kk],  dt);
			metric[1] += Evaluation.Get_loglik_poisson_single(predictedRate[kk], 
					events_test[kk], dt);
		}
		metric[2] = (metric[0] - lik_homo_train) / (VecUtil.sum2d(events) +0.0);
		//	sampler.K / sampler.N_bin_test; 
		// rate_homo_train; //(sampler.N+0.0);
		metric[3] = (metric[1] - lik_homo_test) / (VecUtil.sum2d(events_test) +0.0);
		// bits per second
		metric[4] = (metric[1] - lik_homo_test) / (Tall - T - 1.0);

		//  sampler.K / sampler.N_bin_test;
		//rate_homo_test;//(sampler.N_test+0.0);
		System.out.println("::::::::   LIKELIHOOD   ::::::::");
		System.out.println("Training Loglik :   " + metric[0] + "   Baseline: " + lik_homo_train);
		System.out.println("Testing Loglik  :   " + metric[1] + "   Baseline: " + lik_homo_test);
		System.out.println("Training BPS    :   " + metric[2]);
		System.out.println("Testing BPS     :   " + metric[3]);
		System.out.println("Testing BPSecond:   " + metric[4]);
		return(metric);
	}
	
}
