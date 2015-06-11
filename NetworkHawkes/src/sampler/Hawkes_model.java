package sampler;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import mathutil.VecUtil;
import cern.jet.random.tdouble.Beta;
import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
//import org.apache.commons.math3.special.Gamma;

public class Hawkes_model {
	/*
	 * variables initialized in read_data
	 */
	// number of events
	public int N;
	// number of processes
	public int K;
	// max time
	public double T;
	// max time including held-out events
	public double Tall;
	public int N_test;
	public int N_bin_test;

	// data info
	ArrayList<Event> data = new ArrayList<Event>();
	ArrayList<Event> testset = new ArrayList<Event>();
	Network_model network;
	Impulse_model impulse;
	Baseline_model baseline;

	/*
	 * variables initialized in binning
	 */
	// max impulse memory
	public double tmax;	
	// number of bins
	public int N_bin;
	// bin time interval
	public double dt;
	// bin time max memory
	public int binmax;
	// number of non-zero bins
	public int N2;
	// events in matrix form
	public int[][] events = new int[K][N_bin];		
	// aggregated rates
	double[][] aggRate = new double[K][N_bin];
	// events for test set
	public int[][] events_test = new int[K][N_bin_test];

	/*
	 * variables not initialized
	 */


	/*
	 * parameters specified directly
	 */
	public double alpha_w = 2.0 ;
	public double beta_w = 1.0; 	// 1.0 ----> this gets updated in algorithm
	public double alpha_lam = 1.0;
	public double beta_lam = 1.0;
	public double tau0 = 1;
	public double tau1 = 2;
	public double mu0 = -1.0; 
	public double kap0 = 10.0; // 10.0
	public double alpha0 = 10.0; // 10.0
	public double beta0 = 1.0;
	public int seed = 0;



	/*
	 * Main function
	 */
	public static void main(String[] args) throws IOException{
		// initialization
		int nitr = 1000, seed = 1;
		boolean isER = true, isLatent = false;
		String header = null;
		String runname;
		String[] files = new String[2];

		boolean start_from_zero = false;
		double T = 900;
		double Tall = 1000;
		int K = 30;
		int N_bin = 900;
		double tmax = 10.0;
		boolean allow_self_connect = false;
		int M = 41;
		double Period = 3600;

		double network_ini = 0.0;
		double baseline_ini = -1;
		boolean skip_W = true;
		boolean skip_W_hyper = false;
		double W_fix = 0.8;
		boolean skip_theta = true;
		boolean skip_lambda = false;
		boolean skip_A = false;
		boolean show_diag = false;
		boolean sample_A_block = true;
		int how_many_flip = 1;
		boolean const_baseline = false;
		boolean same_z_per_cell = true;
		boolean skip_non_used_w = false;
		boolean skip_all_network = false;
		double sigma_LGCP = 1;
		
		// default set up in local machine
		if(args.length == 0){
			nitr = 1000;
			isER = true;
			isLatent = false;
			//						header = "/Users/zehangli/Dropbox/UW_grad/518-15spr/ImplicitNetwork/NetworkHawkes/data/";					
			//						runname = "rhawkes-sim/seed3-net1";
			//						files[0] = header + runname + "-time.txt";
			//						files[1] = header + runname + "-proc.txt";
			header = "/Users/zehangli/Dropbox/UW_grad/518-15spr/ImplicitNetwork/NetworkHawkes/data/";					
			runname = "sp100/event";
			files[0] = header + runname + "-time.txt";
			files[1] = header + runname + "-proc.txt";  

			T = 92856;
			Tall = 116071;
			K = 100;
			N_bin = 92856;
			tmax = 12.0;	

			M = 41;
			Period = 23214;
			const_baseline = false;

			// set up in cluster
		}else{
			header = "/homes/lizehang/NetworkHawkes/data/";
			runname = args[0];
			files[0] = header + runname + "-time.txt";
			files[1] = header + runname + "-proc.txt";
			isER = Boolean.parseBoolean(args[1]);
			isLatent = Boolean.parseBoolean(args[2]);
			T = Double.parseDouble(args[3]);
			Tall = Double.parseDouble(args[4]);
			K = Integer.parseInt(args[5]);
			N_bin = Integer.parseInt(args[6]);
			tmax = Double.parseDouble(args[7]);
			allow_self_connect = Boolean.parseBoolean(args[8]);
			network_ini = Double.parseDouble(args[9]);
			baseline_ini = Double.parseDouble(args[10]);

			skip_W = Boolean.parseBoolean(args[11]);
			skip_W_hyper = Boolean.parseBoolean(args[12]);
			W_fix = Double.parseDouble(args[13]);

			skip_theta = Boolean.parseBoolean(args[14]);
			skip_lambda = Boolean.parseBoolean(args[15]);

			skip_A = Boolean.parseBoolean(args[16]);
			nitr = Integer.parseInt(args[17]);
			show_diag = Boolean.parseBoolean(args[18]);

			sample_A_block = Boolean.parseBoolean(args[19]);
			M = Integer.parseInt(args[20]);
			Period = Integer.parseInt(args[21]);

			how_many_flip = Integer.parseInt(args[22]);
			seed = Integer.parseInt(args[23]);
			const_baseline = Boolean.parseBoolean(args[24]);
			same_z_per_cell = Boolean.parseBoolean(args[25]);
			skip_non_used_w = Boolean.parseBoolean(args[26]);
			// for a few added parameters
			if(args.length > 27){
				skip_all_network =  Boolean.parseBoolean(args[27]);
				sigma_LGCP = Double.parseDouble(args[28]);
			}

		}

		DoubleRandomEngine rngEngine=new DoubleMersenneTwister(seed);
		Normal rngN=new Normal(0.0,1.0,rngEngine);
		Gamma rngG=new Gamma(1.0,1.0,rngEngine);
		Beta rngB = new Beta(1.0, 1.0, rngEngine);
		Random rand = new Random(seed);



		Hawkes_model sampler = new Hawkes_model();
		// initialize
		Helper.read_data(sampler, files, T, Tall, K, start_from_zero);
		Helper.binning(sampler, N_bin, tmax);

		System.out.println("Total number of events   : " + sampler.N);
		System.out.println("Total time               : " + sampler.T);
		System.out.println("Total number of processes: " + sampler.K);
		System.out.println("Total number of bins     : " + sampler.N_bin);
		System.out.println("Time width for each bin  : " + sampler.dt);
		System.out.println("Maximum time of influence: " + sampler.tmax);


		// more initializations for each model
		sampler.network.setSelfConnect(allow_self_connect);

		if(network_ini == 0){
			sampler.network.initempty(rngG);			
		}else if(network_ini == 1){
			sampler.network.initfull(rngG, allow_self_connect);
		}else{
			sampler.network.initp(network_ini, rand, rngG);			
		}

		if(!const_baseline){
			// initial to Periodic kernel				
			sampler.baseline.initLGCP(sampler.events, sampler.T, M, Period, sigma_LGCP);
			//			sampler.baseline.lambda_LGCP_diag();			
		}else{
			// initial to homo rate
			//		 sampler.baseline.initconst(sampler.events, sampler.T, baseline_ini);
			// initial to prior sample
			sampler.baseline.initconst(sampler.events, sampler.T, rngG);
		}

		// save variables
		int nburn = (int)(nitr / 2);
		int[][][] A_out = new int[nitr - nburn][sampler.K][sampler.K];
		double[][][] W_out = new double[nitr - nburn][sampler.K][sampler.K];
		double[][][] P_out = new double[nitr - nburn][sampler.K][sampler.K];
		double[][][] baseline_out = new double[nitr - nburn][sampler.K][sampler.baseline.M];

		double[] lik_train = new double[nitr - nburn];
		double[] lik_test = new double[nitr - nburn];
		double[] bits_per_spike_train = new double[nitr - nburn];
		double[] bits_per_spike_test = new double[nitr - nburn];
		double[] bits_per_second_test = new double[nitr - nburn];

		double lik0 = 0;
		double lik1 = 0;
		double bps0 = 0; 
		double bps1 = 0;
		double bps2 = 0;
		// baseline rate
		double lik_homo_train = Evaluation.Get_loglik_homo(sampler.events, sampler.events, sampler.dt);
		double lik_homo_test = Evaluation.Get_loglik_homo(sampler.events, sampler.events_test, sampler.dt);
		//		double lik_homo_train = Evaluation.Get_loglik_homo_stupid(sampler.events, sampler.events, sampler.dt);
		//		double lik_homo_test = Evaluation.Get_loglik_homo_stupid(sampler.events, sampler.events_test, sampler.dt);

		//		double rate_homo_train = Evaluation.Get_rate_homo(sampler.events, T);
		//		double rate_homo_test = Evaluation.Get_rate_homo(sampler.events_test, Tall - T);

		for(int i = 0; i < nitr; i++){
			if(i > 50){
				System.out.println("checkpoint");
			}
			System.out.println("itr " + i);
			if(i >= nitr - nburn){
				VecUtil.fillselected(A_out, i - nitr + nburn, sampler.network.A);
				VecUtil.fillselected(W_out, i - nitr + nburn, sampler.network.W);
				VecUtil.fillselected(P_out, i - nitr + nburn, sampler.network.P);
				VecUtil.fillselected(baseline_out, i - nitr + nburn, sampler.baseline.rate_dense);

				// metrics
				lik_train[i - nitr + nburn] = lik0;
				lik_test[i - nitr + nburn] = lik1;
				bits_per_spike_train[i - nitr + nburn] = bps0;
				bits_per_spike_test[i - nitr + nburn] = bps1;
				bits_per_second_test[i - nitr + nburn] = bps2;

			}

			// re-count configurations
			if(!skip_all_network){
				sampler.network.Count_config(sampler.N, sampler.data);
				if(show_diag) sampler.network.diag_config();	
				// compare with MF configuration
				//			sampler.network.Count_config_MF(sampler.N, sampler.data);
				//			if(show_diag) sampler.network.diag_config();
			}



			// update theta
			if(!skip_all_network){
				sampler.impulse.Sample_theta(sampler.data, rngG, rngN, skip_theta);
				if(show_diag) sampler.impulse.diag_theta();
			}


			// this two steps relies heavily on configuration count
			if(i != 0){
				// update lambda
				if(const_baseline){
					sampler.baseline.Sample_lambda_const(rngG, sampler.network, skip_lambda);
					//					if(show_diag) sampler.baseline.lambda_const_diag();					
				}else{
					sampler.baseline.sampleLGCP(sampler.network, rand, rngN, sampler.N, sampler.data);
					//					sampler.baseline.sampleLGCP_fixed(sampler.network, rand);
					//					sampler.baseline.lambda_LGCP_diag();					
				}

				// update P					
				if(!skip_all_network){
					sampler.network.Sample_P(rngB, isER, isLatent);
					if(show_diag) sampler.network.diag_p();
				}
				// update W
				if(!skip_all_network){		
					sampler.network.Sample_W(rngG, skip_W, skip_W_hyper, W_fix, skip_non_used_w);
					if(show_diag) sampler.network.diag_w();
				}

			}


			// update A
			if(!skip_all_network){

				if(!sample_A_block){
					sampler.network.Sample_A(sampler.N, sampler.data, sampler.tmax, sampler.N_bin,
							sampler.binmax, sampler.baseline, sampler.network, 
							sampler.dt, sampler.impulse, rand, sampler.events, skip_A);				
				}else{
					sampler.network.Sample_A_block(sampler.N, sampler.data, sampler.tmax, sampler.N_bin,
							sampler.binmax, sampler.baseline, sampler.network, 
							sampler.dt, sampler.impulse, rand, sampler.events, skip_A, how_many_flip);				
				}
				if(show_diag) sampler.network.diag_A();
			}

			// update z
			if(!skip_all_network){
				//			Helper.Sample_Z_MAP(sampler, rand);
				Helper.Sample_Z(sampler, rand);
				if(same_z_per_cell) Helper.unify_Z(sampler.data, K, N_bin);
			}
			/*
			 * Calculate rate and likelihood
			 */
			double[][] currentRate = Evaluation.getAggRate(sampler.N, sampler.K, 
					sampler.N_bin, sampler.binmax, sampler.tmax, sampler.data, 
					sampler.baseline, sampler.network, sampler.impulse, sampler.dt);
			double[][] predictedRate = Evaluation.getAggRate_with_history(sampler.N_test, sampler.K, 
					sampler.N_bin_test, sampler.binmax, sampler.tmax, sampler.testset, 
					sampler.baseline, sampler.network, sampler.impulse, sampler.dt, 
					sampler.data, sampler.T);
			double[] metrics = Helper.getMetric(K, currentRate, predictedRate, 
					sampler.events, sampler.events_test, sampler.dt, 
					lik_homo_train, lik_homo_test, T, Tall);
			lik0 = metrics[0];
			lik1 = metrics[1];
			bps0 = metrics[2];
			// bits per rate
			bps1 = metrics[3];
			// bits per second
			bps2 = metrics[4];
		}
		runname = runname + seed;
		// finish and output to file
		Helper.save_output(header, runname, A_out, W_out, P_out, baseline_out, lik_train, lik_test, 
				bits_per_spike_train, bits_per_spike_test, bits_per_second_test);		
	}
}
