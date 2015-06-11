package sampler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import mathutil.MathUtil;
import mathutil.VecUtil;
import cern.jet.random.tdouble.Beta;
import cern.jet.random.tdouble.Gamma;

public class Network_model {
	// number of processes
	public int K;
	// A matrix
	public int[][] A = new int[K][K];
	// W matrix
	public double[][] W = new double[K][K];
	// P matrix
	public double[][] P = new double[K][K];
	// intermediate variables
	public double[][] Nkk = new double[K][K];
	public double[] Nk = new double[K];
	public double[] N0k = new double[K];	

	// mean field for network structure
	// MeanField[1][0] <- expected event from baseline on process 1
	// MeanField[1][5] <- expected event from process 4 to process 1
	double[][] MeanField = new double[K][K + 1];
	boolean self;

	//  hyper-parameters
	double alpha_w;
	double beta_w;
	double tau0;
	double tau1;

	public Network_model(int K){
		this.K = K;
		this.A = new int[K][K];
		this.W = new double[K][K];
		this.P = new double[K][K];
		this.Nkk = new double[K][K];
		this.Nk = new double[K];
		this.N0k = new double[K];	
		this.MeanField = new double[K][K + 1];
	}

	public void pass_hyper(double alpha_w, 
			double beta_w, 
			double tau0,
			double tau1){
		this.alpha_w = alpha_w;
		this.beta_w = beta_w;
		this.tau0 = tau0;
		this.tau1 = tau1;
	}
	public void setSelfConnect(boolean self){
		this.self = self;
	}
	public void resetMF(){
		this.MeanField = new double[this.K][this.K + 1];
	}
	public void initfull(Gamma rngG, boolean self_connect){
		for(int i = 0 ; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				this.P[i][j] = 1;
				this.A[i][j] = 1;
				this.W[i][j] = rngG.nextDouble(this.alpha_w, this.beta_w);
				if( (!self_connect) & i == j){
					this.P[i][j] = 0;
					this.A[i][j] = 0;
				}
			}
		}
	}
	public void initempty(Gamma rngG){
		for(int i = 0 ; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				this.P[i][j] = 0;
				this.A[i][j] = 0;
				this.W[i][j] = rngG.nextDouble(this.alpha_w, this.beta_w);
			}
		}
	}
	public void initp(double p, Random rand, Gamma rngG){
		for(int i = 0 ; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				this.P[i][j] = p;
				this.A[i][j] = (rand.nextDouble() > p)? 0 : 1;
				this.W[i][j] = rngG.nextDouble(this.alpha_w, this.beta_w);
			}
		}

		if(!this.self){
			for(int i = 0 ; i < this.K; i++){
				this.P[i][i] = 0; 
				this.A[i][i] = 0;
			}
		}	
	}

	/*
	 * function to count configurations
	 */
	public void Count_config(int N, ArrayList<Event> data){
		this.Nkk = new double[K][K];
		this.Nk = new double[K];
		this.N0k = new double[K];

		for(int i = 0; i < N; i++){
			// add to Nk
			this.Nk[data.get(i).proc]++;
			// add to N0k and Nkk
			if(data.get(i).parent == -1){
//				System.out.println(data.get(i).proc + "on process " + data.get(i).pc);
				this.N0k[data.get(i).proc]++;
			}else{
//				System.out.println(data.get(i).proc + "on process " + data.get(i).pc);
				this.Nkk[data.get(i).pc][data.get(i).proc]++;	
			}
//			if(data.get(i).pc != -1 & data.get(i).parent == -1){
//				System.out.println("Someone forgot to update parent process");
//			}
		}
	}
	/*
	 * count baseline events by time bin
	 */
	/*
	 * count dense events
	 */
	public int[] countBaseline(int N, ArrayList<Event> data, int M, int N_bin, int ncell, int k){
		int[] events_base = new int[N_bin];
		for(int i = 0; i < N; i++){
			if(data.get(i).parent == -1 & data.get(i).proc == k){
				events_base[i] ++;				
			}
		}
		return(events_base);
	}
	/*
	 * count dense baseline events
	 */
	public int[] countDense(int N, ArrayList<Event> data, int M, int N_bin, int ncell, int k){
		int[] events_dense = new int[M];
		for(int i = 0; i < N; i++){
			if(data.get(i).parent == -1 & data.get(i).proc == k){
				int bin = data.get(i).timebin;
				int m = (int)(bin / (ncell + 0.0));
				events_dense[m] ++;				
			}
		}
		return(events_dense);
	}
	/*
	 * function to count configurations by mean field 
	 */
	public void Count_config_MF(int N, ArrayList<Event> data){
		this.Nkk = new double[K][K];
		//		this.Nk = new double[K];
		this.N0k = new double[K];

		for(int i = 0; i < K; i++){
			this.N0k[i] = this.MeanField[i][0] * this.Nk[i];
			for(int j = 0; j < K ; j++){
				// notice the order of the last equation
				this.Nkk[j][i] = this.MeanField[i][j+1] * this.Nk[j];
			}
		}
	}

	/*
	 * function to sample W
	 */
	public void Sample_W(Gamma rngG, boolean skip_update_para, boolean skip_hyper, double fixW, boolean skip_non_used_w){

		// re-sample hyper parameter
		if(!skip_hyper){
			this.beta_w = rngG.nextDouble(this.K * this.K * this.alpha_w, 
								VecUtil.sum2d(this.W));
			// error in paper???
//			double newshape = VecUtil.sum2d(this.A) * this.alpha_w;
//			double newrate = VecUtil.sum2d(this.W);
//			this.beta_w = rngG.nextDouble(newshape, newrate);	
		}
		if(skip_update_para){
			for(int i = 0; i < this.K; i++){
				for(int j = 0; j < this.K; j++){
					if(fixW > 0){
						this.W[i][j] = fixW;
					}else{
						this.W[i][j] =rngG.nextDouble(this.alpha_w , this.beta_w);
					}
				}
			}
			return;
		}
		// re-sample W
		for(int i = 0; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				// avoid sampling W under A=0
				// it makes W small, then when A flips, very few Z could follow
				if(skip_non_used_w & this.A[i][j] == 0) continue;
				double shape = this.alpha_w + this.Nkk[i][j] * this.A[i][j];
				double rate = this.beta_w + this.Nk[i] * this.A[i][j];
				this.W[i][j] = rngG.nextDouble(shape, rate);
//				if(this.Nkk[i][j] > 0){
//					System.out.println("send :" +this.Nkk[i][j] +" from: " +this.Nk[i]);
//					System.out.println("W :" +this.W[i][j] +" mean: " +this.Nkk[i][j]/(this.Nk[i] + 0.0));
//					System.out.println(" ");
//				}
			}
		}		
	}
	
	
	/*
	 * new function to sample W
	 */
	
	public void Sample_W_test(Gamma rngG, boolean skip_update_para, boolean skip_hyper, double fixW, boolean skip_non_used_w){

		// re-sample hyper parameter
		if(!skip_hyper){
			this.beta_w = rngG.nextDouble(this.K * this.K * this.alpha_w, 
								VecUtil.sum2d(this.W));
		}
		if(skip_update_para){
			for(int i = 0; i < this.K; i++){
				for(int j = 0; j < this.K; j++){
					if(fixW > 0){
						this.W[i][j] = fixW;
					}else{
						this.W[i][j] =rngG.nextDouble(this.alpha_w , this.beta_w);
					}
				}
			}
			return;
		}
		// re-sample W
		for(int i = 0; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				// avoid sampling W under A=0
				// it makes W small, then when A flips, very few Z could follow
				if(skip_non_used_w & this.A[i][j] == 0) continue;
				double shape = this.alpha_w;
				double rate = this.beta_w;
				if(A[i][j] == 0){
					shape += this.Nkk[i][j];
					rate += this.Nk[i];
				}else{
					shape += this.Nk[j] - this.N0k[j];
					for(int kk = 0; kk < this.K; kk++) rate += this.Nk[kk] * this.A[kk][j];
				}
				this.W[i][j] = rngG.nextDouble(shape, rate);
			}
		}		
	}
	/*
	 * function to sample P matrix
	 */
	public void Sample_P(Beta rngB, boolean isER, boolean isLatent){
		// Erdos-Renyi model
		if(isER){
			double non_zero = VecUtil.sum2d(this.A);
			double tau0p = this.tau0 + non_zero;
			double tau1p = this.tau1 + this.K * this.K - non_zero;
			double rho = rngB.nextDouble(tau0p, tau1p);
			// populate all with rho
			for(int i = 0; i < this.K; i++){
				for(int j = 0; j < this.K; j++){				
					this.P[i][j] = rho;
				}
			}
			return;
		}

		// Latent distance model
		if(isLatent){
			for(int i = 0; i < this.K; i++){
				// for smaller index
				for(int j = 0; j < i; j++){
					this.P[i][j] = this.P[j][i];
				}
				// for larger index
				for(int j = i; j < this.K; j++){				
					this.P[i][j] = rngB.nextDouble(this.tau0 + this.A[i][j] + this.A[j][i], 
							this.tau1 + 2 - this.A[i][j] - this.A[j][i]);
				}
			}			
		}
	}


	/*
	 * function to sample A
	 */
	public void Sample_A(int N, ArrayList<Event> data, double tmax, int N_bin, 
			int binmax, Baseline_model baseline, Network_model network,
			double dt, Impulse_model impulse,
			Random rand, int[][] events, boolean skip_reestimation){
		if(skip_reestimation){
			return;
		}
		// starting likelihood
		double start_lik = 0.0;
		for(int k = 0; k < this.K; k++){
			double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, k, dt, this.A);		
			start_lik += Evaluation.Get_loglik_poisson_single(currentRate, events[k], dt);
		}
		// keep track of number of acceptance
		int change = 0;
		int nochange = 0;
		// starting loop
		for(int kp = 0; kp < this.K; kp++){

			
			// than try change every sender
			for(int k = 0; k < this.K; k++){
				// if not self connected
				if((!this.self) & k == kp){
					continue;
				}
				// keep track of current link
				int current_link = this.A[k][kp];
				//double[] newRate = this.test_rate_change(N, data, tmax, N_bin, impulse, currentRate, k, kp, dt);
				// calculate likelihood for a receiver process 
				double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
						data, baseline, network, impulse, kp, dt, this.A);
				double current_lik = Evaluation.Get_loglik_poisson_single(currentRate, events[kp], dt);

				// calculate again with flipped A 
				this.A[k][kp] = 1 - this.A[k][kp];
				double[] newRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
						data, baseline, network, impulse, kp, dt, this.A);

				double lik_new = Evaluation.Get_loglik_poisson_single(newRate, events[kp], dt);

				// get likelihood
				double lik0 = 0, lik1 = 0;
				if(current_link == 1){
					lik0 = lik_new; lik1 = current_lik;
				}else{
					lik0 = current_lik; lik1 = lik_new;
				}

				// add prior information
				lik0 += Math.log(1 - this.P[k][kp]);
				lik1 += Math.log(this.P[k][kp]);

				// calculate test statistics
				//				double test = lik1 - Math.log(Math.exp(lik0) + Math.exp(lik1));
				double test = lik1 - MathUtil.logSumOfExponentials(lik0, lik1);
				double u = Math.log(rand.nextDouble());
 				if(test >  u){
					this.A[k][kp] = 1;
				}else{
					this.A[k][kp] = 0;
				}
				if(this.A[k][kp] == current_link){
					nochange ++;
				}else{
					change ++;
				}
			}
		}
		// calculate final likelihood
		double end_lik = 0.0;
		for(int k = 0; k < this.K; k++){
			double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, k, dt, this.A);		
			end_lik += Evaluation.Get_loglik_poisson_single(currentRate, events[k], dt);
		}
		// print diagnostics
		System.out.println("::::: A update finish :::::");
		System.out.println("Start:  " + start_lik + "  End:  " + end_lik);
		System.out.println("Change: " + (end_lik - start_lik));
		System.out.println("Accept rate: " + (change/(nochange + change + 0.0)));
	}

	/*
	 * function to sample A
	 */
	public void Sample_A_block(int N, ArrayList<Event> data, double tmax, int N_bin, 
			int binmax, Baseline_model baseline, Network_model network,
			double dt, Impulse_model impulse,
			Random rand, int[][] events, boolean skip_reestimation, int how_many_flip){
		if(skip_reestimation){
			return;
		}
		// starting likelihood
		double start_lik = 0.0;
		for(int k = 0; k < this.K; k++){
			double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, k, dt, this.A);		
			start_lik += Evaluation.Get_loglik_poisson_single(currentRate, events[k], dt);
		}
		// keep track of number of acceptance
		int change = 0;
		int nochange = 0;
		how_many_flip = how_many_flip > this.K - 1 ? this.K - 1 : how_many_flip;
		// starting loop for proposing new set of A
		for(int kp = 0; kp < this.K; kp++){

			// keep current A column
			int[] current_column = new int[this.K];

			// calculate likelihood for a receiver process ahead of time before changing
			double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, kp, dt, this.A);
			double current_lik = Evaluation.Get_loglik_poisson_single(currentRate, events[kp], dt);
			// save current column
			for(int k = 0; k < this.K; k++){
				current_column[k] = this.A[k][kp];
			}
			
			// propose new column
//			for(int k = 0; k < this.K; k++){
//				// if not self-connected
//				if((!this.self) & k == kp){
//					continue;
//				}else{
//					this.A[k][kp] = (rand.nextDouble() > this.P[k][kp] ? 0 : 1);							
//				}
//			}
			//propose new column again
			int toflip = -1;
			int count_flip = 0;
			while(count_flip < how_many_flip){
				toflip = rand.nextInt(this.K);
				if(toflip != kp){
					this.A[toflip][kp] = 1 - this.A[toflip][kp];
					count_flip ++;
				}
			}
		
			// calculate again with new set of A
			double[] newRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, kp, dt, this.A);				
			double lik_new = Evaluation.Get_loglik_poisson_single(newRate, events[kp], dt);

			// add prior likelihood

			for(int k = 0; k < this.K; k++){
				lik_new += (this.A[k][kp] == 1 
						? Math.log(this.P[k][kp]): Math.log(1 - this.P[k][kp])) ;
				current_lik += (current_column[k] == 1 
						? Math.log(this.P[k][kp]): Math.log(1 - this.P[k][kp])) ;
			}
			

			// calculate test statistics
			//				double test = lik1 - Math.log(Math.exp(lik0) + Math.exp(lik1));
			double test = lik_new - MathUtil.logSumOfExponentials(current_lik, lik_new);
			double u = Math.log(rand.nextDouble());
			if(test >  u){
				// do nothing, accept this A column
				change ++;
			}else{
				// move back original column
				for(int k = 0; k < this.K; k++){
					this.A[k][kp] = current_column[k];
				}
				nochange ++; 
			}
		}
		// calculate final likelihood
		double end_lik = 0.0;
		for(int k = 0; k < this.K; k++){
			double[] currentRate = Evaluation.getAggRate_onechain_updatedAdj(N, this.K, N_bin, binmax, tmax, 
					data, baseline, network, impulse, k, dt, this.A);		
			end_lik += Evaluation.Get_loglik_poisson_single(currentRate, events[k], dt);
		}
		// print diagnostics
		System.out.println("::::: A update finish :::::");
		System.out.println("Start:  " + start_lik + "  End:  " + end_lik);
		System.out.println("Change: " + (end_lik - start_lik));
		System.out.println("Accept rate: " + (change/(nochange + change + 0.0)));
	}


	/*
	 * function to test event rate change
	 */
	public double[] test_rate_change(int N, ArrayList<Event> data, double tmax, double N_bin,
			Impulse_model impulse, 
			double[] rate, int k, int kp, double dt){
		//		// get the rates for each time point on kp
		//		double[] rate = allrates[kp].clone();
		// see if it is there already
		boolean exist = (this.A[k][kp] == 1);
		int sign;
		if(exist){ sign = -1; }else{  sign = 1;  }

		// change all events happening on k influence on kp
		for(int i = 0; i < N; i++){
			if(data.get(i).proc != k) continue;

			// for all following time points
			for(int l=1; (l * dt) < tmax; l++){
				// if the bin count plus influence is larger than total bin
				int impacted_bin = data.get(i).timebin + l; 
				if(impacted_bin >= N_bin) break;
				double impulse_to_add = impulse.calculateImpulse( (l+0.0)*dt, k, kp, tmax);
				rate[impacted_bin] += impulse_to_add * sign;
				if (Double.isNaN(rate[impacted_bin]) | rate[impacted_bin] < 0) {
					System.out.println("Error in bin " + impacted_bin);
				}
			}
		}
		return(rate);
	}

	public void diag_config(){
		System.out.println("::::network diag::::");
		System.out.println("no. of events on process k               : " + Arrays.toString(this.Nk));
		System.out.println("no. of baseline events from process k    : " + Arrays.toString(this.N0k));
//		double[] parenting = new double[this.K];
//		for(int i = 0; i < this.K; i++) parenting[i] = VecUtil.getSum(this.Nkk[i]);
//		System.out.println("no. of events caused by process k        :" + Arrays.toString(parenting));
//		double[] selfparenting = new double[this.K];
//		for(int i = 0; i < this.K; i++) selfparenting[i] = this.Nkk[i][i];
//		System.out.println("no. of events caused by process k itself :" + Arrays.toString(selfparenting));
		System.out.println();
	}

	public void diag_w(){
		System.out.println("::::W resampling diag::::");
		System.out.println("beta_w     : " + this.beta_w);
		//		System.out.println("all A      : " + Arrays.deepToString(this.A));
		//		System.out.println("all W      : " + Arrays.deepToString(this.W));
		double[][] effW = new double[this.K][this.K];
		for(int i = 0; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				effW[i][j] = this.W[i][j] * this.A[i][j];
			}
		}
		//		System.out.println("effective W: " + Arrays.deepToString(effW));
		System.out.println("W sum      : " + VecUtil.sum2d(this.W));
		System.out.println("A sum      : " + VecUtil.sum2d(this.A));
		System.out.println("W eff sum  : " + VecUtil.sum2d(effW));

		System.out.println();
	}

	public void diag_p(){
		//		System.out.println("::::P resampling diag::::");
		//		System.out.println("all P      : " + Arrays.deepToString(this.P));
	}
	public void diag_A(){
		System.out.println("::::A resampling diag::::");
		//		System.out.println("all A      : " + Arrays.deepToString(this.A));
		int[] give = new int[this.K];
		int[] receive = new int[this.K];
		for(int i = 0; i < this.K; i++){
			for(int j = 0; j < this.K; j++){
				if(this.A[i][j] == 1){
					give[i] ++;
					receive[j] ++;
				}
			}
		}
		System.out.println("Potential child    : " + Arrays.toString(give));
		System.out.println("Potential parent     : " + Arrays.toString(receive));
	}
}
