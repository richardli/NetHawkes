package sampler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import mathutil.MathUtil;
import mathutil.VecUtil;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.Normal;

public class Baseline_model {
	int K;
	double T;
	int N_bin;
	double alpha_lam;
	double beta_lam;
	public double[][] lambda0 = new double[K][N_bin];
	// for LGCP case
	public int M;
	public int ncell;
	public double dt_LGCP;
	public int[][] events_dense;
	public double[][] rate_dense;
	public double[][] y0 = new double[K][M];

	
	public double[] alpha = new double[K];
	public double[] mu = new double[K];
	public double[][] cov = new double[M][M];
	public MultivariateNormalDistribution priorDist;
	double Period;
	
	// hyper parameters for alpha and mu
	// assume alpha ~ LN(ap, as), mu ~ LN(mp, ms)
	public double[] ap = new double[K];
	public double[] mp = new double[K];
	public double[] as = new double[K];
	public double[] ms = new double[K];
	
	public Baseline_model(int K, double T, int N_bin){
		this.K = K;
		this.N_bin = N_bin;
		this.T = T;
		this.lambda0 = new double[K][N_bin];
	}

	public void pass_hyper(double alpha_lam, 
			double beta_lam){
		this.alpha_lam = alpha_lam;
		this.beta_lam = beta_lam;
	}
	
	public void initLGCP(int[][] events, double T, int M, double Period, double sigma_LGCP){
		// calculate number of cells in each block
		this.ncell = (int)(this.N_bin / (M + 0.0) + 0.9999);
		this.M = M;
		this.dt_LGCP = this.T / this.M;
		this.Period = Period;
		
		// get denser version of the events
		this.events_dense = new int[this.K][M];
		for(int i = 0; i < this.N_bin; i++){
			int m = (int)(i / (ncell + 0.0));
			for(int j = 0; j < this.K; j++){
				this.events_dense[j][m] += events[j][i];
			}
		}
		
		// initialize background rate hypers
		this.alpha = new double[K];
		this.mu = new double[K];
		this.ms = new double[K];
		this.as = new double[K];
		this.mp = new double [K];
		this.ap = new double[K];
		this.cov = new double[this.M][this.M];
		this.y0 = new double[this.K][this.M];
		this.rate_dense = new double[this.K][this.M];
		// constants to be used later
		double A = Math.exp(1) - 1;
		double B = Math.exp(2) - 1;
		
		for(int i = 0; i < this.K; i++){
			double max = (double) (MathUtil.max(events_dense[i]));
			double min = (double) (MathUtil.min(events_dense[i]));
			// how to set mean?? actual mean to far away from max count in practice
			double mm = (double) (VecUtil.getMean(events_dense[i]));
			double d = Math.min(mm - min, max - mm);
			// set hyper parameters directly
//			// the following calculated for LN prior, somehow does not work			
//			// calculate expectation for mu and alpha
			double delta = (A*B)*d*d - 4*A*B*mm*mm;
			delta = delta > 0 ? delta : 0;
			double s = ( 2*A*mm + Math.sqrt(delta) ) / (2*A+2);
			double t = mm - s;
			t  = t > 0 ? t : 1e-8;
			// get hyper priors
			this.mp[i] = Math.log(s) - 0.5;
			this.ap[i] = Math.log(t) - 1;
			this.ms[i] = sigma_LGCP;
			this.as[i] =  sigma_LGCP;
			// get first sample
			this.mu[i] = Math.exp(this.mp[i] + 0.5); 
			this.alpha[i] = Math.exp(this.ap[i] + 0.5);
		}
		
		// initialize y0
		this.first_sample_LGCP(this.Period);
	}
	
	public void initLGCP_old(int[][] events, double T, int M, double Period){
		// calculate number of cells in each block
		this.ncell = (int)(this.N_bin / (M + 0.0) + 0.9999);
		this.M = M;
		this.dt_LGCP = this.T / this.M;
		this.Period = Period;
		
		// get denser version of the events
		this.events_dense = new int[this.K][M];
		for(int i = 0; i < this.N_bin; i++){
			int m = (int)(i / (ncell + 0.0));
			for(int j = 0; j < this.K; j++){
				this.events_dense[j][m] += events[j][i];
			}
		}
		
		// initialize background rate hypers
		this.alpha = new double[K];
		this.mu = new double[K];
		this.cov = new double[this.M][this.M];
		this.y0 = new double[this.K][this.M];
		this.rate_dense = new double[this.K][this.M];
		this.ms = new double[K];
		this.as = new double[K];
		this.mp = new double [K];
		this.ap = new double[K];
		
		for(int i = 0; i < this.K; i++){
			double max = (double) (MathUtil.max(events_dense[i]));
			double min = (double) (MathUtil.min(events_dense[i]));
			double mm = (double) (VecUtil.getMean(events_dense[i]));
			// set hyper parameters directly
			double a0 = (-2 * mm + Math.sqrt(4*mm*mm + 3*(max-mm)*(max-mm))) / 6.0;
			double a1 = (-2 * mm + Math.sqrt(4*mm*mm + 3*(mm-min)*(mm-min))) / 6.0;
			this.alpha[i] = a0 > a1 ? a0 : a1;
			this.mu[i] = mm - alpha[i] / 2.0;
			if(this.mu[i] < 0) this.mu[i] = 1E-9;
			
			this.mp[i] = Math.log(this.mu[i]) - 0.5; 
			this.ap[i] = Math.log(this.alpha[i]) - 0.5;
			
		}
		
		// initialize y0
		this.first_sample_LGCP(this.Period);
	}
	
	public void period_kernel(double Period){
		double l = 0.15;
		this.cov = new double[this.M][this.M];
		for(int i = 0; i < this.M; i++){
			for(int j = 0; j < this.M; j++){
				this.cov[i][j] = Math.exp(-2 * Math.pow(Math.sin(Math.abs(i-j)*this.dt_LGCP/(Period+0.0) * Math.PI), 2)/ (l+0.0)) ;
			}
		}
		return;
	}
	
	public void first_sample_LGCP(double Period){
		period_kernel(Period);
		this.priorDist= new MultivariateNormalDistribution(new double[this.M], this.cov);
		for(int i = 0; i < this.K; i++){
			this.y0[i] = priorDist.sample();
			this.rate_dense[i] = y_to_rate(this.y0[i], i);
		}
		interpolate(this.lambda0, this.rate_dense);
	}
	
	// calculate rate using current mu and alpha
	public double[] y_to_rate(double[] sample, int k){
		double[] rate = new double[sample.length];
		for(int i = 0; i < sample.length; i++){
			rate[i] = this.alpha[k] * Math.exp(sample[i]) + this.mu[k];
		}
		return(rate);
	}
	
	// calculate loglik using current mu and alpha	
	public double LGCP_lik(double[] sample,  Network_model network, int k, 
			int N, ArrayList<Event> data){
		double lik = 0;
		// TODO: finish here
		double[] rate = y_to_rate(sample, k);
		
		double[] longrates = new double[this.N_bin];
		interpolate(longrates, rate);
		int[] events_N0 = network.countBaseline(N, data, this.M, this.N_bin, this.ncell, k);
		lik = Evaluation.Get_loglik_poisson_single(longrates, events_N0, 1);
		return(lik);
	}
	
	// calculate rate using new mu and alpha
	public double[] y_to_rate(double[] sample, double mu, double alpha, int k){
		double[] rate = new double[sample.length];
		for(int i = 0; i < sample.length; i++){
			rate[i] = alpha * Math.exp(sample[i]) + mu;
		}
		return(rate);
	}
	// calculate loglik using new mu and alpha
	public double LGCP_lik(double[] sample, double mu, double alpha,  
			Network_model network, int k, int N, ArrayList<Event> data){
		double lik = 0;
		// TODO: finish here
		double[] rate = y_to_rate(sample, mu, alpha, k);
		double[] longrates = new double[this.N_bin];
		
		interpolate(longrates, rate);
		int[] events_N0 = network.countBaseline(N, data, this.M, this.N_bin, this.ncell, k);
		lik = Evaluation.Get_loglik_poisson_single(longrates, events_N0, 1);
		return(lik);
	}
	
	
	public void sampleLGCP(Network_model network, Random rand, Normal rngN, 
			int N, ArrayList<Event> data){
		double[] prior_sample = this.priorDist.sample();
//		System.out.println(Arrays.toString(prior_sample));
		for(int i = 0; i < this.K; i++){
			// sampler from prior 0-mean normal 
			double prior_sample_centered_lnmu =  rngN.nextDouble() * this.ms[i];
			double prior_sample_centered_lnalpha = rngN.nextDouble() * this.as[i];
//			if(i == 0){
//				System.out.println(prior_sample_centered_lnmu);
//				System.out.println(prior_sample_centered_lnalpha);
//			}
			// perform elliptical slice sampling here
			// input: current y(t) or rate, prior sample of y(t) or rate
			// output: new y(t) or rate
		
			double current_lik = LGCP_lik(this.y0[i], network, i, N, data);
			
			double test = Math.log(rand.nextDouble()) + current_lik;
//			double test = current_lik;
			// set up bracket
			double phi = rand.nextDouble() * 2 * Math.PI;
			double phi_min = phi - 2 * Math.PI;
			double phi_max = phi;
			double[] ynew = new double[this.M];
			double centered_lnmu_new =  0;
			double centered_lnalpha_new = 0;
			double mu_new;
			double alpha_new;

			double new_lik = 0;
			while(true){
				double max_dist = 0;
				for(int j = 0; j < this.M; j++){
					ynew[j] = this.y0[i][j] * Math.cos(phi) + prior_sample[j] * Math.sin(phi);				
					double diff = Math.abs(ynew[j] - this.y0[i][j]);
					if(diff > max_dist) max_dist = diff;
				}
				centered_lnmu_new = (Math.log(this.mu[i]) - this.mp[i])* Math.cos(phi) + prior_sample_centered_lnmu* Math.sin(phi);
//				mu_new =  Math.exp(this.mp[i] + (Math.log(this.mu[i]) - this.mp[i]) * Math.cos(phi) 
//						+ (Math.log(prior_sample_mu) - this.mp[i]) * Math.sin(phi));
				double diff = Math.abs(centered_lnmu_new - prior_sample_centered_lnmu);
				if(diff > max_dist) max_dist = diff;
				
				centered_lnalpha_new = (Math.log(this.alpha[i]) - this.ap[i])* Math.cos(phi) + prior_sample_centered_lnalpha* Math.sin(phi);
//				alpha_new = Math.exp(this.ap[i] + (Math.log(this.alpha[i]) - this.ap[i]) * Math.cos(phi) 
//						+ (Math.log(prior_sample_alpha) - this.ap[i]) * Math.sin(phi));
				diff = Math.abs(centered_lnalpha_new - prior_sample_centered_lnalpha);
				if(diff > max_dist) max_dist = diff;
				
				// get mu and alpha
				mu_new = Math.exp(centered_lnmu_new + this.mp[i]);
				alpha_new = Math.exp(centered_lnalpha_new + this.ap[i]); 
//				if(i == 0){
//					System.out.println(mu_new);
//					System.out.println(alpha_new);
//				}
				new_lik = LGCP_lik(ynew, mu_new, alpha_new, network, i, N, data);
				if(new_lik >= test){
					break;
				}else if(max_dist < 1e-9){
					break;
				}else{
					if(phi > 0) phi_max = phi;
					if(phi < 0) phi_min = phi;
					phi = rand.nextDouble() * (phi_max - phi_min) + phi_min;
				}
			}		

			if(new_lik >= test){
				// calculate rates
				double[] ratenew = y_to_rate(ynew, mu_new, alpha_new, i);
				// update y0
				this.mu[i] = mu_new;
				this.alpha[i] = alpha_new;
				for(int j = 0; j < this.M; j++){
					this.y0[i][j] = ynew[j];
					this.rate_dense[i][j] = ratenew[j];
				}
			}	
		}
		// interpolate
		interpolate(this.lambda0, this.rate_dense);
	}
	
	public void sampleLGCP_fixed(Network_model network, Random rand, 
			int N, ArrayList<Event> data){
		for(int i = 0; i < this.K; i++){
			double[] prior_sample = this.priorDist.sample();
			// perform elliptical slice sampling here
			// input: current y(t) or rate, prior sample of y(t) or rate
			// output: new y(t) or rate
		
			double current_lik = LGCP_lik(this.y0[i], network, i, N, data);
			
			double test = Math.log(rand.nextDouble()) + current_lik;
			// set up bracket
			double phi = rand.nextDouble() * 2 * Math.PI;
			double phi_min = phi - 2 * Math.PI;
			double phi_max = phi;
			double[] ynew = new double[this.M];
			
			while(true){
				double max_dist = 0;
				for(int j = 0; j < this.M; j++){
					ynew[j] = this.y0[i][j] * Math.cos(phi) + prior_sample[j] * Math.sin(phi);				
					double diff = Math.abs(ynew[j] - this.y0[i][j]);
					if(diff > max_dist) max_dist = diff;
				}
				
				double new_lik = LGCP_lik(ynew, network, i, N, data);
				if(new_lik >= test){
					break;
				}else if(max_dist < 1e-9){
					break;
				}else{
					if(phi > 0) phi_max = phi;
					if(phi < 0) phi_min = phi;
					phi = rand.nextDouble() * (phi_max - phi_min) + phi_min;
				}
			}		

			// calculate rates
			double[] ratenew = y_to_rate(ynew, i);

			// update y0
			for(int j = 0; j < this.M; j++){
				this.y0[i][j] = ynew[j];
				this.rate_dense[i][j] = ratenew[j];
			}
		}
		// interpolate
		interpolate(this.lambda0, this.rate_dense);
	}
	
	public void interpolate(double[][] rates, double[][] rate_dense){
		// update lambda0 from rate_dense
		for(int k = 0; k < this.K; k++){
			interpolate(rates[k], rate_dense[k]);
		}
	}
	
	public void interpolate(double[] rates, double[] rate_dense){
		for(int i = 0; i < this.N_bin; i++){
			int below = (int)(i % this.ncell);
			int index = (int)(i / this.ncell);
			// first cell
			if(index == 0 & below < this.ncell/2){
				rates[i] = rate_dense[0] - 
							(rate_dense[0] - rate_dense[1]) * (0.5 - below/(this.ncell+0.0));
				if(rates[i] < 0) rates[i] = 1e-9;
			// last cell
			}else if(index == this.M - 1 & below >= this.ncell/2){
				rates[i] = rate_dense[this.M - 1] + 
						(rate_dense[this.M - 1] - rate_dense[this.M - 2]) * (below/(this.ncell+0.0) - 0.5);
				if(rates[i] < 0) rates[i] = 1e-9;
			// in the middle first half cell
			}else if(below < this.ncell / 2){
				rates[i] = rate_dense[index] -
						(rate_dense[index] - rate_dense[index-1]) * (0.5 - below/(this.ncell+0.0));
			// in the middle second half cell	
			}else if(below >= this.ncell / 2){
				rates[i] = rate_dense[index] -
						(rate_dense[index] - rate_dense[index+1]) * (below/(this.ncell+0.0) - 0.5);
			}
			
			// take into account of splitting cell
			if(rates[i] < 0){
				System.out.println("Negative!");
			}
			rates[i] /= this.ncell;
		}
	}
	
	public void initconst(int[][] events, double T){
		for(int i = 0; i < this.K; i++){
			int total = VecUtil.getSum(events[i]);
			rate_dense[i] = VecUtil.fill(this.lambda0[i], (total + 0.0) / T);
		}
	}
	
	public void initconst(int[][] events, double T, double const_rate){
		if(const_rate < 0){
			this.initconst(events, T);
		}
		for(int i = 0; i < this.K; i++){
			this.lambda0[i] = VecUtil.fill(this.lambda0[i], const_rate);
		}
	}
	
	public void initconst(int[][] events, double T, Gamma rngG){
		for(int i = 0; i < this.K; i++){
			double const_rate = rngG.nextDouble(this.alpha_lam, this.beta_lam);
			this.lambda0[i] = VecUtil.fill(this.lambda0[i], const_rate);
		}
	}
	/*
	 * function to sample lambda0, constant case
	 */
	public void Sample_lambda_const(Gamma rngG, Network_model network, boolean skip_reestimation){
		if(skip_reestimation) return;
		for(int i =0; i < this.K; i++){
			double newlambda = rngG.nextDouble(this.alpha_lam + network.N0k[i], this.beta_lam + this.T);
			this.lambda0[i] = VecUtil.fill(this.lambda0[i], newlambda);
		}
	}
	
	public void lambda_const_diag(){
		System.out.println("::::lambda diag::::");
		double[] lambda0_const = new double[this.K];
		for(int i = 0; i < this.K; i++) lambda0_const[i] = this.lambda0[i][0];
		System.out.println("lambda0     : " + Arrays.toString(lambda0_const));
		System.out.println();
	}
	
	public void lambda_LGCP_diag(){
		System.out.println("::::lambda diag::::");
		for(int i = 0; i < this.K; i++) {
			System.out.println("lambda0 process " + i + "    : " 
					+ Arrays.toString(this.lambda0[i]));
		}
		System.out.println();
	}

}
