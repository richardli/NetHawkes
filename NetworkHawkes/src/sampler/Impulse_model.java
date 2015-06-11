package sampler;

import java.util.ArrayList;
import java.util.Arrays;

import mathutil.VecUtil;
import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.Normal;

public class Impulse_model {
	int K;
	double kap0;
	double mu0; 
	double alpha0;
	double beta0;
	// mu matrix
	public double[][] mu = new double[K][K];
	// tau matrix
	public double[][] tau = new double[K][K];
	// sufficient statistics
	private double[][] xbar = new double[K][K];
	private double[][] s2 = new double[K][K];
	private double[][] m = new double[K][K];
	
	public Impulse_model(int K){
		this.K = K;
	}

	public void pass_hyper(double alpha0, double beta0, double kap0, double mu0){
		this.alpha0 = alpha0;
		this.beta0 = beta0; 
		this.kap0 = kap0; 
		this.mu0 = mu0;
		this.mu = new double[K][K];
		this.tau = new double[K][K];
		this.xbar = new double[K][K];
		this.s2 = new double[K][K];
		this.m = new double[K][K];
	}

	/*
	 * function to calculate impulse function (k: parent, kp: child)
	 * note t and tmax are actual time, not bin number!!
	 */
	public double calculateImpulse(double t, int k, int kp, double tmax){
		double g = tmax / (t * (tmax - t)) * Math.sqrt(this.tau[k][kp] / 2 / Math.PI);
		g *= Math.exp(-this.tau[k][kp]/2 * Math.pow((Math.log(t/(tmax-t)) - this.mu[k][kp]), 2));
		return(g);
	}


	/*
	 * extract k-k' pair index (k0: parent, k1: event)
	 */
	public ArrayList<Integer> getKKP(ArrayList<Event> data, 
			int k0, int k1){
		ArrayList<Integer> list = new ArrayList<Integer>();
		for(int i = 0; i < data.size(); i++){
			if(data.get(i).proc == k1 && 
					data.get(i).pc == k0) list.add(i);
		}
		return(list);
	}

	/*
	 * function to sample theta
	 */
	public void Sample_theta(ArrayList<Event> data, 
			Gamma rngG, Normal rngN, boolean dont_update){
		
		if(dont_update){
			for(int i = 0;  i < this.K; i++){
				for(int j = 0; j < this.K; j++){
					// update step
//					this.tau[i][j] = rngG.nextDouble(alpha0, beta0);
//					this.mu[i][j] = rngN.nextDouble(mu0, 1/(kap0 * this.tau[i][j]));
					
					// mean update
					this.tau[i][j] = alpha0 / beta0;
					this.mu[i][j] = mu0;
				}
			}
			return;
		}
		
		for(int i = 0;  i < this.K; i++){
			for(int j = 0; j < this.K; j++){

				// calculate sufficient statistics
				ArrayList<Integer> inter_events = getKKP(data, i,  j);
				if(inter_events.size() == 0){
					// update step
					this.tau[i][j] = rngG.nextDouble(this.alpha0, this.beta0);
					this.mu[i][j] = rngN.nextDouble(this.mu0, 1/(this.kap0 * this.tau[i][j]));
				}else{
					double[] logit_inter_times = new double[inter_events.size()];
					for(int ii = 0; ii < inter_events.size(); ii++){
						logit_inter_times[ii] = data.get(inter_events.get(ii)).logit_deltat;
					}
					xbar[i][j] = VecUtil.getMean(logit_inter_times);
					s2[i][j] = VecUtil.getVar(logit_inter_times);
					m[i][j] = logit_inter_times.length;

					// calculate new paramters
					double mu_new, kap_new, alpha_new, beta_new;
					mu_new = (kap0 * mu0 + m[i][j] * xbar[i][j]) / (kap0 + m[i][j]);
					kap_new = kap0 + m[i][j];
					alpha_new = alpha0 + 1/2;
					beta_new = beta0 + 1/2 * (m[i][j] * s2[i][j] + (kap0 * m[i][j] * (xbar[i][j] - mu0)*(xbar[i][j] - mu0))/(kap0 + m[i][j]));
					// update step
					this.tau[i][j] = rngG.nextDouble(alpha_new, beta_new);
					this.mu[i][j] = rngN.nextDouble(mu_new, 1/(kap_new * this.tau[i][j]));
					
				}
			}
		}
	}
	
	public void diag_theta(){
		System.out.println("::::impulse diag::::");
		System.out.println("xbar   : " + Arrays.deepToString(this.xbar));
		System.out.println("s^2    : " + Arrays.deepToString(this.s2));
		System.out.println("count  : " + Arrays.deepToString(this.m));		
		System.out.println("mu     : " + Arrays.deepToString(this.mu));
		System.out.println("tau    : " + Arrays.deepToString(this.tau));
		System.out.println();
	}

}
