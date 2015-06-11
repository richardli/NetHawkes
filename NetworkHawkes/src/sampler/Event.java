package sampler;

import java.util.HashSet;

public class Event {
	public double time;
	public int proc;
	public int timebin;
	
	public int parent;
	public int pc;
	public double deltat; 		// time since parent event
	public double logit_deltat;

	public HashSet<Integer> plist; 
	
	Event(double time){
		this.time = time;
		this.plist = new HashSet<Integer>();
		this.updateParent();
	}
	
//	Event(double time, int proc){
//		this.time = time;
//		this.proc = proc;
//		this.parent = -1;  // -----> background rate
//		this.plist =  new HashSet<Integer>();
//	}
	
	public void setProc(int proc){
		this.proc = proc;
	}
	
	public void addPotentialParent(int parent){
		this.plist.add(parent);
	}
	
	// update parent using the parent index and parent itself
	public void updateParent(int i, Event parent_event, double tmax){
		this.parent = i;
		this.pc = parent_event.proc;
		this.deltat = this.time - parent_event.time;
		this.logit_deltat = calculate_logit_delta(this.deltat, tmax);
	}
	
	// update parent without any information --> set to background
	public void updateParent(){
		this.parent = -1;
		this.pc = -1;
		this.deltat = 0.0;
		this.logit_deltat = 0.0;
	}
	// update parent to be the same to another event
	public void copyParent(Event another){
		this.parent = another.parent;
		this.pc = another.pc; 
		this.deltat = another.deltat;
		this.logit_deltat = another.logit_deltat;
	}
//	public void updateParent(int parent, int pc, double deltat, double tmax){
//		this.parent = parent;
//		this.pc = pc;
//		this.deltat = deltat;
//		this.logit_deltat = calculate_logit_delta(deltat, tmax);
//	}
	
	public void binning(double bin_width, double T){
		// stay with the original, avoid time = T
		this.timebin = (int)(this.time / bin_width);
		if(this.timebin > 0 & (this.time % bin_width == 0)) this.timebin--;
	}
	
	public void binning(double bin_width, double T, double offset){
		// stay with the original, avoid time = T
		this.timebin = (int)((this.time - offset) / bin_width);
		if(this.timebin > 0 & ((this.time - offset) % bin_width == 0)) this.timebin--;
//		if(this.timebin == T){
//			System.out.println("oops");
//		}

	}
	
	public double calculate_logit_delta(double t, double tmax){
		return( Math.log(t / (tmax - t) ));
	}
}
