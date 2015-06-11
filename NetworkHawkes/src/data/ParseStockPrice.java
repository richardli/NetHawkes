package data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

public class ParseStockPrice {
	public static void main(String[] args) throws IOException, ParseException {
		String header = "/Users/zehangli/Dropbox/UW_grad/518-15spr/ImplicitNetwork/NetworkHawkes/";
		String[] spfile = {header + "data/sp100/2015-04-27-5s.txt", 
				header + "data/sp100/2015-04-28-5s.txt", 
				header + "data/sp100/2015-04-29-5s.txt", 
				header + "data/sp100/2015-04-30-5s.txt", 
				header + "data/sp100/2015-05-01-5s.txt"}; 
		String timefile = header + "data/sp100/log5s-w2.txt";
		SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd|HH:mm:ss");
		
		/*
		 *  Output file format:
		 *  Name, price(t0), price(t1), ....
		 */
		HashMap<String, ArrayList<Double>> price = new HashMap<String, ArrayList<Double>>();
		ArrayList<String> time_points = new  ArrayList<String>();
		/*
		 * Read prices
		 */
				for(int i = 0; i < spfile.length; i++){
					BufferedReader sc = new BufferedReader(new FileReader(spfile[i])); 
					String line;
					int count = 0;
					while((line = sc.readLine())!= null){
						count++;
						String[] field = line.split(" ");
						// remove quotation mark
						String stock_name = field[0].substring(1, field[0].length()-1);
						// Date time_saved = format.parse(field[1] + "|" + field[2]);
						// System.out.println(time_saved);
						double lastprice = Double.parseDouble(field[5]);
						if(price.get(stock_name) != null){
							price.get(stock_name).add(lastprice);
						}else{
							price.put(stock_name, new ArrayList<Double>());
							price.get(stock_name).add(lastprice);
						}
					}
					sc.close();
					System.out.println(count + " lines read.");
				}
				// System.out.println(price.get("MDLZ").toString().replace("[", "").replace("]", ""));
				for(String name : price.keySet()){
					BufferedWriter sc = new BufferedWriter(new FileWriter(header + "/data/sp100/" + name + ".txt")); 
					sc.write(price.get(name).toString().replace("[", "").replace("]", ""));
					sc.close();
				}

		/*
		 * Read Times
		 * 
		 */
		BufferedReader sc = new BufferedReader(new FileReader(timefile));
		String line;
		int count = 0;
		while((line = sc.readLine())!=null){
			count ++;
			String[] field = line.split(" ");
			String date_read = field[1].substring(1, field[1].length());
			String time_read = field[2];
			time_points.add(date_read + "|" + time_read);
			//			System.out.println(date_read + "|" + time_read);
			//			System.out.println(format.parse(date_read + "|" + time_read));
		}
		sc.close();

		// System.out.println(price.get("MDLZ").toString().replace("[", "").replace("]", ""));
		BufferedWriter sw = new BufferedWriter(new FileWriter(header + "/data/sp100/time_points.txt")); 
		sw.write(time_points.toString().replace("[", "").replace("]", ""));
		sw.close();
	}
}





