/**@author Tyler Peryea
 * 
 * This is a really naive and simple program to calculate the synergy scores
 * used in the AZ drug combination dream challenge.
 * 
 * Given matrix response values, and hill parameters for monotherapies,
 * it will reproduce the provided synergy scores to <1E-5 accuracy
 * for all but 21 non-replicate values. 
 * 
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class SynergyCalculator {
	public static class HillParams {
		double hillSlope;
		double ic50;
		double emax;
		double emin;
		double maxC;

		public HillParams(double ic, double s, double emax, double emin) {
			this.hillSlope = s;
			this.ic50 = ic;
			this.emax = emax;
			this.emin = emin;
		}

		public double predictPartial(double x, double y) {
			double ef1 = (y - emin) / (emax - y);
			if (ef1 < 0)
				return 0;
			double p1 = x
					/ (ic50 * Math.signum(ef1) * Math.pow(Math.abs(ef1),
							1 / hillSlope));
			return p1;
		}

		public double sample(double dose) {
			double resp = emin
					+ (emax - emin)
					/ (1 + Math.pow(10, (Math.log10(ic50) - Math.log10(dose))
							* hillSlope));
			return resp;
		}
	}

	public static class Experiment {
		public HillParams hpA;
		public HillParams hpB;
		public String cmpA;
		public String cmpB;
		public String cell;
		public double syn;
		public Object qa;

		public String getID() {
			return cell + "" + cmpA + "." + cmpB;
		}

		public double[][] predictLoewes(double[] c1, double[] c2) {
			double[][] predict = new double[c1.length][c2.length];
			for (int i = 0; i < c1.length; i++) {
				for (int j = 0; j < c2.length; j++) {
					double t = findBest(c1[i], c2[j], hpA, hpB);
					predict[i][j] = t;
				}
			}
			return predict;
		}

		public double getSynergy(double[] c1, double[] c2, double[][] matrix) {
			double[][] lm = predictLoewes(c1, c2);

			for (int i = 0; i < lm.length; i++) {
				for (int j = 0; j < lm.length; j++) {
					matrix[i][j] = Math.max(Math.min(100 - matrix[i][j], 200),
							-100);
				}
			}
			return getDiffSynergy(c1, c2, matrix, lm);
		}

		public static double getDiffSynergy(double[] c1, double[] c2,
				double[][] real, double[][] predicted) {
			double sum = 0;
			double sumArea = 0;

			for (int i = 0; i < predicted.length; i++) {
				for (int j = 0; j < predicted[0].length; j++) {
					if (i > 0 && j > 0 && j < 5 && i < 5) {
						double dx = Math.log(c1[i + 1]) - Math.log(c1[i]);
						double dy = Math.log(c2[j + 1]) - Math.log(c2[j]);

						double dval1 = ((real[i][j] - predicted[i][j]));
						double dval2 = ((real[i + 1][j] - predicted[i + 1][j]));
						double dval3 = ((real[i][j + 1] - predicted[i][j + 1]));
						double dval4 = ((real[i + 1][j + 1] - predicted[i + 1][j + 1]));
						double dsum = (dval1 + dval2 + dval3 + dval4);

						// this should happen for correct approx,
						// but does not happen in provided values
						// uncomment for more accuracy
						//
						// dsum=dsum/4;

						sum += dsum * dx * dy;
						sumArea += dx * dy;

					}
				}
			}

			// This is wrong, as it turns out some
			// areas are larger than others
			// However, it does make the values agree with those
			// provided.

			// comment for better accuracy
			sumArea = 21.20759244;

			return sum / sumArea;
		}
	}

	

	public static void readRaw(String fname, double[] c1, double[] c2,
			double[][] matrix) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(fname)));
		String line = null;
		int i = 0;
		while ((line = br.readLine()) != null) {
			String[] cols = line.replace("\"", "").split(",");
			if (i == 0) {
				Arrays.toString(parseStringArray(cols, c2, 1, 6));
			}
			if (i > 0 && i <= 6) {
				c1[i - 1] = Double.parseDouble(cols[0]);
				parseStringArray(cols, matrix[i - 1], 1, 6);
			}
			i++;
		}
		br.close();
	}

	public static double[] parseStringArray(String[] s, double[] out, int off,
			int max) {
		if (out == null)
			out = new double[Math.min(s.length - off, max)];
		for (int i = off; i < Math.min(s.length, max + off); i++) {
			double d = Double.parseDouble(s[i]);
			out[i - off] = d;
		}
		return out;
	}

	public static double findBest(double x1, double x2, HillParams hpA,
			HillParams hpB) {
		double guess = findBestInterp(x1, x2, hpA, hpB);
		return guess;
	}

	//Search for the root by linear interpolation
	public static double findBestInterp(double x1, double x2, HillParams hpA,
			HillParams hpB) {

		double max = 100;
		double min = 0;
		for (int i = 0; i < 100; i++) {
			double guess = (max + min) / 2;
			double e = evaluate(x1, x2, hpA, hpB, guess);
			if (e > 0) {
				max = guess;
			} else if (e < 0) {
				min = guess;
			} else {
				break;
			}
		}

		return (max + min) / 2;
	}

	public static double evaluate(double x1, double x2, HillParams hpA,
			HillParams hpB, double y) {
		double p1 = hpA.predictPartial(x1, y);
		double p2 = hpB.predictPartial(x2, y);
		double err = 1 - (p1 + p2);
		return err;
	}
	
	public static void main(String[] args) throws IOException {
		Map<String, Experiment> begin = new HashMap<String, Experiment>();

		if(args.length<1){
			System.out.println("Please provide path to dream folder");
		}
		
		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(args[0]
						+ "/data/ch1_train_combination_and_monoTherapy.csv")));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {

			String[] cols = line.replace("\"", "").split(",");
			Experiment exp = new Experiment();
			exp.cell = cols[0];
			exp.cmpA = cols[1];
			exp.cmpB = cols[2];
			exp.syn = Double.parseDouble(cols[11]);
			exp.qa = cols[12];
			HillParams hpA = new HillParams(Double.parseDouble(cols[5]),
					Double.parseDouble(cols[6]),
					100 - Double.parseDouble(cols[7]), 0);
			hpA.maxC = Double.parseDouble(cols[3]);
			HillParams hpB = new HillParams(Double.parseDouble(cols[8]),
					Double.parseDouble(cols[9]),
					100 - Double.parseDouble(cols[10]), 0);
			hpB.maxC = Double.parseDouble(cols[4]);
			exp.hpA = hpA;
			exp.hpB = hpB;
			begin.put(exp.getID(), exp);

		}
		br.close();

		File folder = new File(args[0]
				+ "/data/Raw_Data_csv/ch1_training_combinations");
		File[] listOfFiles = folder.listFiles();
		Map<String, String> output = new HashMap<String, String>();

		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile()) {
				String fname = listOfFiles[i].getName();
				String[] nameP = fname.split("\\.");
				String cmpA = nameP[0];
				String cmpB = nameP[1];
				String cell = nameP[2];

				Experiment exp = begin.get(cell + "" + cmpA + "." + cmpB);
				double[] c1 = new double[6];
				double[] c2 = new double[6];
				double[][] matrix = new double[6][6];

				readRaw(listOfFiles[i].getAbsolutePath(), c1, c2, matrix);
				if (exp != null && exp.qa.equals("1")) {
					double[] r2 = matrix[0];
					double[] r1 = new double[r2.length];
					for (int k = 0; k < r1.length; k++) {
						r1[k] = matrix[k][0];
					}

					//Ignore those with >1 replicate for now
					if (output.containsKey(exp.getID())) {
						output.put(exp.getID(), "NULL");
					} else {
						try {
							double calc = exp.getSynergy(c1, c2, matrix);
							output.put(exp.getID(), exp.getID() + "\t" + calc
									+ "\t" + exp.syn + "\t" + fname + "\t"
									+ (calc - exp.syn));
						} catch (Exception e) {

						}
					}
				}
			}
		}
		
		for (String key : output.keySet()) {
			String o = output.get(key);
			if (!o.equals("NULL"))
				System.out.println(output.get(key));
		}

	}
}
