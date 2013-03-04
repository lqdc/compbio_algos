package RabinKarp;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.math.BigInteger;
import java.util.Random;

public class Rk_redo {

	public static long largeModulus = 3;
	public static int pattern_len;
	public static void main(String[] args) {
		File aFile = new File("/home/roman/Desktop/sequence_match.fasta");
		String longText = getContents(aFile);
		String patternToMatch = "GTTACTTGTTTGTAGAATAAGATGTGACCATTC";
		
		long execTimes = 0;// total execution time with each modulus
		long execTime1;// a single execution time
		int runs = 5;// number of runs to get average time for each modulus
		double[] execTime = new double[30];//array with average execution times for each modulus
		long[] modulusArray = new long[30];//array with moduluses that were used
		//Run 30 times. First 15 runs = small modulus. Then from 15 to 30, large random
		long startTime;
		//populate modulus arrays
		for (int i =0; i < 30; i++){
			if (i < 15)
				modulusArray[i] = i+2;
			else
				modulusArray[i] = new BigInteger(27, new Random()).longValue();
		}
		//time the execution
		for (int j = 0; j < 30; j++) {
			execTimes = 0;
			largeModulus = modulusArray[j];
			if (j == 15)
				System.out.println("\nNow using a large modulus\n");
			for (int i = 0; i < runs; i++) {
				execTime1 = 0;
				startTime = System.currentTimeMillis();
				
				int offset = rksearch(patternToMatch, longText);
				if (offset == -1) {
					System.out.println("Pattern not found");
				}
				// else
				// System.out.println("Pattern found! Offset is "+ offset);
				long endTime = System.currentTimeMillis();
				execTime1 = endTime - startTime;
				execTimes += execTime1;
				// System.out.println("Execution took " + (execTime1) +
				// " milliseconds");
			}

			double avgTime = execTimes / (double) runs;
			execTime[j] = avgTime;
			//System.out.println("Average execution time for modulus of "
				//	+ largeModulus + " was " + (avgTime) / 1000.0 + " seconds");
		}
		// graph the execution time averages
		XYChart chart = new XYChart(modulusArray, execTime);
		// estimate average times for small and large moduluses
		double total_large = 0;
		double total_small = 0;
		double avg_large = 0;
		double avg_small = 0;
		for (int i = 0; i < 30; i++)
			if (i < 15)
				total_small += execTime[i];
			else
				total_large += execTime[i];
		avg_large = total_large / 15.0;
		avg_small = total_small / 15.0;
		// print results
		System.out.println("\nAverage for large moduluses = " + avg_large
				+ " milliseconds");
		System.out.println("Average for small moduluses = " + avg_small
				+ " milliseconds");
		double percentage = (avg_small / avg_large - 1) * 100;
		System.out
				.printf("\nSmall moduluses take  %.2f percent longer to find the string "
						+ "than large moduluses in a text of length %d bytes.\n",
						percentage, longText.length());
	}

	public static int rksearch(String pattern, String text){
	    long patternHash, textHash; 
	    int hashesToDo, i; 
	    
	    // get the lengths of the strings and the number of iterations
	    pattern_len = pattern.length();
	    hashesToDo = text.length() - pattern_len ;
	
	    // Do initial hashes
	    patternHash = hash(pattern, pattern_len);
	    textHash = hash(text, pattern_len);
	
	    //compare pattern hash with text hash.  If it fails, use rolling hash
	    //to rehash the text one character ahead
	    for (i = 0; i < hashesToDo; i++) {
	        if (patternHash == textHash && pattern.equals(text.substring(i, i+pattern_len))) 
	        	return i;
	        textHash = hash_increment(text, i, textHash, pattern_len);
	    } 
	    
	    // Pattern not found so return -1
	    return -1;
	}

	// initial hashing
	public static long hash(String string, int strLength){
	    int sum = 0; 
	    for (int i =0; i < strLength; i++)
	    	sum += string.charAt(i); 
	    return sum % largeModulus;
	}

	//rolling hash function
	public static long hash_increment(String string, int prevIndex, long prevHash, int patLength){
		long newHash = (prevHash - ((int) string.charAt(prevIndex)) + ((int) string.charAt(prevIndex + patLength))) % largeModulus;
		if (newHash < 0)
			return newHash + largeModulus;
		else return newHash;
	}
	/**
	 * Method gets contents of a file and converts them to a String
	 * 
	 * @param aFile
	 *            File that has been opened for reading
	 * @return String with the contents of the file
	 */
	public static String getContents(File aFile) {
		StringBuilder contents = new StringBuilder();

		try {

			BufferedReader input = new BufferedReader(new FileReader(aFile));
			try {
				String line = null; // not declared within while loop
				while ((line = input.readLine()) != null) {
					contents.append(line);
					contents.append(System.getProperty("line.separator"));
				}
			} finally {
				input.close();
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return contents.toString();
	}
}
package RabinKarp;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.general.DefaultPieDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.File;
import java.io.IOException;
public class XYChart {
	public XYChart(long[] modulusArray, double[] execTime){
		plot(modulusArray, execTime);
	}
    public static void plot(long[] modulusArray, double[] execTime) {
        //         Create a simple XY chart
        XYSeries series1 = new XYSeries("Small Modulus");
        XYSeries series2 = new XYSeries("Large Modulus");
        for (int i = 0; i < 15; i++){
        	series1.add(modulusArray[i], execTime[i]);
        }
        for (int i = 15; i < 30; i++)
        	series2.add(i-13, execTime[i]);
        //         Add the series to your data set
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series1);
        dataset.addSeries(series2);
        //         Generate the graph
        JFreeChart chart = ChartFactory.createXYLineChart("Figure 3 (3.3 MB DNA seq)", // Title
                "Modulus if small or index of a large modulus", // x-axis Label
                "Execution Time (millis)", // y-axis Label
                dataset, // Dataset
                PlotOrientation.VERTICAL, // Plot Orientation
                true, // Show Legend
                true, // Use tooltips
                false // Configure chart to generate URLs?
            );
        try {
            ChartUtilities.saveChartAsJPEG(new File("/home/roman/Desktop/chart.jpg"), chart, 500,
                300);
        } catch (IOException e) {
            System.err.println("Problem occurred creating chart.");
        }
    }
}



