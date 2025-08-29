package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class Maxent {

	public Maxent() throws IOException
	{
		loadMatrix5();
		loadMatrix3();
	}
	
	double[] bgd_5 = {0.27, 0.23, 0.23, 0.27};
	double[] cons1_5 = {0.004, 0.0032, 0.9896, 0.0032};
	double[] cons2_5 = {0.0034, 0.0039, 0.0042, 0.9884};

	double[] bgd_3 = {0.27, 0.23, 0.23, 0.27};
	double[] cons1_3 = {0.9903, 0.0032, 0.0034, 0.0030};
	double[] cons2_3 = {0.0027, 0.0037, 0.9905, 0.0030};

	
	private Map<String, Double> matrix5 = new HashMap<>();
	
	
	private void loadMatrix5() throws IOException {
		File file = new File("data/score5_matrix.txt");
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("[\t ]+");
			String seq = parts[0];
			double val = Double.parseDouble(parts[1]);
			matrix5.put(seq, val);
		}
		reader.close();
	}

	//Calculate 5' splice site strength
	//(exon)XXX|XXXXXX(intron)
    //          **

	public double score5(byte[] seq)
	{
		if(seq.length != 9)
			throw new IllegalArgumentException();
		byte c1 = Fasta.code(seq[3]);
		byte c2 = Fasta.code(seq[4]);
		double score = cons1_5[c1]*cons2_5[c2]/(bgd_5[c1]*bgd_5[c2]);
		byte[] rest = new byte[seq.length - 2];
		for(int i = 0; i < 3; i++)
			rest[i] = seq[i];
		for(int i = 5; i < seq.length; i++)
			rest[i-2] = seq[i];
		String key = new String(rest).toUpperCase();
		double restScore = matrix5.get(key);
		return Math.log(score*restScore)/Math.log(2);
	}
	
	double[][] matrix3 = new double[9][16384];
	
	private void loadMatrix3() throws IOException
	{
		File file = new File("data/score3_matrix.txt");
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("[\t ]+");
			int i = Integer.parseInt(parts[0]);
			int j = Integer.parseInt(parts[1]);
			double val = Double.parseDouble(parts[2]);
			matrix3[i][j] = val;
		}
		reader.close();
	}
	
	public double score3(byte[] seq)
	{
		if(seq.length != 23)
	        throw new IllegalArgumentException();
	    byte c1 =  Fasta.code(seq[18]);
	    byte c2 = Fasta.code(seq[19]);
	    double score = cons1_3[c1] * cons2_3[c2] / (bgd_3[c1] * bgd_3[c2]);

	    byte[] rest = new byte[seq.length - 2];
		for(int i = 0; i < 18; i++)
			rest[i] = seq[i];
		for(int i = 20; i < seq.length; i++)
			rest[i-2] = seq[i];
		
	    double rest_score = 1;
	    rest_score *= matrix3[0][hashseq(rest,0,7)];
	    rest_score *= matrix3[1][hashseq(rest,7,7)];
	    
	    //This was in the original autopvs1, will be out of bounds ???
	    //rest_score *= matrix3[2][hashseq(rest,14,9)]; 
	    
	    //This is my version
	    rest_score *= matrix3[2][hashseq(rest,14,7)];
	    
	    rest_score *= matrix3[3][hashseq(rest,4,7)];
	    rest_score *= matrix3[4][hashseq(rest,11,7)];
	    rest_score /= matrix3[5][hashseq(rest,4,3)];
	    rest_score /= matrix3[6][hashseq(rest,7,4)];
	    rest_score /= matrix3[7][hashseq(rest,11,3)];
	    rest_score /= matrix3[8][hashseq(rest,14,4)];
	    return Math.log(score * rest_score)/Math.log(2);

	}
	
	int hashseq(byte[] seq, int offset, int length)
	{
		int res = 0;
		int C = 1;
		for(int i = offset+length-1; i >= offset; i--)
		{
			res += Fasta.code(seq[i]) * C;
			C *= 4;
		}
		return res;
	}
}
