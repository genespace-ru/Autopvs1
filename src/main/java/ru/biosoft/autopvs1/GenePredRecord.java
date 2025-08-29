package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GenePredRecord {
	String bin;
	String name;
	String chrom;
	String strand;
	int txStart, txEnd;
	int cdsStart, cdsEnd;
	int exonCount;
	
	
	//zero based half-open (deduced by comparing gpe file with gff file)
	int[] exonStarts;
	int[] exonEnds;
	
	int score;
	String name2;
	String cdsStartStat;
	String cdsEndStat;
	String exonFrames;
	
	//GenePred extension format:
	public static List<GenePredRecord> readGenePredFile(File gpeFile) throws IOException
	{
		List<GenePredRecord> result = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(gpeFile));
		String line;
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith("#"))
				continue;
			String[] parts = line.split("\t");
			
			GenePredRecord t = new GenePredRecord();
			t.bin = parts[0];
			t.name = parts[1];
			t.chrom = parts[2];
			t.strand = parts[3];
			t.txStart = Integer.parseInt(parts[4]);
			t.txEnd = Integer.parseInt(parts[5]);
			t.cdsStart = Integer.parseInt(parts[6]);
			t.cdsEnd = Integer.parseInt(parts[7]);
			t.exonCount = Integer.parseInt(parts[8]);
			
			String[] xs = parts[9].split(",");
			t.exonStarts = new int[xs.length];
			for(int i = 0; i < xs.length; i++)
				t.exonStarts[i] = Integer.parseInt(xs[i]);
			
			xs = parts[10].split(",");
			t.exonEnds = new int[xs.length];
			for(int i = 0; i < xs.length; i++)
				t.exonEnds[i] = Integer.parseInt(xs[i]);
			
			t.score = Integer.parseInt(parts[11]);
			t.name2 = parts[12];
			
			t.cdsStartStat = parts[13];
			t.cdsEndStat = parts[14];
			t.exonFrames = parts[15];
			result.add(t);
		}
		return result;
	}
}
