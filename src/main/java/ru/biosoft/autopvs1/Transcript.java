package ru.biosoft.autopvs1;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class Transcript {
	String name,version;
	String gene;
	String chrom;
	boolean forwardStrand;
	int txStart, txEnd;
	int cdsStart, cdsEnd;//If no CDS cdsStart==cdsEnd
	int[] exonStarts, exonEnds;//zero based, half open
	String[] exonFrames;
	
	public static Transcript fromGenePred(GenePredRecord r)
	{
		Transcript t= new Transcript();
		String[] nameParts = r.name.split("[.]");
		if(nameParts.length == 1)
			t.name = nameParts[0];
		else if(nameParts.length == 2)
		{
			t.name = nameParts[0];
			t.version = nameParts[1];
		}else
			throw new IllegalArgumentException("Invalid transcript name: " + r.name);
		t.gene = r.name2;
		t.chrom = r.chrom;
		t.forwardStrand = r.strand.equals("+");
		t.txStart = r.txStart;
		t.txEnd = r.txEnd;
		t.cdsStart= r.cdsStart;
		t.cdsEnd = r.cdsEnd;
		
		t.exonEnds = r.exonEnds;
		t.exonStarts = r.exonStarts;
		/*
		t.exonStarts = Arrays.copyOf( r.exonStarts, r.exonStarts.length );
		t.exonEnds = Arrays.copyOf( r.exonEnds, r.exonEnds.length );
		t.exonFrames = r.exonFrames.split(",");
		if(!t.forwardStrand)
		{
			reverse(t.exonStarts);
			reverse(t.exonEnds);
			reverse(t.exonFrames);
		}
		*/
		return t;
	}
	
	public String getFullName()
	{
		if(version != null)
			return name + "." + version;
		else 
			return name;
	}

	private static <T> void reverse(T[] xs) {
		for(int i = 0; i < xs.length/2; i++)
		{
			int j = xs.length - i - 1;
			T tmp = xs[i];
			xs[i] = xs[j];
			xs[j] = tmp;
		}
	}
	private static void reverse(int[] xs) {
		for(int i = 0; i < xs.length/2; i++)
		{
			int j = xs.length - i - 1;
			int tmp = xs[i];
			xs[i] = xs[j];
			xs[j] = tmp;
		}
	}

	public List<Integer> getCDSSizes() {
		List<Integer> result = new ArrayList<>();
		for(int i = 0; i < exonStarts.length; i++)
		{
			if(exonEnds[i] > cdsStart && exonStarts[i] < cdsEnd)
			{
				int start = exonStarts[i];
				int end = exonEnds[i];
				if(start < cdsStart)
					start = cdsStart;
				if(end > cdsEnd)
					end = cdsEnd;
				result.add(end-start);
			}
			else
				result.add(0);
		}
		if(!forwardStrand)
			Collections.reverse(result);
		return result;
	}
	
	public int getCDSLength()
	{
		int res = 0;
		for(Integer s : getCDSSizes())
			res += s;
		return res;
	}

	//Comment from python code: return a list of intron position [(start, end], (start, end], ... ]
	
	//Actually: zero based [start,end)
	public List<Interval> getIntrons() {
		List<Interval> results = new ArrayList<>();
		for(int i = 0; i < exonStarts.length-1; i++)
		{
			Interval intron = new Interval(exonEnds[i], exonStarts[i+1]);
			results.add(intron);
		}
		return results;
	}

	public List<Interval> getCodingExons() {
		List<Interval> result = new ArrayList<>();
		for(int i = 0; i < exonStarts.length; i++)
		{
			int start = exonStarts[i];
			int end = exonEnds[i];
			if(start >= cdsEnd || end <= cdsStart)
			{
				result.add(null);
				continue;
			}
			if(start < cdsStart)
				start = cdsStart;
			if(end > cdsEnd)
				end = cdsEnd;
			result.add(new Interval(start, end));
		}
		if(!forwardStrand)
			Collections.reverse(result);
		return result;
	}

	public byte[] getCodingSeq(Map<String, byte[]> genome) {
		byte[] chrSeq = genome.get(chrom);
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		for(Interval cExon : getCodingExons())
		{
			if(cExon == null)
				continue;
			byte[] seq = Arrays.copyOfRange(chrSeq, cExon.start, cExon.end);
			if(!forwardStrand)
				seq = Fasta.rc(seq);
			try {
				buffer.write(seq);
			} catch (IOException e) {
				throw new AssertionError();
			}
		}
		return buffer.toByteArray();
	}

	public List<Interval> getCDSList() {
		List<Interval> result = new ArrayList<>();
		for(int i = 0; i < exonStarts.length; i++)
		{
			int eStart = exonStarts[i];
			int eEnd = exonEnds[i];
			if(eStart < cdsEnd && eEnd > cdsStart)//exon overlaps CDS
			{
				if(eStart < cdsStart)
					eStart = cdsStart;
				if(eEnd > cdsEnd)
					eEnd = cdsEnd;
				result.add(new Interval(eStart, eEnd));
			}
		}
		return result;
	}
}
