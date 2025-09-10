package ru.biosoft.autopvs1;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class PVS1 {

	String genomeVersion;
	Config config;
	
	Data data;
	Splicing splicing;
	
	public PVS1(String genomeVersion, Config config) throws IOException {
		this.genomeVersion = genomeVersion;
		this.config = config;
		loadData();
		splicing = new Splicing(data);
	}

	private void loadData() throws IOException {
		data = new Data();
		data.load(config, genomeVersion);
	}

	public PVS1Result run(VEPResult r) {
		PVS1Result res = new PVS1Result();
		verify(r, res);
		adjust(r, res);
		return res;
	}

	private void adjust(VEPResult r, PVS1Result res) {
		if(r.transcript == null)
		{
			res.strength = ACMGStrength.UNSET;
			return;
		}
		String geneName = r.transcript.gene;
		if(geneName.equals("MYH7"))
		{
			if(res.strength_raw.compareTo(ACMGStrength.MODERATE) >= 0)
				res.strength = ACMGStrength.MODERATE;
			else
				res.strength = res.strength_raw;
		}
		else if(data.pvs1Levels.containsKey(geneName))
		{
			String level = data.pvs1Levels.get(geneName);
			switch(level)
			{
			case "L0": res.strength = res.strength_raw; break;
			case "L1": res.strength = res.strength_raw.downgrade(1); break;
			case "L2": res.strength = res.strength_raw.downgrade(2); break;
			case "L3": res.strength = ACMGStrength.UNMET; break;
			default://L1E,L2E,L3E not handled ???
				res.strength = ACMGStrength.UNSET;
			}
		}else
		{
			res.strength = ACMGStrength.UNSET;
		}
	}

	private void verify(VEPResult r, PVS1Result res) {
		if(r.transcript == null || r.transcript.cdsStart == r.transcript.cdsEnd)
		{
			res.strength_raw = ACMGStrength.UNMET;
			res.criterion = "NF0";
			return;
		}
		switch(r.conseqence)
		{
		case "nonsense": case "frameshift":
			if(r.transcript.gene.equals("PTEN") && get_pHGVS_termination(r.hgvs_p) < 374)//Any nonsense or truncation to shorter then 374?
			{
				res.criterion = "PTEN";
				res.strength_raw = ACMGStrength.VERY_STRONG;
				return;
			}
			if(isNMDTartget(r))
			{
				if(isBiologicallyRelevant(r))
				{
					res.criterion = "NF1";
					res.strength_raw = ACMGStrength.VERY_STRONG;
				}
				else
				{
					res.criterion = "NF2";
					res.strength_raw = ACMGStrength.UNMET;
				}
			}else
			{
				if(isCriticalToProteinFunction(r))
				{
					res.criterion = "NF3";
					res.strength_raw = ACMGStrength.STRONG;
				}else
				{
					if(exonLOFAreFrequentInPop(r) || !isBiologicallyRelevant(r))
					{
						res.criterion = "NF4";
						res.strength_raw = ACMGStrength.UNMET;
					}else if(removesMoreThen10PercentOfProtein(r))
					{
						res.criterion = "NF5";
						res.strength_raw = ACMGStrength.STRONG;
					}else
					{
						res.criterion = "NF6";
						res.strength_raw = ACMGStrength.MODERATE;
					}
				}
			}
			break;
		case "splice-5": case "splice-3":
			if(r.transcript.gene.equals("PTEN"))
			{
				Pattern p = Pattern.compile("c[.]([0-9]+)");
				Matcher m = p.matcher(r.hgvs_c);
				if(m.find())
				{
					int pos = Integer.parseInt(m.group(1));
					pos /= 3;
					if(pos < 374)
					{
						res.criterion = "PTEN";
						res.strength_raw = ACMGStrength.VERY_STRONG;
						return;
					}
				}
			}
			if(r.transcript.gene.equals("CDH1"))
			{
				res.criterion = "CDH1";
				res.strength_raw = ACMGStrength.STRONG;
				return;
			}
			splicing(r, res);
			break;
		case "init-loss":
			startCodon(r, res);
			break;
		default:
				res.criterion = "IC5";
				res.strength_raw = ACMGStrength.UNMET;
		}
	}

	
	private void splicing(VEPResult r, PVS1Result res) {
		splicing.calculate(r);
		if(splicing.isPreserveReadingFrame())
		{
			if(splicing.isCriticalToProteinFunction())
			{
				res.criterion = "SS10";
				res.strength_raw = ACMGStrength.STRONG;
			}else if(exonLOFAreFrequentInPop(r) || !isBiologicallyRelevant(r))
			{
				res.criterion = "SS7";
				res.strength_raw = ACMGStrength.UNMET;
			}else if(splicing.variantRemoves10PercentOfProtein())
			{
				res.criterion = "SS8";
				res.strength_raw = ACMGStrength.STRONG;
				
			}else
			{
				res.criterion = "SS9";
				res.strength_raw = ACMGStrength.MODERATE;
			}
		}else if(splicing.isUndergoNMD())
		{
			if(isBiologicallyRelevant(r))
			{
				res.criterion = "SS1";
				res.strength_raw = ACMGStrength.VERY_STRONG;
			}else
			{
				res.criterion = "SS2";
				res.strength_raw = ACMGStrength.UNMET;
			}
		}
		else
		{
			if(splicing.isCriticalToProteinFunction())
			{
				res.criterion = "SS3";
				res.strength_raw = ACMGStrength.STRONG;
			}else if(exonLOFAreFrequentInPop(r) || !isBiologicallyRelevant(r))
			{
				res.criterion = "SS4";
				res.strength_raw = ACMGStrength.UNMET;
			}else if(splicing.variantRemoves10PercentOfProtein())
			{
				res.criterion = "SS5";
				res.strength_raw = ACMGStrength.STRONG;
				
			}else
			{
				res.criterion = "SS6";
				res.strength_raw = ACMGStrength.MODERATE;
			}
		}
	}

	private boolean removesMoreThen10PercentOfProtein(VEPResult r) {
		 Pattern pattern = Pattern.compile("p[.][^0-9]+([0-9]+)([^0-9]+fs)?([*]|X|Ter)([0-9]+)?");
		 Matcher matcher = pattern.matcher(r.hgvs_p);
		 int codon_offset = -1;
		 if(matcher.find())
		 {
			 codon_offset = Integer.parseInt(matcher.group(1));
		 }
         int codon_length = r.transcript.getCDSLength() / 3;
         
         return (codon_offset > 0 && ((double)codon_length - codon_offset) / codon_length > 0.1);
	}


	static class ExonLOFPopmax
	{
		boolean value;
		String desc;
	}
	
	private boolean exonLOFAreFrequentInPop(VEPResult r) {
		return getExonLOFPopmaxInfo(r).value;
	}

	private ExonLOFPopmax getExonLOFPopmaxInfo(VEPResult r) {
		ExonLOFPopmax res = new ExonLOFPopmax();
		int start = r.variant.pos - 5;
		int end = r.variant.pos + 5;
		String chrom = r.variant.chr;
		if(!chrom.startsWith("chr"))
			chrom = "chr"+ r.variant.chr;
		List<BedFile.Item> items = data.exonLOFPopmax.index.queryOverlapping(chrom, start, end-1);
		if(items.size() > 0)
		{
			BedFile.Item item = items.get(0);//Use only one match???
			String[] lof_list = item.name.split("[|]");
		    int lof_num = lof_list.length - 1;
		    double sum_freq = 0;
		    for(int i = 1; i < lof_list.length; i++)
		    	sum_freq += Double.parseDouble( lof_list[i].split(":")[1]);
		    String[] parts = lof_list[1].split(":");
		    String max_lof = parts[0];
		    String max_freq = parts[1];
		    parts = lof_list[0].split("[.]");
		    String transcript = parts[0];
		    String version = parts[1];
		    String exon = parts[2];
		    if(Double.parseDouble(max_freq) > 0.001)
		    {
		    	res.value = true;
		    	res.desc = "Maximum LOF population frequency in exon "+exon+" of "+transcript+"."+version+" is "
		    			+ "<a href=\"https://gnomad.broadinstitute.org/variant/"+max_lof+"\">"+max_freq+"</a>, "
		    			+ "higher than the threshold (0.1%) we pre-defined.";
		    	return res;
		    }else
		    {
		    	res.value = false;
		    	res.desc = "Maximum LOF population frequency in exon "+exon+" of "+transcript+"."+version+" is "
		    			+ "<a href=\"https://gnomad.broadinstitute.org/variant/"+max_lof+"\">"+max_freq+"</a>, "
		    			+ "lower than the threshold (0.1%) we pre-defined.";
		    }
		}
		else
		{
			res.value = false;
			res.desc = "No LOF variant is found or the LOF variant dosen't exist in gnomAD.";
		}

		return res;
	}

	// Truncated/altered region is critical to protein function.
	private boolean isCriticalToProteinFunction(VEPResult r) {
		
		return getFunctionalRegion(r).isFunctional;
	}
	
	static class FunctionalRegion
	{
		boolean isFunctional;
		String description;
	}
	
	private FunctionalRegion getFunctionalRegion(VEPResult r)
	{
		FunctionalRegion res = new FunctionalRegion();
		if(r.transcript.gene.equals("CDH1"))
		{
			res.isFunctional = get_pHGVS_termination(r.hgvs_p) <= 836;
			if(res.isFunctional)
				res.description = "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/30311375\">CDH1 gene-specific criteria</a>: " +
                "Truncations in NMD-resistant zone located upstream the most 3′ well-characterized " +
                "pathogenic variant c.2506G>T (p.Glu836Ter). ";
			else
				res.description = "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/30311375\">CDH1 gene-specific criteria</a>: " + 
				"Truncations in NMD-resistant zone located downstream the most 3′ well-characterized " +
				"pathogenic variant c.2506G>T (p.Glu836Ter).";
			return res;

		}
		
		String chrom = r.variant.chr;
		if(chrom.startsWith("chr"))
			chrom = chrom.substring("chr".length());
		int start,end;
		if(r.transcript.forwardStrand) {
			start = r.variant.pos - 1;//make it zero-based as in bed file
			end = r.transcript.cdsEnd;
		}else
		{
			start = r.transcript.cdsStart;
			end = r.variant.pos;
		}
		return getFunctionalRegion(chrom, start, end, data);
	}

	public static FunctionalRegion getFunctionalRegion(String chrom, int start, int end, Data data) {
		if(!chrom.startsWith("chr"))
			chrom = "chr"+chrom;
		FunctionalRegion res = new FunctionalRegion();
		List<BedFile.Item> domains = data.domain.index.queryOverlapping(chrom, start, end-1);
		List<BedFile.Item> hotspots = data.hotspot.index.queryOverlapping(chrom, start, end-1);
		List<BedFile.Item> curatedRegions = data.curatedRegion.index.queryOverlapping(chrom, start, end-1);
		
		if(curatedRegions.size()>0)
		{
			res.isFunctional = true;
			res.description = "";
			for(BedFile.Item item : curatedRegions)
				res.description += "Expert curated region: " + item.name + "\n";
		}else if(hotspots.size() > 0)
		{
			res.isFunctional = true;
			res.description = "";
			for(BedFile.Item item : hotspots)
			{
				String[] parts = item.name.split("[|]");
				String missensePLP = parts[3];
				String missenseBLB = parts[4];
				String genomicPosition = parts[0];
				res.description += "Mutation hotspot: <b>" +missensePLP + "</b> pathogenic missense variants and <b>" + missenseBLB + "</b> benign "
				+"missense variant in <b>"+ genomicPosition +"</b>\n";
			}
				
		}
		if(domains.size() > 0)
		{
			//Should check all matches, but in original autopvs1 only one match used ???
			Collections.sort(domains, Comparator.comparing(item->item.start));
			BedFile.Item item = domains.get(0);
			//for(BedFile.Item item : domains)
			{
				String[] parts = item.name.split("[|]");
				String domain_name = parts[0];
				String amino_acids = parts[1];
				String genomic_position = parts[2];
				String tag = parts[3];
				Integer missense_total = Integer.parseInt(parts[4]);
				int missense_PLP = Integer.parseInt(parts[5]);
				int missense_BLB = Integer.parseInt(parts[6]);
				String block_position = parts[7];
				if ((missense_BLB == 0 && missense_PLP >= 5) || (missense_BLB > 0 && missense_PLP / missense_BLB >= 10) )
					res.isFunctional = true;
				if(missense_total == 0)
					res.description += "No ClinVar missense variant is found in <b>"+domain_name+"</b> domain ("+amino_acids+").\n";
				else
				{
					String[] aaParts = amino_acids.split(" ");
					String uniport_id = aaParts[aaParts.length-1];
					String amino_acid = String.join(" ", Arrays.copyOf(aaParts, aaParts.length-1));
			        res.description += "<b>"+missense_PLP+"</b> ClinVar pathogenic missense variant(s) and <b>"+missense_BLB+"</b> benign missense variant(s) "+
			                        "are found in <b>"+domain_name+"</b> domain ("+amino_acid + " <a href=\"https://www.uniprot.org/uniprot/"+uniport_id+"\">"+uniport_id+"</a>).\n";
				}
			}
		}
		if(hotspots.size()==0 && domains.size() == 0)
			res.description += "Neither mutational hotspot nor functional domain is found.\n";
		return res;
	}

	private boolean isBiologicallyRelevant(VEPResult r) {
		return r.transcript != null; // always true ???
	}

	/**
	 * Function: Nonsense-mediated decay (NMD) classification
	 * See more information about NMD: https://en.wikipedia.org/wiki/Nonsense-mediated_decay 
	 * NMD classify: not occurring in the 3′ most exon or the 3′-most 50 bp of the penultimate exon
	 * @param r
	 * @return
	 */
	private boolean isNMDTartget(VEPResult r) {
		if(r.transcript.gene.equals("GJB2"))
			return true;//Any mutation in this gene leads to NMD??? See https://pmc.ncbi.nlm.nih.gov/articles/PMC11021044/. We should check at least that mutation cause premature stop codon. 
		int newStopCodon = get_pHGVS_termination(r.hgvs_p);//???can be -1s
		List<Integer> cdsSizes = r.transcript.getCDSSizes().stream().filter(l->l>0).collect(Collectors.toList());
		if(cdsSizes.size() <= 1)
			return false;
		int nmdCutoff = 0;
		for(int i = 0; i < cdsSizes.size() - 1; i++)//except last exon
			nmdCutoff += cdsSizes.get(i);
		//subtract 50bp of penultimate exon (or the whole exon if shorter then 50bp)
		nmdCutoff -= Math.min(cdsSizes.get(cdsSizes.size()-2), 50);
		
		return newStopCodon * 3 <= nmdCutoff;
	}

	private int get_pHGVS_termination(String pHGVS) {
		//Парсинг hgvs запутанный, не проще ли посчитать длину нового CDS напрямую?
		if(pHGVS.contains("fs"))
		{
			Pattern pattern1 = Pattern.compile("p[.][^0-9]+([0-9]+)[^0-9]+fs([*]|X|Ter)([0-9]+)");
			Pattern pattern2 = Pattern.compile("p[.][^0-9]+([0-9]+)fs");
			Matcher m1 = pattern1.matcher(pHGVS);
			Matcher m2 = pattern2.matcher(pHGVS);
			if(m1.find())
			{
				return Integer.parseInt(m1.group(1)) + Integer.parseInt(m1.group(3));
			}else if(m2.find())
			{
				return Integer.parseInt(m2.group(1));
			}
			else
				return -1;
		}else if(pHGVS.contains("*") || pHGVS.contains("X") || pHGVS.contains("Ter"))
		{
			Pattern pattern = Pattern.compile("p[.][^0-9]+([0-9]+)([*]|X|Ter)");
		    Matcher matcher = pattern.matcher(pHGVS);
		    if(matcher.find())//???find or matches
		    	return Integer.parseInt(matcher.group(1));
		    return -1;

		}
		else
			return -1;
	}

	private void startCodon(VEPResult r, PVS1Result res) {
		AltStartCodon altStart = findClosestPotentialStartCodon(r.transcript);
		String chr = r.variant.chr;
		if(chr.startsWith("chr"))
			chr = chr.substring("chr".length());
		if(altStart.cDNAPos == -1)
		{
			res.strength_raw = ACMGStrength.VERY_STRONG;
			res.criterion = "IC0";
			return;
		}
		double variantScore = 0;
		for(Interval i : altStart.intervalsBefore)
			for(int x = i.start; x < i.end; x++)
			{
				String key = chr + ":" + x;
				Double count = data.pathogenicDict.get("count").get(key);
				if(count != null)
					variantScore += count;
			}
		
		if(variantScore < 1)
		{
			res.strength_raw = ACMGStrength.SUPPORTING;
			res.criterion = "IC4";
		}else if(variantScore <= 3)
		{
			res.strength_raw = ACMGStrength.MODERATE;
			res.criterion = "IC3";
		}else if(variantScore <= 6)
		{
			res.strength_raw = ACMGStrength.STRONG;
			res.criterion = "IC2";
		}else {
			res.strength_raw = ACMGStrength.VERY_STRONG;
			res.criterion = "IC1";
		}
	}
	
	static class AltStartCodon
	{
		int cDNAPos = -1;
		int genomicPos = -1;
		List<Interval> intervalsBefore = new ArrayList<>();
	}

	private AltStartCodon findClosestPotentialStartCodon(Transcript transcript) {
		byte[] codingSeq = transcript.getCodingSeq(data.genome);
		int altStartCodon = -1;
		int genomicPos = -1;
		List<Interval> altStartInterval = new ArrayList<>();
		
		for(int i = 3; i+2 < codingSeq.length; i+= 3)
		{
			String codon = new String(codingSeq, i, 3);
			if(codon.toUpperCase().equals("ATG"))
			{
				altStartCodon = i;
				break;
			}
		}
		if(altStartCodon != -1)
		{
			List<Interval> cdsList = transcript.getCDSList();
			int cdsSizeTmp = 0;
		
			if(transcript.forwardStrand)
			{
				for(Interval i : cdsList)
				{
					cdsSizeTmp += i.end - i.start;
					if(cdsSizeTmp > altStartCodon)//???error: should be altStartCodon > cdsSizeTmp
					{
						genomicPos = i.end - (cdsSizeTmp - altStartCodon);
						if(i.start < genomicPos)
							altStartInterval.add(new Interval(i.start, genomicPos));
						break;
					}else
						altStartInterval.add(i);
				}
				
			}else { //reverse strand
				Collections.reverse(cdsList);
				for(Interval i : cdsList)
				{
					cdsSizeTmp += i.end - i.start;
					if( cdsSizeTmp > altStartCodon)//???error: should be altStartCodon > cdsSizeTmp
					{
						genomicPos = i.start + (cdsSizeTmp - altStartCodon)-1;
						if(genomicPos+1 < i.end)
							altStartInterval.add(new Interval(genomicPos+1, i.end));
						break;
					}else
						altStartInterval.add(i);
				}
				Collections.reverse(altStartInterval);
			}
		}
		AltStartCodon res = new AltStartCodon();
		res.cDNAPos = altStartCodon;
		res.genomicPos = genomicPos;
		res.intervalsBefore = altStartInterval;
		return res;
	}
	

}
