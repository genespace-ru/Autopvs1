package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import ru.biosoft.autopvs1.Config.GenomeData;

public class Data {
	
	public Map<String, byte[]> genome;
	public BedFile domain,hotspot,curatedRegion,exonLOFPopmax;
	public Map<String, Map<String, Double>> pathogenicDict;
	public Map<String, String> pvs1Levels;//key is a gene name


	public void load(Config config, String genomeVersion) throws IOException
	{
		GenomeData data;
		if(genomeVersion.equals("hg19"))
		{
			data = config.hg19;
		}else if(genomeVersion.equals("hg38"))
		{
			data = config.hg38;
		}else
			throw new IllegalArgumentException("Invalid genome version: " + genomeVersion);
		
		genome = Fasta.loadFasta(new File(data.genome));
		domain = BedFile.load(new File(data.domain));
		hotspot = BedFile.load(new File(data.hotspot));
		curatedRegion = BedFile.load(new File(data.curatedRegion));
		exonLOFPopmax = BedFile.load(new File(data.exonLofPopmax));
		loadPathogenicFile(new File(data.pathogenicSite));
		loadPVS1Levels(new File(config.pvs1Levels));
	}


	private void loadPVS1Levels(File file) throws IOException {
		pvs1Levels = new HashMap<>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while((line = reader.readLine()) != null) {
			String[] parts = line.split("\t");
			String geneName = parts[0];
			String level = parts[1];//What is L1E,L2E,L3E ???
			if(parts.length > 2)
				throw new IllegalArgumentException();
			pvs1Levels.put(geneName, level);
		}
		reader.close();
	}


	private void loadPathogenicFile(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		pathogenicDict = new HashMap<>();
		pathogenicDict.put("score", new HashMap<>());
		pathogenicDict.put("count", new HashMap<>());
		while((line = reader.readLine()) != null) {
			if(line.startsWith("#"))
				continue;
			String[] parts = line.split("\t");
			String key = parts[0] + ":" + parts[1];
			double score;
			switch(parts[6])
			{
			case "4": case "3": case "2": score = 1; break;
			case "1": score = 0.5; break;
			default: score = 1/3;
			}
			pathogenicDict.get("score").compute(key, (k,v)->v==null?score:v+score);
			pathogenicDict.get("count").compute(key, (k,v)->v==null?1:v+1);

		}
		reader.close();
	}
}
