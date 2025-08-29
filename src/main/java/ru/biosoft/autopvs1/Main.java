package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class Main {
	String inVCF;
	String genomeVersion;
	String vepAssembly;
	String outFolder;
	String providedVEPResult;
	
	Options options;
	Config config;

	public static void main(String[] args) throws Exception {
		Main main = new Main();
		try {
			main.parseArgs(args);
		} catch (ParseException e) {
			main.usage();
			System.exit(1);
		}
		main.run();
	}

	private void parseArgs(String[] args) throws ParseException {
             CommandLineParser parser = new DefaultParser();
             options = new Options();

             Option opt = new Option("g", "genome", true, "genome (hg38,hg19)");
             opt.setRequired(true);
             opt.setOptionalArg(false);
             options.addOption(opt);

             opt = new Option("i", "input-vcf", true, "input vcf file");
             opt.setRequired(true);
             opt.setOptionalArg(false);
             options.addOption(opt);
             
             opt = new Option("o", "out", true, "output folder");
             opt.setRequired(true);
             opt.setOptionalArg(false);
             options.addOption(opt);
             
             opt = new Option("v", "vep", true, "VEP result file, optional. Will run VEP if not provided");
             opt.setRequired(false);
             opt.setOptionalArg(false);
             options.addOption(opt);
             


             CommandLine parsed = parser.parse(options, args);
             String[] cmd = parsed.getArgs();
             if (cmd.length > 0) {
                     usage();
                     System.exit(1);
             }

             genomeVersion = parsed.getOptionValue('g');
             if(genomeVersion.equals("hg19") || genomeVersion.equals("GRCh37"))
             {
            	 genomeVersion = "hg19";
            	 vepAssembly = "GRCh37";
             }else if(genomeVersion.equals("hg38") || genomeVersion.equals("GRCh38"))
             {
            	 genomeVersion = "hg38";
            	 vepAssembly = "GRCh38";
             }else
            	 throw new ParseException("Invalid genome version: " + genomeVersion);
             inVCF = parsed.getOptionValue('i');
             outFolder = parsed.getOptionValue('o');
             
             providedVEPResult = parsed.getOptionValue('v');
     }

	private void usage() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("autopvs1", options);
	}
	
	List<Variant> snps;
	Map<String, Variant> snpsById;
	File vcfForVEP;
	File vepOutput;
	

	private void run()  throws Exception {
		loadConfig();
		readData();
		loadVCF();
		Files.createDirectories(Paths.get(outFolder));
		prepareVCFforVEP();
		runVEP();
		filterVEP();
		
		Map<String, Transcript> transcripts = genomeVersion.equals("hg19")?hg19Transcripts:hg38Transcripts;
		Iterator<VEPResult> it = filteredVEPResults.iterator();
		while(it.hasNext())
		{
			VEPResult r = it.next();
			r.transcript = transcripts.get(r.transcriptId);
			if(r.transcript == null)
			{
				System.err.println("No exon/intron structure for transcript " + r.transcriptId + ", will be skipped");
				it.remove();
			}
		}
		
		
		PVS1 pvs1 = new PVS1(genomeVersion, config);
		
		
		File out = new File(outFolder, "results.txt");
		BufferedWriter writer = new BufferedWriter(new FileWriter(out));
		writer.append("AUTO_ID\tVARIANT_ID\tHGVS_C\tHGVS_p\tCONSEQUENCE\tCRITERION\tSTRENGTH_RAW\tSTRENGTH\n");
		
		//TODO: check other consequences (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences)
		//For example transcript_ablation, stop_lost ???
		Set<String> lof_type = new HashSet<>(Arrays.asList("frameshift", "nonsense", "splice-5", "splice-3", "init-loss"));
		
		for(VEPResult r : filteredVEPResults) {
			try {
				r.conseqence = translateVEPConsequence(r.vepConsequence);
				r.isLOF = lof_type.contains(r.conseqence);
				if(r.isLOF) {
					r.pvs1 = pvs1.run(r);
					//Should we output strength=UNMET and UNSET ???
					writer.append(r.variant.autoId + "\t" + r.variant.id + "\t" + r.hgvs_c + "\t" + r.hgvs_p + "\t"
						+ r.conseqence + "\t" + r.pvs1.criterion + "\t" + r.pvs1.strength_raw + "\t" + r.pvs1.strength + "\n");
				}
			} catch(Exception e)
			{
				System.err.println("Error processing " + r.variant.id + " (autoId=" + r.variant.autoId + ")");
				e.printStackTrace();
			}
		}
		writer.close();

	}

	private String translateVEPConsequence(String consequence) {
		if(consequence.contains("frameshift"))
			return "frameshift";
		if(consequence.contains("stop_gained"))
			return "nonsense";
		//From VEP documentation: A splice variant that changes the 2 base region at the 3' end of an intron
		//But in AutoPVS1 we analyze broader region ???
		if(consequence.contains("splice_donor"))
			return "splice-5";
		if(consequence.contains("splice_acceptor"))
			return "splice-3";
		if(consequence.contains("start_lost"))
			return "init-loss";
					
		return consequence;
	}

	private void readData() throws IOException {
		readGeneTrans();
		readGeneAlias();
		hg19Transcripts = readTranscripts(new File(config.hg19.transcript));
		hg38Transcripts = readTranscripts(new File(config.hg38.transcript));
	}
	
	Map<String, String> geneAlias = new HashMap<>();
	private void readGeneAlias() throws IOException {
		// data/hgnc.symbol.previous.tsv
		BufferedReader reader = new BufferedReader(new FileReader(config.geneAlias));
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			geneAlias.put(parts[1], parts[0]);
		}
		reader.close();	
		
	}
	
	Map<String,Transcript> hg19Transcripts;
	Map<String,Transcript> hg38Transcripts;
	private Map<String,Transcript> readTranscripts(File gpeFile) throws IOException
	{
		List<GenePredRecord> genePredList = GenePredRecord.readGenePredFile(gpeFile);
		Map<String, Transcript> result = new HashMap<>();
		for(GenePredRecord r : genePredList)
		{
			Transcript t = Transcript.fromGenePred(r);
			result.put(t.name, t);
			result.put(t.getFullName(), t);
		}
		return result;
	}

	Map<String, String> geneTrans = new HashMap<>();
	Map<String, String> transGene = new HashMap<>();
	private void readGeneTrans() throws IOException
	{
		//gene_trans = data/clinvar_trans_stats.tsv
		//Странный источник соответствия ген-транскрипт??? откуда от получен?
		//This file gives one main transcript per gene
		BufferedReader reader = new BufferedReader(new FileReader(config.geneTrans));
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			String gene = parts[0];
			String transcript = parts[1];
			transGene.put(transcript, gene);
			geneTrans.put(gene, transcript);
		}
		reader.close();	
	}

	List<VEPResult> filteredVEPResults = new ArrayList<>();
	private void filterVEP() throws IOException {
		//VEP outputs all consequences for all transcripts overlapping variant
		//Here we choose only one transcript and one consequence for further analysis.
		BufferedReader reader = new BufferedReader(new FileReader(vepOutput));
		String line;
		String[] header;
		Map<String, Integer> colIdx = new HashMap<>();
		Map<String, List<VEPResult>> vepResultsByVarId = new HashMap<>();
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith("##"))
				continue;
			if(line.startsWith("#"))
			{
				header = line.split("\t");
				header[0] = header[0].substring(1);//remove #
				for(int i = 0; i < header.length; i++)
					colIdx.put(header[i], i);
				continue;
			}
			String[] parts = line.split("\t");
			VEPResult vep = new VEPResult();
			vep.variationId = parts[colIdx.get("Uploaded_variation")];
			vep.geneSymbol = parts[colIdx.get("SYMBOL")];
			vep.transcriptId = parts[colIdx.get("Feature")];
			vep.canonical = parts[colIdx.get("CANONICAL")].equals("YES");
			vep.pick = parts[colIdx.get("PICK")];
			vep.vepConsequence = parts[colIdx.get("Consequence")].replace("_variant", "");//???
			vep.hgvs_c = parts[colIdx.get("HGVSc")];
			vep.hgvs_p = parts[colIdx.get("HGVSp")].replace("%3D", "=");//???
			vep.hgvs_g = parts[colIdx.get("HGVSg")];
			vep.exon = parts[colIdx.get("EXON")];
			vep.intron = parts[colIdx.get("INTRON")];
			vep.variant = snpsById.get(vep.variationId);
			if(vep.variant == null)
				throw new AssertionError();
			if(vep.transcriptId.equals("-"))//ignore intergenic variants reported by VEP
				continue;
			fixGeneSymbol(vep);
			vepResultsByVarId.computeIfAbsent(vep.variationId, k->new ArrayList<>()).add(vep);
		}
		reader.close();
		
		for(List<VEPResult> rs : vepResultsByVarId.values())
		{
			List<VEPResult> finalChoose = new ArrayList<>();
			for(VEPResult r : rs)
			{
				String mainTranscript = geneTrans.get(r.geneSymbol);
				if(mainTranscript == null)
					continue;
				String transcriptNoVersion = r.transcriptId.split("[.]")[0];
				String mainTranscriptNoVersion = mainTranscript.split("[.]")[0];
				if(r.transcriptId.equals( mainTranscript ) 
				  || transcriptNoVersion.equals( mainTranscriptNoVersion ))
					finalChoose.add(r);
			}
			VEPResult finalVEPResult = null;
			if(finalChoose.size() > 1)
			{
				for(VEPResult r : finalChoose)
					if(r.pick.equals("1"))
						finalVEPResult = r;//What if many picked transcripts? Better to check it then selecting arbitrary result
				if(finalVEPResult == null)
					finalVEPResult = finalChoose.get(0);//Why there is no pick? Better to throw exception.
			}else if(finalChoose.size() == 1)
			{
				finalVEPResult = finalChoose.get(0);
			}else
			{
				finalVEPResult = rs.get(0);//arbitrary selection???
			}
			filteredVEPResults.add(finalVEPResult);
		}
	}

	public void fixGeneSymbol(VEPResult vep) {
		if(vep.geneSymbol.equals("-"))
			vep.geneSymbol = transGene.getOrDefault(vep.transcriptId, "-");
		if(vep.geneSymbol.equals("-"))
			vep.geneSymbol = transGene.getOrDefault(vep.transcriptId.split("[.]")[0], "-");
		if(geneAlias.containsKey(vep.geneSymbol))
			vep.geneSymbol = geneAlias.get(vep.geneSymbol);
	}

	private void loadConfig() {
		config = new Config();
		config.loadConfig();
	}

	private void runVEP() throws Exception {
		vepOutput = new File(outFolder, "vep_output.txt");
		
		if(providedVEPResult != null)
		{
			Files.copy(Paths.get(providedVEPResult), vepOutput.toPath());
			return;
		}
		
		List<String> vepCommand = new ArrayList<>();
		vepCommand.add("vep");
		vepCommand.add("--offline");
		vepCommand.add("--refseq");
		vepCommand.add("--use_given_ref");
		vepCommand.add("--dir_cache");
		vepCommand.add(config.vepCache);
		vepCommand.add("--species");
		vepCommand.add("homo_sapiens");
		vepCommand.add("--assembly");
		vepCommand.add(vepAssembly);
		vepCommand.add("--fork");
		vepCommand.add("4");
		vepCommand.add("--canonical");
		vepCommand.add("--flag_pick");
		vepCommand.add("--hgvs");
		vepCommand.add("--hgvsg");
		vepCommand.add("--symbol");
		vepCommand.add("--distance");
		vepCommand.add("500");
		vepCommand.add("--exclude_predicted");//When using the RefSeq or merged cache, exclude predicted transcripts (i.e. those with identifiers beginning with "XM_" or "XR_")
		vepCommand.add("--lookup_ref");//I think --check-ref will be more robust
		vepCommand.add("--numbers");
		vepCommand.add("--force");
		vepCommand.add("--input_file");
		vepCommand.add(vcfForVEP.getPath());
		vepCommand.add("--output_file");
		vepCommand.add(vepOutput.getPath());
		vepCommand.add("--no_stats");
		vepCommand.add("--tab");
		vepCommand.add("--fields");
		vepCommand.add("Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,Consequence,HGVSc,HGVSp,HGVSg,EXON,INTRON");
		
		System.out.println("Starting vep");
		System.out.println(String.join(" ", vepCommand));
		
		ProcessBuilder pb = new ProcessBuilder(vepCommand);
		pb.redirectOutput(new File(outFolder, "vep.stdout"));
		pb.redirectError(new File(outFolder, "vep.stderr"));
		Process proc = pb.start();
		int exitCode = proc.waitFor();
		if(exitCode != 0)
			throw new RuntimeException("vep exit with " + exitCode);
	}

	private void prepareVCFforVEP() throws IOException {
		vcfForVEP = new File(outFolder, "vep_input.vcf");
		BufferedWriter writer = new BufferedWriter(new FileWriter(vcfForVEP));
		writer.append("##fileformat=VCFv4.2\n");
		writer.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		for(Variant r : snps)
		{
			writer.append(r.chr).append('\t')
			.append(String.valueOf(r.pos)).append('\t')
			.append(String.valueOf(r.autoId)).append('\t')
			.append(r.ref).append('\t')
			.append(r.alt).append('\t')
			.append('.').append('\t') //QUAL
			.append("PASS").append('\t') // FILTER
			.append('.') //INFO
			.append('\n');
		}
		writer.close();
	}

	private void loadVCF() throws IOException {
		snps = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(inVCF));
		String line;
		int lastId = 0;
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith("#"))
				continue;
			String[] parts = line.split("\t");
			if(parts.length < 8)
				throw new RuntimeException("Unexpected line in the input VCF: " + line);
			String chr = parts[0];
			int pos = Integer.parseInt(parts[1]);//1-based in VCF
			String id = parts[2];
			String ref = parts[3];
			String alt = parts[4];
			String qual = parts[5];
			String filter = parts[6];
			String info = parts[7];
			
			for(String singleAlt : alt.split(","))
			{
				Variant v = new Variant();
				v.chr = chr;
				v.pos = pos;
				v.ref = ref;
				v.alt = singleAlt;
				v.id = id;
				v.autoId = ++lastId;
				String msg = v.validate();
				if(msg != null)
				{
					System.err.println("Can not process variant, will be skipped: " + line);
					System.err.println(msg);
					continue;
				}
				snps.add(v);
			}
		}
		
		snpsById = new HashMap<>();
		for(Variant v : snps)
			snpsById.put(String.valueOf(v.autoId), v);
	}
	

}
