package ru.biosoft.autopvs1;

public class Config {

	public String vepCache;
	public String pvs1Levels;
	public String geneAlias;
	public String geneTrans;
	
	public GenomeData hg19;
	public GenomeData hg38;
	
	public static class GenomeData
	{
		public String genome;
		public String transcript;
		public String domain;
		public String hotspot;
		public String curatedRegion;
		public String exonLofPopmax;
		public String pathogenicSite;
	}
	
	public void loadConfig()
	{
		vepCache = System.getProperty("user.home") + "/.vep";
		pvs1Levels = "data/PVS1.level";
		geneAlias = "data/hgnc.symbol.previous.tsv";
		
		
		/*
To annotate the transcript, we first retrieved the most prevalent
transcripts used in ClinVar (Landrum et al., 2018). The ClinVar 20200106 release was
used throughout this study. If no result was found in ClinVar, we retrieved transcripts
from RefSeqGene/the Locus Reference Genomic collaborations
(https://www.lrg-sequence.org). If both retrievals were unsuccessful, we used the
longest available transcript in the Reference Sequence (RefSeq) database (O'Leary et
al., 2016).
		 */
		geneTrans = "data/clinvar_trans_stats.tsv";
		
		
		hg19 = new GenomeData();
		hg19.genome = "data/hg19.fa";
		
		//The source of this file is unknown, format also not clear. Better to replace with https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gff.gz
		hg19.transcript = "data/ncbiRefSeq_hg19.gpe";
		hg19.domain = "data/functional_domains_hg19.bed";
		hg19.hotspot = "data/mutational_hotspots_hg19.bed";
		hg19.curatedRegion = "data/expert_curated_domains_hg19.bed";
		hg19.exonLofPopmax = "data/exon_lof_popmax_hg19.bed";
		hg19.pathogenicSite = "data/clinvar_pathogenic_GRCh37.vcf";
		
		hg38 = new GenomeData();
		hg38.genome = "data/hg38.fa";
		hg38.transcript = "data/ncbiRefSeq_hg38.gpe";
		hg38.domain = "data/functional_domains_hg38.bed";
		hg38.hotspot = "data/mutational_hotspots_hg38.bed";
		hg38.curatedRegion = "data/expert_curated_domains_hg38.bed";
		hg38.exonLofPopmax = "data/exon_lof_popmax_hg38.bed";
		hg38.pathogenicSite = "data/clinvar_pathogenic_GRCh38.vcf";
	}
}
