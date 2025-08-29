package ru.biosoft.autopvs1;

/*
#Uploaded_variation  SYMBOL  Feature         CANONICAL  PICK  Consequence          HGVSc                    HGVSp                       HGVSg              EXON  INTRON
1                    F10     NM_000504.4     YES        1     stop_gained          NM_000504.4:c.1043G>A    NP_000495.1:p.Trp348Ter     13:g.113149093G>A  8/8   -
1                    F10     NM_001312674.2  -          -     stop_gained          NM_001312674.2:c.911G>A  NP_001299603.1:p.Trp304Ter  13:g.113149093G>A  7/7   -
1                    F10     NM_001312675.2  -          -     3_prime_UTR_variant  NM_001312675.2:c.*34G>A  -                           13:g.113149093G>A  8/8   -

 * 
            self.vep_variation = final.record['Uploaded_variation']
            self.vep_symbol = final.record['SYMBOL']
            self.vep_trans = final.record['Feature']
            self.vep_canonical = final.record['CANONICAL']
            self.vep_pick = final.record['PICK']
            self.vep_consequence = final.record['Consequence'].replace('_variant', '')
            self.hgvs_c = final.record['HGVSc']
            self.hgvs_p = final.record['HGVSp'].replace('%3D', '=')
            self.hgvs_g = final.record['HGVSg']
            self.vep_exon = final.record['EXON']
            self.vep_intron = final.record['INTRON']

 */
public class VEPResult {
	String variationId;
	String geneSymbol;
	String transcriptId;
	boolean canonical;//a flag indicating if the transcript is denoted as the canonical transcript for this gene
	String pick;
	String vepConsequence;
	String hgvs_c;
	String hgvs_p;
	String hgvs_g;
	String exon;
	String intron;
	
	//extra data
	Variant variant;//input variant
	Transcript transcript;
	String conseqence;
	boolean isLOF;
	public PVS1Result pvs1;
}
