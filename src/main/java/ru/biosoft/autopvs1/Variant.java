package ru.biosoft.autopvs1;

public class Variant {
	String chr;// chromosome name
	int pos;// 1-based
	
	//Correspond to REF field of VCF 4.2 specification
	//Upper case
	String ref;
	
	//Single alternative allele (multiple alt alleles from VCF should be split into separate variants)
	//Each base must be one of [ACGTN]
	//VCF missing value '*' is not allowed
	//VCF angle-bracketed string is not allowed
	String alt;
	
	String id;//original id (from VCF)
	
	int autoId; //record number

	public String validate() {
		if(!chr.matches("[A-Za-z0-9_.]+"))
			return "Invalid chromosome name: " + chr;
		if(pos < 1)
			return "Invalid position: " + pos;
		if(!ref.matches("[ACGTNacgtn]+"))
			return "Invalid ref: " + ref;
		if(!alt.matches("[ACGTNacgtn]+"))
			return "Invalid alt: " + alt;
		if(id.matches("[;]"))
			return "Invalid id: " + id;
		return null;
	}
}
