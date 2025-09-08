package ru.biosoft.autopvs1;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ru.biosoft.autopvs1.PVS1.FunctionalRegion;

public class Splicing {
	private VEPResult r;
	private Data data;

	String type;
	String index;

	byte[] chrSeq;

	byte[] refseq;
	int refseq_start, refseq_end;// zero-based half open

	byte[] altseq;

	Maxent maxent;
	double maxentscore_ref;
	double maxentscore_alt;
	double maxentfoldchange;

	public Splicing(Data data) throws IOException {
		this.data = data;
		maxent = new Maxent();
	}

	public void calculate(VEPResult r) {
		this.r = r;
		init();
		parse();
		calcMaxEntScore();
	}

	private void init() {
		type = null;
		index = null;
		refseq = altseq = new byte[0];
		refseq_start = refseq_end = 0;
		maxentscore_ref = maxentscore_alt = -1;
		maxentfoldchange = 1;
	}

	private void calcMaxEntScore() {
		// What about alt.length() != ref.length() ???
		if (type.equals("donor")) {
			maxentscore_ref = maxent.score5(refseq);
			if(altseq != null)
				maxentscore_alt = maxent.score5(altseq);
		} else if (type.equals("acceptor")) {
			maxentscore_ref = maxent.score3(refseq);
			if(altseq != null)
				maxentscore_alt = maxent.score3(altseq);
		} else
			throw new IllegalArgumentException();
		if(altseq != null)
			maxentfoldchange = maxentscore_alt / maxentscore_ref;

		maxentscore_ref = Math.round(maxentscore_ref * 100) / 100d;
		maxentscore_alt = Math.round(maxentscore_alt * 100) / 100d;
		maxentfoldchange = Math.round(maxentfoldchange * 100) / 100d;
	}

	// DONOR - 5' site, ACCEPTOR - 3' site
	// Define splice site boundaries
	static final int DONOR_EXON = 3;
	static final int DONOR_INTRON = 6;
	static final int ACCEPTOR_EXON = 3;
	static final int ACCEPTOR_INTRON = 20;
	static final double PERCENT_THRESHOLD = 0.7;
	static final double DONOR_THRESHOLD = 3;
	static final double ACCEPTOR_THRESHOLD = 3;

	private void parse() {
		// extract refseq and altseq of splice site
		String chr = r.variant.chr;
		if (!chr.startsWith("chr"))
			chr = "chr" + chr;// We use ucsc fasta file with chr1,... names
		chrSeq = data.genome.get(chr);

		List<Interval> introns = r.transcript.getIntrons();
		for (int intronIdx = 0; intronIdx < introns.size(); intronIdx++) {
			Interval intron = introns.get(intronIdx);
			for (int pos = r.variant.pos; pos < r.variant.pos + r.variant.ref.length(); pos++) {//in the case of insertion/deletion first pos actually not changed (according to VCF specification) ???
				// pos is one-based, intron is zero based [start,end)
				int distanceToIntronStart = (pos - 1) - intron.start;
				int distanceToIntronEnd = (pos - 1) - (intron.end - 1);
				if (r.transcript.forwardStrand) {
					//Why we search for overlap with a whole splicing site, when PVS1 should only consider mutation of 2-base region of intron???
					if (distanceToIntronStart > -DONOR_EXON && distanceToIntronStart <= DONOR_INTRON)// [-2,-1,0,1,2,3,4,5,6]
					{
						type = "donor";
						refseq_start = intron.start - DONOR_EXON;//zero based, including
						refseq_end = intron.start + DONOR_INTRON;//zero based, excluding
						refseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);

						if(r.variant.ref.length() > 1 || r.variant.alt.length() > 1)//ins or del overlap splice site
						{
							altseq = null;//assume that splice site is always disrupted in this case. We should search for nearby cryptic splice site in alternative sequence later.
						}else
						{
							//single nucleotide variant
							altseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
							altseq[r.variant.pos-1-refseq_start] = (byte) r.variant.alt.charAt(0);
						}

						if (distanceToIntronStart >= 0)
							index = "IVS" + (intronIdx + 1) + "+" + (distanceToIntronStart+1);
						else
							index = "EX" + (intronIdx + 1) + "" + distanceToIntronStart;

					}
					if (-ACCEPTOR_INTRON < distanceToIntronEnd && distanceToIntronEnd <= ACCEPTOR_EXON) {
						type = "acceptor";
						refseq_start = intron.end - ACCEPTOR_INTRON;
						refseq_end = intron.end + ACCEPTOR_EXON;
						refseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
						
						if(r.variant.ref.length() > 1 || r.variant.alt.length() > 1)//ins or del overlap splice site
						{
							altseq = null;//assume that splice site is always disrupted in this case. We should search for nearby cryptic splice site in alternative sequence later.
						}else
						{
							//single nucleotide variant
							altseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
							altseq[r.variant.pos-1-refseq_start] = (byte) r.variant.alt.charAt(0);
						}

						if (distanceToIntronEnd > 0)
							index = "EX" + (intronIdx + 2) + "+" + distanceToIntronEnd;
						else
							index = "IVS" + (intronIdx + 1) + "" + (distanceToIntronEnd-1);
						
						//System.err.println(type + " " + pos + " " + new String(refseq) + " " + new String(altseq));
					}

				} else// reverse strand, TODO: rewrite by using relative coordinates
				{
					if (-ACCEPTOR_EXON <= distanceToIntronStart && distanceToIntronStart < ACCEPTOR_INTRON) {
						type = "acceptor";
						refseq_start = intron.start - ACCEPTOR_EXON;
						refseq_end = intron.start + ACCEPTOR_INTRON;
						refseq = Fasta.rc(Arrays.copyOfRange(chrSeq, refseq_start, refseq_end));
						if(r.variant.ref.length() > 1 || r.variant.alt.length() > 1)//ins or del overlap splice site
						{
							altseq = null;//assume that splice site is always disrupted in this case. We should search for nearby cryptic splice site in alternative sequence later.
						}else
						{
							//single nucleotide variant
							altseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
							altseq[r.variant.pos-1-refseq_start] = (byte) r.variant.alt.charAt(0);
							altseq = Fasta.rc(altseq);
						}

						if (distanceToIntronStart >= 0)
							index = "IVS" + (introns.size() - intronIdx) + "-" + (distanceToIntronStart+1);
						else
							index = "EX" + (introns.size() - intronIdx + 1) + "+" + (-distanceToIntronStart);
						
						//System.err.println(type + " " + pos + " " + new String(refseq) + " " + new String(altseq));
					}
					if (-DONOR_INTRON < distanceToIntronEnd && distanceToIntronEnd <= DONOR_EXON) {
						type = "donor";
						refseq_start = intron.end - DONOR_INTRON;
						refseq_end = intron.end + DONOR_EXON;
						refseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
						refseq = Fasta.rc(refseq);
						if(r.variant.ref.length() > 1 || r.variant.alt.length() > 1)//ins or del overlap splice site
						{
							altseq = null;//assume that splice site is always disrupted in this case. We should search for nearby cryptic splice site in alternative sequence later.
						}else
						{
							//single nucleotide variant
							altseq = Arrays.copyOfRange(chrSeq, refseq_start, refseq_end);
							altseq[r.variant.pos-1-refseq_start] = (byte) r.variant.alt.charAt(0);
							altseq = Fasta.rc(altseq);
						}

						if (distanceToIntronEnd > 0)
							index = "EX" + (introns.size() - intronIdx) + "-" + distanceToIntronEnd;
						else
							index = "IVS" + (introns.size() - intronIdx) + "+" + (1 - distanceToIntronEnd);
					}

				}

			}
		}
		
		//what if type == null ??? if this is not a splice variant according to autopvs1 criteria, but splice variant according to VEP. See clinvar:928280 for example.

		if (type.equals("donor")) {
			formatDonor(refseq);
			formatDonor(altseq);
		} else if (type.equals("acceptor")) {
			formatAcceptor(refseq);
			formatAcceptor(altseq);
		}

	}

	// Is this formatting correct when ref.length != alt.length???
	private void formatAcceptor(byte[] seq) {
		for (int i = 0; i < seq.length; i++) {
			byte c = seq[i];
			if (i >= ACCEPTOR_INTRON - 2 && i < ACCEPTOR_INTRON)
				c = (byte) Character.toUpperCase((char) c);
			else
				c = (byte) Character.toLowerCase((char) c);
			seq[i] = c;
		}
	}

	private void formatDonor(byte[] seq) {
		for (int i = 0; i < seq.length; i++) {
			byte c = seq[i];
			if (i >= DONOR_EXON && i < DONOR_EXON + 2)
				c = (byte) Character.toUpperCase((char) c);
			else
				c = (byte) Character.toLowerCase((char) c);
			seq[i] = c;
		}
	}

	public boolean isPreserveReadingFrame() {
		if (isExonSkipping())
			return getSkippedExonLength() % 3 == 0;
		
		if (isSpliceSiteDisrupted() && hasCrypticSpliceSite()) {
			return (getCrypticSpliceSite().pos - refseq_start) % 3 == 0;
		}
		
		//!spliceSiteDisrupted
		return true;
	}

	private boolean isExonSkipping() {
		//Why we assume that absence of splice site leads to exon skipping? what if intron retention???
		if (isSpliceSiteDisrupted() && !hasCrypticSpliceSite())
			return true;
		return false;
	}

	private boolean isSpliceSiteDisrupted() {
		return altseq == null || (maxentfoldchange < PERCENT_THRESHOLD 
				&& maxentscore_alt < 3);
	}

	static class SpliceSite {
		int pos;
		byte[] context;
		double maxentscore;

	}

	public SpliceSite getCrypticSpliceSite() {
		//Why we searach cryptic splice site in refseq, it should work in altseq???
		for (int d = 1; d <= 50; d++)
			for (int o : new int[] { -1, 1 }) {
				SpliceSite ss = new SpliceSite();
				ss.pos = refseq_start + o * d;// zero based
				if (type.equals("donor")) {
					ss.context = Arrays.copyOfRange(chrSeq, ss.pos, ss.pos + 9);
					int alt_index = r.variant.pos - ss.pos - 1;
					if (0 < alt_index && alt_index < 9) {
						// ???incorrect if alt.length != ref.length or alt.length > 1
						String spliceContextAlt = new String(ss.context, 0, alt_index) + r.variant.alt;
						int offset = alt_index + r.variant.alt.length();
						int length = 10 - 2 * r.variant.alt.length() - alt_index;
						if(offset < 0)
							offset = 0;
						if(offset >= ss.context.length)
							offset = ss.context.length-1;
						if(length > ss.context.length - offset)
							length = ss.context.length-offset;
						if(length < 0)
							length = 0;
						spliceContextAlt += new String(ss.context, offset, length);
						ss.context = spliceContextAlt.getBytes();
					}
					if(!r.transcript.forwardStrand)
						ss.context = Fasta.rc(ss.context);

					formatDonor(ss.context);
					if (ss.context.length == 9)
						ss.maxentscore = maxent.score5(ss.context);
					if(ss.context.length < 5)
						continue;
					String core = new String(ss.context, 3, 2);
					String refCore = new String(refseq, 3, 2);
					if ((core.equals("GT") || core.equals(refCore)) && ss.maxentscore > 1
							&& (ss.maxentscore > DONOR_THRESHOLD
									|| ss.maxentscore / maxentscore_ref >= PERCENT_THRESHOLD)) {
						return ss;
					}
				} else if (type.equals("acceptor")) {
					ss.context = Arrays.copyOfRange(chrSeq, ss.pos, ss.pos + 23);
					int alt_index = r.variant.pos - ss.pos - 1;
					if (0 < alt_index && alt_index < 23) {
						// ???incorrect if alt.length != ref.length or alt.length > 1
						String spliceContextAlt = new String(ss.context, 0, alt_index) + r.variant.alt;
						int offset = alt_index + r.variant.alt.length();
						int length = 24 - 2 * r.variant.alt.length() - alt_index;
						if(offset < 0)
							offset = 0;
						if(offset >= ss.context.length)
							offset = ss.context.length-1;
						if(length > ss.context.length - offset)
							length = ss.context.length-offset;
						if(length < 0)
							length = 0;
						spliceContextAlt += new String(ss.context, offset, length);
						ss.context = spliceContextAlt.getBytes();
					}
					if (!r.transcript.forwardStrand)
						ss.context = Fasta.rc(ss.context);
					formatAcceptor(ss.context);
					if (ss.context.length == 23)
						ss.maxentscore = maxent.score3(ss.context);
					if(ss.context.length < 20)
						continue;
					String core = new String(ss.context, 18, 2);
					String refCore = new String(refseq, 18, 2);
					if ((core.equals("AG") || core.equals(refCore)) && ss.maxentscore > 1
							&& (ss.maxentscore > ACCEPTOR_THRESHOLD
									|| ss.maxentscore / maxentscore_ref >= PERCENT_THRESHOLD)) {
						return ss;
					}
				}
			}
		return null;
	}

	public boolean hasCrypticSpliceSite() {
		if (type != null) {
			SpliceSite ss = getCrypticSpliceSite();
			if (ss != null) {
				int exonId = getSkippedExonId();
				exonId--;// make it zero based
				List<Interval> exons = getCrypticCodingExons();
				Interval crypticExon = exons.get(exonId);
				if (crypticExon != null && crypticExon.start < crypticExon.end)
					return true;
			}
		}
		return false;
	}

	Pattern pat1 = Pattern.compile("IVS(\\d+)([+|-])(\\d+)");
	Pattern pat2 = Pattern.compile("EX(\\d+)([+|-])(\\d+)");

	// ???I don't understand how exon skipping occurs
	public int getSkippedExonId() {
		// better to save these values during parsing, rather then parse index
		Matcher m1 = pat1.matcher(index);
		Matcher m2 = pat2.matcher(index);
		if (m1.matches()) {
			int intron_id = Integer.parseInt(m1.group(1));// one-based
			if (m1.group(2).charAt(0) == '+')
				return intron_id;
			else
				return intron_id + 1;
		}
		if (m2.matches()) {
			int exon_id = Integer.parseInt(m2.group(1));// one-based
			return exon_id;
		} else
			return 0;

	}

	private int getSkippedExonLength() {
		Matcher m1 = pat1.matcher(index);
		Matcher m2 = pat2.matcher(index);
		if (m1.matches()) {
			int intron_id = Integer.parseInt(m1.group(1));// one-based
			if (m1.group(2).charAt(0) == '+')
				return r.transcript.getCDSSizes().get(intron_id - 1);
			else
				return r.transcript.getCDSSizes().get(intron_id);// ???error: cdsSizes reversed if negative strand
		}
		if (m2.matches()) {
			int exon_id = Integer.parseInt(m2.group(1));// one-based
			return r.transcript.getCDSSizes().get(exon_id - 1);
		}
		return 0;
	}

	public List<Interval> getCrypticCodingExons() {
		List<Interval> codingExons = r.transcript.getCodingExons();
		int exonId = getSkippedExonId();// Skipped and changed exon are the same???
		exonId--;// make it zero based
		SpliceSite crSpliceSite = getCrypticSpliceSite();// what if crSpliceSite == null???
		if (codingExons.get(exonId) == null)// skipped exon is non-coding
			return codingExons;
		if ((r.transcript.forwardStrand && type.equals("donor"))
				|| (!r.transcript.forwardStrand && type.equals("acceptor")))
			codingExons.get(exonId).end = crSpliceSite.pos + 3;
		else if (type.equals("acceptor"))
			codingExons.get(exonId).start = crSpliceSite.pos + 20;
		else
			codingExons.get(exonId).start = crSpliceSite.pos + 6;
		return codingExons;
	}

	public boolean isCriticalToProteinFunction() {
		return isCriticalToProteinFunctionDetail().isFunctional;
	}

	private FunctionalRegion isCriticalToProteinFunctionDetail() {
		String chr = r.variant.chr;
		if (chr.startsWith("chr"))
			chr = chr.substring("chr".length());
		FunctionalRegion result = new FunctionalRegion();
		int start, end;
		if (hasCrypticSpliceSite() && isPreserveReadingFrame()) {
			int exonId = getSkippedExonId() - 1;
			SpliceSite crSpliceSite = getCrypticSpliceSite();
			if (r.transcript.forwardStrand) {
				if (type.equals("acceptor")) {
					start = r.transcript.exonStarts[exonId];
					end = crSpliceSite.pos;
				} else {
					start = crSpliceSite.pos;
					end = r.transcript.exonEnds[exonId];
				}
			} else {
				if (type.equals("acceptor")) {
					start = r.transcript.exonStarts[r.transcript.exonStarts.length - 1 - exonId];// ???error: exonStarts
																									// are not reversed
					end = crSpliceSite.pos;
				} else {
					start = crSpliceSite.pos;
					end = r.transcript.exonEnds[r.transcript.exonStarts.length - 1 - exonId];// ???error: exonEnds are
																								// not reversed
				}
			}
			if (start >= end)// we faild to interpret cryptic splice site
			{
				result.isFunctional = false;
				result.description = "NA";
				return result;
			}
		} else if (hasCrypticSpliceSite() && !isPreserveReadingFrame()) {
			// ???I don't understand this part
			// TODO: new stop codon position
			int exonId = getSkippedExonId() - 1;
			SpliceSite crSpliceSite = getCrypticSpliceSite();
			if (r.transcript.forwardStrand) {
				start = crSpliceSite.pos;
				end = r.transcript.cdsEnd;
			} else {
				start = r.transcript.cdsStart;
				end = crSpliceSite.pos;
			}
		} else if (isExonSkipping()) {
			int exonId = getSkippedExonId() - 1;
			if (r.transcript.forwardStrand) {
				start = r.transcript.exonStarts[exonId];
				end = r.transcript.exonEnds[exonId];
			} else {
				start = r.transcript.exonStarts[r.transcript.exonStarts.length - 1 - exonId];// ???error: exonStarts are
																								// not reversed
				end = r.transcript.exonEnds[r.transcript.exonEnds.length - 1 - exonId];// ???error: exonEnds are not
																						// reversed
			}

		} else {
			// ??? what is here? intron retention?
			result.isFunctional = false;
			result.description = "NA";
			return result;
		}

		return PVS1.getFunctionalRegion(chr, start, end, data);
	}

	public boolean variantRemoves10PercentOfProtein() {
		if (hasCrypticSpliceSite()) {
			int start, end;// ???should depend on donor/acceptor
			SpliceSite ss = getCrypticSpliceSite();
			if (r.transcript.forwardStrand) {
				start = refseq_start;
				end = ss.pos;
			} else {
				start = ss.pos;
				end = refseq_start;
			}
			return ((double) start - end) / r.transcript.getCDSLength() > 0.1;// ??? should be end-start
		} else if ((double) getSkippedExonLength() / r.transcript.getCDSLength() > 0.1)
			return true;
		else
			return false;
	}

	public boolean isUndergoNMD() {
		if (isPreserveReadingFrame())
			return false;
		if (hasCrypticSpliceSite())
			return getTransSeqInfo().isNMDTarget;
		else if (getSkippedExonId() >= r.transcript.exonStarts.length - 1)// last 2 exons
			return false;
		else
			return true;
	}

	static class TransSeqInfo {
		byte[] transSeq;
		int stopCodon;
		boolean isNMDTarget;
	}

	private TransSeqInfo getTransSeqInfo() {
		TransSeqInfo res = new TransSeqInfo();

		List<Integer> cdsSizes = new ArrayList<>();
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		for (Interval exon : getCrypticCodingExons()) {
			if (exon == null) {
				// non-coding exon
				cdsSizes.add(0);
			} else {
				cdsSizes.add(exon.end - exon.start);
				byte[] exonSeq = Arrays.copyOfRange(chrSeq, exon.start, exon.end);
				if (!r.transcript.forwardStrand)
					exonSeq = Fasta.rc(exonSeq);
				try {
					buffer.write(exonSeq);
				} catch (IOException e) {
					throw new AssertionError();
				}
			}
		}

		res.transSeq = buffer.toByteArray();
		res.stopCodon = 0;// ???what if stop codon not found?
		for (int x = 0; x + 2 < res.transSeq.length; x += 3) {
			String codon = new String(res.transSeq, x, 3).toUpperCase();
			if(codon.equals("TAA") || codon.equals("TAG") || codon.equals("TGA")) {
				res.stopCodon = x;
				break;
			}
		}

		if (cdsSizes.size() == 1 || cdsSizes.stream().filter(x -> x > 0).count() == 1)
			res.isNMDTarget = true;
		else {
			int nmdCutoff = 0;
			for (int i = 0; i < cdsSizes.size() - 1; i++)// except last exon
				nmdCutoff += cdsSizes.get(i);
			// subtract 50bp of penultimate exon (or the whole exon if shorter then 50bp)
			nmdCutoff -= Math.min(cdsSizes.get(cdsSizes.size() - 2), 50);

			res.isNMDTarget = res.stopCodon <= nmdCutoff;
		}
		return res;
	}

}
