package ru.biosoft.autopvs1;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TestMaxent {

	@Test
	public void testScore3() throws Exception {
		Maxent maxent = new Maxent();
		
		double score = maxent.score3("ttccaaacgaacttttgtAGgga".getBytes());
		assertEquals(2.89, score, 1e-2);
		
		score = maxent.score3("ctctactactatctatctagatc".getBytes());
		assertEquals(6.71, score, 1e-2);
		
		score = maxent.score3("tgtctttttctgtgtggcAGtgg".getBytes());
		assertEquals(8.19, score, 1e-2);
		
		score = maxent.score3("ttctctcttcagacttatAGcaa".getBytes());
		assertEquals(-0.08, score, 1e-2);
	}
	
	@Test
	public void testScore5() throws Exception {
		Maxent maxent = new Maxent();
		
		double score = maxent.score5("acggtaagt".getBytes());
		assertEquals(11.81, score, 1e-2);
		
		score = maxent.score5("cagGTAAGT".getBytes());
		assertEquals(10.86, score, 1e-2);
		
		score = maxent.score5("gagGTAAGT".getBytes());
		assertEquals(11.08, score, 1e-2);
		
		score = maxent.score5("taaATAAGT".getBytes());
		assertEquals(-0.12, score, 1e-2);
	}
}
