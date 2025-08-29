package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class Fasta {
	public static byte[] subSeq(byte[] x, int from, int length, boolean rc, byte stub)
    {
    	byte[] result = new byte[length];
    	for(int i = from; i < from + length; i++)
    		result[i-from] = i < 0 || i >= x.length ? stub : x[i];
    	if(rc)
    		result = rc(result);
    	return result;
    }
    
    public static byte[] rc(byte[] seq)
    {
    	byte[] rc = new byte[seq.length];
    	for(int i = 0; i < seq.length; i++)
    	{
    		rc[rc.length - i - 1] = comp(seq[i]);
    	}
    	return rc;
    }

	public static byte comp(byte b) {
		switch(b){
		case 'A': return 'T';
		case 'a': return 't';
		case 'C': return 'G';
		case 'c': return 'g';
		case 'G': return 'C';
		case 'g': return 'c';
		case 'T': return 'A';
		case 't': return 'a';
		default: return b;
		}
	}
	
	public static byte code(byte letter)
	{
		switch(letter)
		{
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		default: return 4;
		}
	}

	public static Map<String, byte[]> loadFasta(File fastaFile) throws IOException
    {
    	Map<String, byte[]> result = new HashMap<>();
    	InputStream is = new FileInputStream(fastaFile);
    	if(fastaFile.getName().endsWith(".gz"))
    		is = new GZIPInputStream(is);
    	try(BufferedReader reader = new BufferedReader(new InputStreamReader(is)))
    	{
    		String line = reader.readLine();
    		while(line != null)
    		{
    			if(!line.startsWith(">"))
    				throw new IOException("Expecting >");
    			String name = line.substring(1);
    			int idx = name.indexOf(' ');
    			if(idx != -1)
    				name = name.substring(0, idx);
    			ByteArrayOutputStream data = new ByteArrayOutputStream();
    			while((line = reader.readLine()) != null)
    			{
    				if(line.startsWith(">"))
    					break;
    				data.write(line.getBytes());
    			}
    			result.put(name, data.toByteArray());
    		}
    	}
    	return result;
    }
}
