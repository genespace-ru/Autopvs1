package ru.biosoft.autopvs1;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ru.biosoft.rtree.SiteIndex;

public class BedFile {
	static class Item
	{
		String chrom;
		int start, end;
		String name;
	}
	public List<Item> items = new ArrayList<>();
	public Map<String, Item> byName = new HashMap<>();
	public SiteIndex<Item> index = new SiteIndex<>();
	
	
	public static BedFile load(File file) throws IOException
	{
		BedFile res = new BedFile();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith("#"))
				continue;
			String[] parts = line.split("\t");
			String chrom = parts[0];
			int start = Integer.parseInt(parts[1]);
			int end = Integer.parseInt(parts[2]);
			if(parts.length >= 12)
			{
				int blockCount = Integer.parseInt(parts[9]);
				String[] blockSizes = parts[10].split(",");
				String[] blockStarts = parts[11].split(",");
				for(int i = 0; i < blockCount; i++)
				{
					int itemStart = start + Integer.parseInt(blockStarts[i]);
					int itemEnd = itemStart + Integer.parseInt(blockSizes[i]);
					
					Item item = new Item();
					item.name = parts[3] + "|" + itemStart + "-" + itemEnd;
					item.chrom = chrom;
					item.start = itemStart;
					item.end = itemEnd;
					res.addItem(item);
				}
			}else
			{
				Item item = new Item();
				item.chrom = chrom;
				item.start = start;
				item.end = end;
				if(parts.length > 3)
					item.name = parts[3];
				else
					item.name = chrom+":"+start+"-"+end;
				res.addItem(item);
			}
			
					
		}
		reader.close();
		res.index.build();
		return res;
	}
	
	private void addItem(Item item)
	{
		items.add(item);
		byName.put(item.name, item);
		index.addSite(item, item.chrom, item.start, item.end-1);
	}
}
