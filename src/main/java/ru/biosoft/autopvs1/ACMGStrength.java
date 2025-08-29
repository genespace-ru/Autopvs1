package ru.biosoft.autopvs1;

public enum ACMGStrength {
	UNSET, UNMET, SUPPORTING, MODERATE, STRONG, VERY_STRONG;
	
	public ACMGStrength add(int count)
	{
		int newOrdinal = this.ordinal() + count;
		if(newOrdinal < 1)//down to UNMET
			newOrdinal = 1;
		if(newOrdinal >= values().length)
			newOrdinal = values().length-1;
		return values()[newOrdinal];
	}
	
	public ACMGStrength upgrade(int count)
	{
		return add(count);
	}
	public ACMGStrength downgrade(int count)
	{
		return add(-count);
	}
}
