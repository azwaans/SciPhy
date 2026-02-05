package sciphy.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType.Base;

@Description("Datatype for Sciphy Data")
public class SciPhyData extends Base {

    public SciPhyData() {
        stateCount = 14;
        //for now, we use a integer representation of a Sciphy sequence.
        //todo change this mapping to actual nucleotide insert
        codeLength = 2; // with length 3 nucleotide inserts input: codeLength = 3
        codeMap = "0001020304050607080910111213" ; // with length 3 nucleotide inserts input codeMap = ATGCCTTAGTT etc
        //with both codeLength 2 and 3 the resulting representation should be an integer array.
        mapCodeToStateSet = new int[14][];
        for (int i = 0; i < 14; i++) {
            mapCodeToStateSet[i] = new int[1];
            mapCodeToStateSet[i][0] = i;
        }


    }

    public void initAndValidate() {
        //TODO split the string into list of strings
        stateCount = -1;
        mapCodeToStateSet = null;
        codeLength = -1;
        codeMap = null;
    }

    @Override
    public String getTypeDescription() {
        return "SciphyData";
    }
    
    @Override
    public boolean isAmbiguousCode(int code) {
    	return code < 0;
    }

    @Override
    public String getCharacter(int code) {
    	if (code < 0) {
    		return "?";
    	}
    	return code + "";
    }
    
	@Override
	public int[] getStatesForCode(int code) {
		return new int[]{code};
	}
}

