package sciphy;

import beast.base.core.Log;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import sciphy.evolution.likelihood.SciPhyParsimony;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;
import org.junit.Test;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;

public class SciPhyParsimonyTest {

    @Test
    public void testAncestralSetsUneditedSequencesUneditedAncestorsParsimony() {
        // Testing the ancestral state reconstruction for a cherry
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;

        //first calculate states dictionary
        //manually create states:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));

        //parsimony
        double parsimonyExpected = 0;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);


    }



    @Test
    public void testAncestralSetsEditedSequencesEditedAncestorParsimony() {
        // Testing the ancestral state reconstruction for a cherry
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,3,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();

        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");
        RealParameter orig = new RealParameter("6");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength, "origin", orig);

        //test ancestral states sets calculations
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele123 = Arrays.asList(1, 2, 3, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        //the general dictionary
        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 4);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele123));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //parsimony
        double parsimonyExpected = 3;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);



    }



    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestorsParsimony() {
        // Testing the ancestral state reconstruction at internal nodes
        //CHERRY
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "2,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();

        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele21));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele2));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));

        //parsimony
        double parsimonyExpected = 4;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);


    }


    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestorsLengthsParsimony() {
        // Testing the ancestral state reconstruction at internal nodes
        //CHERRY
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "2,1,1,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();

        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele211 = Arrays.asList(2, 1, 1, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 4);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele211));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele21));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele2));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));
        //parsimony
        double parsimonyExpected = 5;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);

    }


    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestorsPositionsParsimony() {
        // Testing the ancestral state reconstruction at internal nodes
        //CHERRY
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "3,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();

        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;


        //first calculate states dictionary
        //manually create states:

        List<Integer> allele31 = Arrays.asList(3, 1, 0, 0, 0);
        List<Integer> allele3 = Arrays.asList(3, 1, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 2);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele31));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele3));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));

        //parsimony
        double parsimonyExpected = 3;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);

    }



    @Test
    public void testAncestralSets3LeavesParsimony() {

        // Testing the ancestral state reconstruction at internal nodes
        String newick = "((A:1,B:1):1,C:2)";
        Sequence a = new Sequence("A", "1,2,0,0,0");
        Sequence b = new Sequence("B", "1,1,0,0,0");
        Sequence c = new Sequence("C", "2,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");
        alignment.initByName("sequence", c, "dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();

        RealParameter editprobs = new RealParameter("0.5 0.5");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);


        //calculate states dictionary
        parsimony.calculateParsimony();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;
        assertEquals(5, statesDictionary.size());

        //Manually create states
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele11 = Arrays.asList(1, 1, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        //check the ancestral dictionaries
        //todo find a way to extract node numbers in a way that we know their position in the tree
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele11));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        // node c
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele21));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele2));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));

        //internal node between a and b
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(3)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(3)).contains(allele1));

        //root node a/b/c
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(4)).contains(allele0));

        //parsimony
        double parsimonyExpected = 5;
        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);

    }




    @Test
    public void testCherryIdenticalEditsParsimony() {


        //Testing the ancestral state reconstruction at internal nodes
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);


        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arraylength = new IntegerParameter("5");


        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);


        //initialise partialLikelihoods
        parsimony.partialParsimonies = new double[2][tree1.getNodeCount()][];

        //Manually calc the parsimony for that tree:

        //logparsimony
        double parsimonyExpected = 2;


        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);


    }



    @Test
    public void testCherry2Shared1DifferentInsertParsimony() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        //Using Origin Time: total height is 6


        Sequence a = new Sequence("CHILD1", "1,2,2,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        RealParameter editprobs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arraylength = new IntegerParameter("5");


        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);


        //initialise partialLikelihoods
        parsimony.partialParsimonies = new double[2][tree1.getNodeCount()][];

        //Manually calc the parsimony for that tree:

        //logparsimony
        double parsimonyExpected = 3;


        double parsimonyCalc = parsimony.calculateParsimony();

        Log.info.println(Arrays.deepToString(parsimony.partialParsimonies));

        assertEquals(parsimonyCalc, parsimonyExpected);


    }



    @Test
    public void testCherryNoInsertParsimony() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        RealParameter editprobs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);


        //initialise partialLikelihoods
        parsimony.partialParsimonies = new double[2][tree1.getNodeCount()][];

        //Manually calc the parsimony for that tree:
        double parsimonyExpected = 0;


        double parsimonyCalc = parsimony.calculateParsimony();

        assertEquals(parsimonyCalc, parsimonyExpected);


    }


    @Test
    public void test3LeavesSameLengthInsertsParsimony() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);
        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("4");
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);


        //initialise partialLikelihoods
        parsimony.partialParsimonies = new double[2][tree1.getNodeCount()][];

        //Manually calc the parsimony for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);

        //logparsimony
        double parsimonyCalc = parsimony.calculateParsimony();

        Hashtable<Integer, List<List<Integer>>> statesDictionary = parsimony.ancestralStates;

        //first calculate states dictionary
        //manually create states:
        int blue_sky = 0;
//        assertEquals(3, statesDictionary.size());
//
        //1st leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele12));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele1));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(1)).contains(allele0));

        //internal node leaf
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(3)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(3)).contains(allele12));


        // root
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(4)).contains(allele0));


        //3rd leaf node
        assertEquals(statesDictionary.get(parsimony.makeCachingIndexStates(2)).size(), 3);
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele21));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele2));
        assertTrue(statesDictionary.get(parsimony.makeCachingIndexStates(2)).contains(allele0));


        double parsimonyExpected = 4;


        assertEquals(parsimonyExpected, parsimonyCalc);


    }




    @Test
    public void test3LeavesSameLengthInsertsAllDifferentParsimony() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,1,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        SciPhyParsimony parsimony = new SciPhyParsimony();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //parsimony class

        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("4");
        IntegerParameter arrayLength = new IntegerParameter("5");

        parsimony.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);

        //logparsimony
        double parsimonyExpected = 5;

        double parsimonyCalc = parsimony.calculateParsimony();
        assertEquals(parsimonyExpected, parsimonyCalc);

    }

}