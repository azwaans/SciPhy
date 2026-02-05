package sciphy;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import sciphy.evolution.likelihood.SciPhyTreeLikelihood;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;

public class SciPhyLikelihoodTest {

    SciPhyTreeLikelihood sciphyLikelihood;
    SciPhySubstitutionModel substModel;
    SiteModel siteM;
    StrictClockModel clockModel;
    RealParameter origin;
    IntegerParameter arraylength;

    // TODO test for input of an empty arrayLength 

    @Before
    public void setUp() {
        sciphyLikelihood = new SciPhyTreeLikelihood();
        substModel = new SciPhySubstitutionModel();
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0.0 0.0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model

        siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);


        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);

        origin = new RealParameter("6");
        arraylength = new IntegerParameter("5");

    }

    @Test(expected = RuntimeException.class)
    public void testExceptionForSingleBranchInput() {

        String newick = "(CHILD1:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        // this should fail
        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);
    }

    @Test(expected = RuntimeException.class)
    public void testExceptionForSingleNodeInput() {

        String newick = "(CHILD1:0)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);
    }

    @Test(expected = RuntimeException.class)
    public void testExceptionForOriginNegativeInput() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        RealParameter originNegative = new RealParameter("-1");
        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", originNegative, "arrayLength", arraylength);

    }

    @Test
    public void testMinusInfinityForOriginSmallerThanTreeHeightInput() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        RealParameter originNegative = new RealParameter("4");
        // this should fail
        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", originNegative, "arrayLength", arraylength);
        double LogPCalc = sciphyLikelihood.calculateLogP();
        assertEquals(Double.NEGATIVE_INFINITY, LogPCalc);


    }

    @Test(expected = RuntimeException.class)
    public void testExceptionForNegativeTargetLengthInput() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        IntegerParameter arrayNegative = new IntegerParameter("-1");
        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayNegative);

    }

    //TODO maybe automatically set the array length based on input data?
    @Test(expected = RuntimeException.class)
    public void testExceptionForMismatchedTargetLengthInput() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        IntegerParameter arrayTooLarge = new IntegerParameter("6");

        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayTooLarge);

    }


    @Test
    public void testGetPossibleAncestorsSetsEdited() {
        Sequence a = new Sequence("cell1", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = SciPhyTreeLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 3, 1e-5);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

        assertTrue(ancs_sequence_a.contains(allele12));
        assertTrue(ancs_sequence_a.contains(allele1));
        assertTrue(ancs_sequence_a.contains(allele0));

    }

    @Test
    public void testGetPossibleAncestorsSetsUnedited() {
        Sequence a = new Sequence("cell1", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = SciPhyTreeLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 1);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited
        assertTrue(ancs_sequence_a.contains(allele0));

    }

    @Test
    public void testAncestralSetsUneditedSequencesUneditedAncestors() {
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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

        //create a sub model with values
        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        RealParameter editprobs = new RealParameter("1.0 0 0 0");
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;

        //first calculate states dictionary
        //manually create states:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));

    }

    @Test
    public void testAncestralSetsEditedSequencesEditedAncestors() {
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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele123 = Arrays.asList(1, 2, 3, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 4);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele123));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));


    }

    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestors() {
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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele21));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele2));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));


    }

    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestorsLengths() {
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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


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
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 4);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele211));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele21));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele2));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));


    }

    @Test
    public void testAncestralSetsEditedSequencesUneditedAncestorsPositions() {
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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:

        List<Integer> allele31 = Arrays.asList(3, 1, 0, 0, 0);
        List<Integer> allele3 = Arrays.asList(3, 1, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 2);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));


        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele31));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele3));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));


    }


    @Test
    public void testAncestralSets3Leaves() {

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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);


        //calculate states dictionary
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;
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
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele12));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele11));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        // node c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 3);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele21));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele2));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));

        //internal node between a and b
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(3)).size(), 2);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(3)).contains(allele1));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(3)).contains(allele0));

        //root node a/b/c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(4)).contains(allele0));


    }


    @Test
    public void testLikelihoodCherryIdenticalEdits() {


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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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


        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arraylength = new IntegerParameter("5");

     
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);
        

        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        double clockRate = 0.5;

        double partial0000Internal = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 5 * clockRate, 5);
        double partial1000Internal = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele12, 5 * clockRate, 5);
        double partial1200Internal = substitutionModel.getSequenceTransitionProbability(allele12, allele12, 5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele12, allele12, 5 * clockRate, 5);

        //root node
        double partialOrigin = partial0000Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5) + partial1200Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }

    @Test
    public void testLikelihoodCherry2Shared1DifferentInsert() {


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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arraylength = new IntegerParameter("5");


        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele122 = Arrays.asList(1, 2, 2, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        double clockRate = 0.5;

        double partial0000Internal = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 5 * clockRate,5) * substitutionModel.getSequenceTransitionProbability(allele0, allele122, 5 * clockRate, 5);
        double partial1000Internal = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele122, 5 * clockRate, 5);
        double partial1200Internal = substitutionModel.getSequenceTransitionProbability(allele12, allele12, 5 * clockRate,5) * substitutionModel.getSequenceTransitionProbability(allele12, allele122, 5 * clockRate, 5);

        //root node
        double partialOrigin = partial0000Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5) + partial1200Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }

    @Test
    public void testLikelihoodCherryNoInsert() {


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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("6");
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        double clockRate = 0.5;

        double partial0000Internal = substitutionModel.getSequenceTransitionProbability(allele0, allele0, 5 * clockRate,5) * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 5 * clockRate,5);

        //root node
        double partialOrigin = partial0000Internal * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate,5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }


    @Test
    public void testLikelihood3LeavesSameLengthInserts() {


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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("4");
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        double clockRate = 0.5;


        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5);
        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5);

        double partial1200Internal1 = substitutionModel.getSequenceTransitionProbability(allele12, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele12, allele12, 1 * clockRate, 5);

        double partial0000Internal2 = (partial0000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5) + partial1200Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)) * (substitutionModel.getSequenceTransitionProbability(allele0, allele21, 2 * clockRate, 5));

        //root node
        double partialOrigin = partial0000Internal2 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2 * clockRate, 5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected, LogPCalc);


    }


    @Test
    public void testLikelihood3LeavesSameLengthInsertsAllDifferent() {


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

        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();

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

        //likelihood class

        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("4");
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele11 = Arrays.asList(1, 1, 0, 0, 0);
        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        double clockRate = 0.5;

        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele0, allele11, 1 * clockRate, 5);
        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele11, 1 * clockRate, 5);

        double partial0000Internal2 = (partial0000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)) * (substitutionModel.getSequenceTransitionProbability(allele0, allele21, 2 * clockRate, 5));


        //root node
        double partialOrigin = partial0000Internal2 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2 * clockRate, 5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected, LogPCalc);

    }

}