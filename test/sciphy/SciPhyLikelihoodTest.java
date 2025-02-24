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
        sciphyLikelihood = new SciPhyTreeLikelihood();

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
        List<List<Integer>> ancs_sequence_a;
        ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

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
    public void testGetPossibleAncestorsSetsBothMissingProcesses() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b, "dataType", "integer");
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0.0 0.0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);

        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.1");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);



        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);



        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 2, 1e-5);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

        assertTrue(ancs_sequence_a.contains(alleleLost));
        assertTrue(ancs_sequence_a.contains(alleleWC));

    }

    @Test
    public void testGetPossibleAncestorsSetsOnlyDropout() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b, "dataType", "integer");
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0.0 0.0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);

        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        //missingness parameters
        RealParameter missRate = new RealParameter("0.0");
        RealParameter missProbability = new RealParameter("0.1");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);



        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);



        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 1, 1e-5);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

        assertTrue(ancs_sequence_a.contains(alleleWC));

    }

    @Test
    public void testGetPossibleAncestorsSetsOnlyLoss() {
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b, "dataType", "integer");
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);
        RealParameter editprobs = new RealParameter("0.8 0.2");
        RealParameter stateFrequencies = new RealParameter("1.0 0.0 0.0 ");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);

        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
        //missingness parameters
        RealParameter missRate = new RealParameter("0.1");
        RealParameter missProbability = new RealParameter("0.0");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);



        sciphyLikelihood.initByName("data", alignment, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);



        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 1, 1e-5);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

        assertTrue(ancs_sequence_a.contains(alleleLost));

    }

    @Test(expected = RuntimeException.class)
    public void testGetPossibleAncestorsMissingState() {
        //dummy tree and alignment
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence b = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignmentDummy = new Alignment();
        alignmentDummy.initByName("sequence", b, "dataType", "integer");
        alignmentDummy.initByName("sequence", c, "dataType", "integer");
        //this won't be used for this test

        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignmentDummy, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        sciphyLikelihood.initByName("data", alignmentDummy, "tree", tree, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);



        Sequence a = new Sequence("cell1", "?,?,?,?,?");


        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

        //manually create ancestral states
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
        assertEquals(ancs_sequence_a.size(), 2, 1e-5);
        //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

        assertTrue(ancs_sequence_a.contains(alleleWC));
        assertTrue(ancs_sequence_a.contains(alleleLost));



    }

    @Test
    public void testGetPossibleAncestorsSetsUnedited() {
        Sequence a = new Sequence("cell1", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        //ancestral sequences
        List<List<Integer>> ancs_sequence_a = sciphyLikelihood.getPossibleAncestors(sequence_a);

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
    public void testAncestralSetsMissingUneditedAncestors() {
        // Testing the ancestral state reconstruction for a cherry
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
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


        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.1");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];
//
//        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;

        //first calculate states dictionary
        //manually create states:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 2);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleWC));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleLost));


        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));
        assertEquals(1,statesDictionary.get(likelihood.makeCachingIndexStates(2)).size());

    }

    @Test
    public void testAncestralSetsMissingUneditedAncestorsOnlyRate() {
        // Testing the ancestral state reconstruction for a cherry
        String newick = "(CHILD1:1,CHILD2:1)";
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
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
        RealParameter missRate = new RealParameter("1.0");
        RealParameter missProb = new RealParameter("0.0");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);

        //likelihood class
        RealParameter meanRate = new RealParameter("1.0");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        IntegerParameter arrayLength = new IntegerParameter("5.0");

        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "arrayLength", arrayLength);

        //test ancestral states sets calculations
        likelihood.calculateLogP();
        Hashtable<Integer, List<List<Integer>>> statesDictionary = likelihood.ancestralStates;

        //first calculate states dictionary
        //manually create states:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleLost));


        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));

    }

    @Test
    public void testAncestralSetsMissingUneditedAncestorsOppositeOrder() {
        // Testing the ancestral state reconstruction for a cherry
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "?,?,?,?,?");

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
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProb = new RealParameter("0.5");
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingProbability",missProb,"missingRate",missRate);        //site model
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
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 2);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(alleleWC));
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(alleleLost));


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
    public void testAncestralSets3LeavesMISSINGRATE1() {

        // Testing the ancestral state reconstruction at internal nodes
        String newick = "((A:1,B:1):1,C:2)";
        Sequence a = new Sequence("A", "?,?,?,?,?");
        Sequence b = new Sequence("B", "0,0,0,0,0");
        Sequence c = new Sequence("C", "0,0,0,0,0");

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
        RealParameter missProb = new RealParameter("0.0");
        RealParameter missRate = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, -0);


        //check the ancestral dictionaries
        //todo find a way to extract node numbers in a way that we know their position in the tree
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleLost));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(allele0));

        // node c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));

        //internal node between a and b
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(3)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(3)).contains(allele0));

        //root node a/b/c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(4)).contains(allele0));


    }

    @Test
    public void testAncestralSets3LeavesMISSINGRATE2() {

        // Testing the ancestral state reconstruction at internal nodes
        String newick = "((A:1,B:1):1,C:2)";
        Sequence a = new Sequence("A", "?,?,?,?,?");
        Sequence b = new Sequence("B", "?,?,?,?,?");
        Sequence c = new Sequence("C", "0,0,0,0,0");

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
        RealParameter missProb = new RealParameter("0.0");
        RealParameter missRate = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, -0);


        //check the ancestral dictionaries
        //todo find a way to extract node numbers in a way that we know their position in the tree
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleLost));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(alleleLost));

        // node c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(allele0));

        //internal node between a and b
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(3)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(3)).contains(alleleLost));

        //root node a/b/c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(4)).contains(allele0));


    }

    @Test
    public void testAncestralSets3LeavesMISSINGRATE3() {

        // Testing the ancestral state reconstruction at internal nodes
        String newick = "((A:1,B:1):1,C:2)";
        Sequence a = new Sequence("A", "?,?,?,?,?");
        Sequence b = new Sequence("B", "?,?,?,?,?");
        Sequence c = new Sequence("C", "?,?,?,?,?");

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
        RealParameter missProb = new RealParameter("0.0");
        RealParameter missRate = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, -0);


        //check the ancestral dictionaries
        //todo find a way to extract node numbers in a way that we know their position in the tree
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(0)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(0)).contains(alleleLost));

        //2nd leaf
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(1)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(1)).contains(alleleLost));

        // node c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(2)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(2)).contains(alleleLost));

        //internal node between a and b
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(3)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(3)).contains(alleleLost));

        //root node a/b/c
        assertEquals(statesDictionary.get(likelihood.makeCachingIndexStates(4)).size(), 1);
        assertTrue(statesDictionary.get(likelihood.makeCachingIndexStates(4)).contains(alleleLost));


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
    public void testLikelihoodCherryOneMissingOnlyRate() {


        //Testing the ancestral state reconstruction at internal nodes
        String newick = "(CHILD1:1,CHILD2:1)";
        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "?,?,?,?,?");

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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.0");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);


        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);


        //likelihood class
        RealParameter meanRate = new RealParameter("1.0");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("2");
        IntegerParameter arraylength = new IntegerParameter("5");


        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];



        double LogPCalc = likelihood.calculateLogP();
        double calculatedOnR = -3.932752;
        assertEquals(LogPCalc, calculatedOnR,0.00001);


    }

    @Test
    public void testLikelihoodCherryOneMissingOnlyRateOnlyLost() {


        //Testing the ancestral state reconstruction at internal nodes
        String newick = "(CHILD1:1,CHILD2:1)";
        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Sequence b = new Sequence("CHILD2", "?,?,?,?,?");

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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.0");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);


        //likelihood class
        RealParameter meanRate = new RealParameter("1.0");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", meanRate);
        RealParameter origin = new RealParameter("2");
        IntegerParameter arraylength = new IntegerParameter("5");


        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arraylength);


        //initialise partialLikelihoods
        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];



        double LogPCalc = likelihood.calculateLogP();
        double calculatedOnR = -0.9327521;


        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);
        List<Integer> allele00 = Arrays.asList(0, 0, 0, 0, 0);

        double partial00000Internal1 = substitutionModel.getSequenceTransitionProbability(allele00, alleleLost, 1, 5) ;
        assertEquals(Math.log(partial00000Internal1), calculatedOnR,0.00001);


        assertEquals(LogPCalc, calculatedOnR,0.00001);



    }


    @Test
    public void TestLikelihoodCherryIdenticalEditsWithAMissing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";

        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");
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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.5");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);

        double clockRate = 0.5;


        //internal1 is the node connecting Child1 and Child3
        double partial12000Internal1Left = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 1 * clockRate, 5);
        double partial12000Internal1Right = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, alleleWC, 1 * clockRate, 5);

        double partial12000Internal1 = partial12000Internal1Left * partial12000Internal1Right;


        double partial00000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 1 * clockRate, 5);
        double partial10000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleWC, 1 * clockRate, 5);

        double partial00000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                                        substitutionModel.getSequenceTransitionProbability(allele0, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                                        substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                                        (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 5 * clockRate, 5));

        double partial10000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele1, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 5 * clockRate, 5));

        double partial12000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele12, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 5 * clockRate, 5));


        double partial00000Root =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Internal2);

        double LogPExpected = Math.log(partial00000Root);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void TestLikelihoodCherryIdenticalEditsWithAMissingNoMissRate() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";

        Sequence a = new Sequence("CHILD1", "1,0,0,0,0");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,0,0,0,0");
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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.0");
        RealParameter missProbability = new RealParameter("0.5");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);

        double clockRate = 0.5;
        double missingRate = 0.0;

        //internal1 is the node connecting Child1 and Child3
        double partial00000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele1, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 1 * clockRate, 5);
        double partial10000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele1, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleWC, 1 * clockRate, 5);

        double partial00000Internal2 = (substitutionModel.getSequenceTransitionProbability(allele0, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele1, 5 * clockRate, 5));

        double partial10000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele1, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele1, 5 * clockRate, 5));


        double partial00000Root = (substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Internal2);

        double LogPExpected = Math.log(partial00000Root);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }


    @Test
    public void TestLikelihoodCherryIdenticalEditsWithAMissingNoMissProb() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";

        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");
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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.0");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        double clockRate = 0.5;

        //internal1 is the node connecting Child1 and Child3
        double partial12000Internal1Left = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 1 * clockRate, 5);
        double partial12000Internal1Right = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, alleleLost, 1 * clockRate, 5);

        double partial12000Internal1 = partial12000Internal1Left * partial12000Internal1Right;

        double partial00000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleLost, 1 * clockRate, 5);
        double partial10000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleLost, 1 * clockRate, 5);

        double partial00000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 5 * clockRate, 5));

        double partial10000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele1, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 5 * clockRate, 5));

        double partial12000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele12, allele12, 4 * clockRate, 5)* partial12000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele1, 4 * clockRate, 5)* partial10000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele0, 4 * clockRate, 5)* partial00000Internal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 5 * clockRate, 5));


        double partial00000Root =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Internal2);

        double LogPExpected = Math.log(partial00000Root);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void TestLikelihoodCherryIdenticalEditsWith2Missing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";

        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");
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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProbability = new RealParameter("0.5");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        double clockRate = 0.5;
        double missingRate = 0.5;
        double missingProbability = 0.5;

        //internal1 is the node connecting Child1 and Child3
        double partialWCInternal1Left = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);
        double partialWCInternal1Right = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);

        double manualPartialWCInternal1 = ((1- Math.exp(- missingRate * clockRate)) + Math.exp(- missingRate * clockRate) * missingProbability) * ((1- Math.exp(- missingRate * clockRate)) + Math.exp(- missingRate * clockRate) * missingProbability);
        double partialmanualWCInternal1 = partialWCInternal1Left * partialWCInternal1Right;
        assertEquals(manualPartialWCInternal1, partialmanualWCInternal1);


        double partialWCInternal1 =  (substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5)) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5));

        double partialLostInternal1 =  (substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleLost, alleleWC, 1 * clockRate, 5)) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleLost, alleleWC, 1 * clockRate, 5));



        double partial00000Root =  (substitutionModel.getSequenceTransitionProbability(allele0, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, alleleLost, 4 * clockRate, 5)* partialLostInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 5 * clockRate, 5));

        double partial10000Root =  (substitutionModel.getSequenceTransitionProbability(allele1, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, alleleLost, 4 * clockRate, 5)* partialLostInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 5 * clockRate, 5));

        double partial12000Root =  (substitutionModel.getSequenceTransitionProbability(allele12, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, alleleLost, 4 * clockRate, 5)* partialLostInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 5 * clockRate, 5));


        double partial00000Origin =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)* partial12000Root +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Root +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Root);

        double LogPExpected = Math.log(partial00000Origin);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void TestLikelihoodCherryIdenticalEditsWith2MissingNoRate() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";

        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");
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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.0");
        RealParameter missProbability = new RealParameter("0.5");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);

        double clockRate = 0.5;

        //internal1 is the node connecting Child1 and Child3
        double partialWCInternal1Left = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);
        double partialWCInternal1Right = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);

        double partialWCInternal1 = partialWCInternal1Left * partialWCInternal1Right;


        double partial00000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele0, alleleWC, 4 * clockRate, 5)* partialWCInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 5 * clockRate, 5));

        double partial10000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele1, alleleWC, 4 * clockRate, 5)* partialWCInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 5 * clockRate, 5));

        double partial12000Internal2 =  (substitutionModel.getSequenceTransitionProbability(allele12, alleleWC, 4 * clockRate, 5)* partialWCInternal1) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 5 * clockRate, 5));


        double partial00000Root =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Internal2);

        double LogPExpected = Math.log(partial00000Root);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void TestLikelihoodCherryIdenticalEditsWith2MissingNoRate4tips() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,(CHILD2:1,CHILD4:1):4.0)";

        Sequence a = new Sequence("CHILD1", "?,?,?,?,?");
        Sequence b = new Sequence("CHILD3", "?,?,?,?,?");
        Sequence c = new Sequence("CHILD2", "1,2,0,0,0");
        Sequence d = new Sequence("CHILD4", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "dataType", "integer");


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

        //missingness parameters
        RealParameter missRate = new RealParameter("0.0");
        RealParameter missProbability = new RealParameter("0.5");

        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
                "missingRate",missRate,"missingProbability",missProbability);
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

        ///Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        List<Integer> alleleLost = Arrays.asList(-2, -2, -2, -2, -2);

        double clockRate = 0.5;
        double missingRate = 0.0;
        double missingProbability = 0.5;

        //Manual checks for the transition probability calculations
        double manualPartialWCInternal1 = ((1-Math.exp(-missingRate * clockRate)) + (Math.exp(-missingRate * clockRate)*missingProbability)) * (((1-Math.exp(-missingRate * clockRate)) + (Math.exp(-missingRate * clockRate)*missingProbability))) ;
        double manualpartialWCInternal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);
        assertEquals(manualPartialWCInternal1, manualpartialWCInternal1);

        //Manual likelihood calculation relying on the sequenceTransitionProbability
        double partialWCInternal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleWC, alleleWC, 1 * clockRate, 5);
        double partialLostInternal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleLost, alleleWC, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(alleleLost, alleleWC, 1 * clockRate, 5);


        double partial00000Internal2 =  (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1 * clockRate, 5))*
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1 * clockRate, 5));

        double partial10000Internal2 =  (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1 * clockRate, 5))*
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1 * clockRate, 5));

        double partial12000Internal2 =  (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 1 * clockRate, 5))*
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 1 * clockRate, 5));



        double partial00000Root =  (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 4 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 4 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate, 5)* partial00000Internal2) *
                (substitutionModel.getSequenceTransitionProbability(allele0, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, alleleLost, 4 * clockRate, 5)* partialLostInternal1);

        double partial10000Root =  (substitutionModel.getSequenceTransitionProbability(allele1, allele12, 4 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele1, 4 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele0, 4 * clockRate, 5)* partial00000Internal2)  *
                (substitutionModel.getSequenceTransitionProbability(allele1, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                        substitutionModel.getSequenceTransitionProbability(allele1, alleleLost, 4 * clockRate, 5)* partialLostInternal1);

        double partial12000Root =  (substitutionModel.getSequenceTransitionProbability(allele12, allele12, 4 * clockRate, 5)* partial12000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele1, 4 * clockRate, 5)* partial10000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele0, 4 * clockRate, 5)* partial00000Internal2)  *
                (substitutionModel.getSequenceTransitionProbability(allele12, alleleWC, 4 * clockRate, 5)* partialWCInternal1 +
                        substitutionModel.getSequenceTransitionProbability(allele12, alleleLost, 4 * clockRate, 5)* partialLostInternal1);

        double partial00000Origin =  (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5)* partial00000Root) +
                (substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)* partial10000Root) +
                        (substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)* partial12000Root)  ;

        double LogPExpected = Math.log(partial00000Origin);
        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


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
    public void testLikelihoodCherry2Shared1DifferentInsertMissing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";
        //Using Origin Time: total height is 6


        Sequence a = new Sequence("CHILD1", "1,2,2,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD3", "?,?,?,?,?");


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
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        RealParameter editprobs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);

        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProb = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        double clockRate = 0.5;

        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele122, 1 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 1 * clockRate, 5);
        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele122, 1* clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleWC, 1 * clockRate, 5);
        double partial1200Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele122, 1 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, alleleWC, 1 * clockRate, 5);
        double partial1220Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele122, allele122, 1 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele122, alleleWC, 1 * clockRate, 5);

        //root node
        double partial0000Root = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate,5) *  partial0000Internal1 +
                                 substitutionModel.getSequenceTransitionProbability(allele0, allele1, 4 * clockRate,5) *  partial1000Internal1 +
                                 substitutionModel.getSequenceTransitionProbability(allele0, allele12, 4 * clockRate,5) *  partial1200Internal1 +
                                 substitutionModel.getSequenceTransitionProbability(allele0, allele122, 4 * clockRate,5) *  partial1220Internal1) *
                substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 5 * clockRate,5);

        double partial1000Root = (substitutionModel.getSequenceTransitionProbability(allele1, allele0, 4 * clockRate,5) *  partial0000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele1, 4 * clockRate,5) *  partial1000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele12, 4 * clockRate,5) *  partial1200Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele122, 4 * clockRate,5) *  partial1220Internal1) *
                substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 5 * clockRate,5);

        double partial1200Root = (substitutionModel.getSequenceTransitionProbability(allele12, allele0, 4 * clockRate,5) *  partial0000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele1, 4 * clockRate,5) *  partial1000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele12, 4 * clockRate,5) *  partial1200Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele122, 4 * clockRate,5) *  partial1220Internal1) *
                substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 5 * clockRate,5);


        double partialOrigin = substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate,5) *  partial0000Root +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate,5) *  partial1000Root +
                substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate,5) *  partial1200Root ;


        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);


        double LogPCalc = likelihood.calculateLogP();
        Log.info.println("Log_lik:" + LogPCalc);

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
    public void testLikelihoodCherryNoInsertMissing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:4,CHILD2:5.0)";
        //Using Origin Time: total height is 6


        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "0,0,0,0,0");
        Sequence c = new Sequence("CHILD3", "?,?,?,?,?");


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
        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
        RealParameter editprobs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProb = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleWC = Arrays.asList(-1, -1, -1, -1, -1);
        double clockRate = 0.5;

        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele0, 1 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 1 * clockRate,5);

        //root node
        double partial0000Internal2 = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 4 * clockRate,5) * partial0000Internal1 )* substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele0, 5 * clockRate,5);

        //loglikelihood
        double partial0000Origin = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate,5) * partial0000Internal2 );

        double LogPExpected = Math.log(partial0000Origin);


        double LogPCalc = likelihood.calculateLogP();
        Log.info.println("Log_lik:" + LogPCalc);

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
    public void testLikelihood3LeavesSameLengthInsertsMissing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "(((CHILD1:0.5,MISS:0.5)INTERNAL1:0.5,CHILD3:1)INTERNAL2:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Sequence d = new Sequence("MISS", "?,?,?,?,?");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "dataType", "integer");


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


        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProb = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);
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
        List<Integer> alleleWC = Arrays.asList(-1,-1,-1,-1,-1);
        double clockRate = 0.5;


        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 0.5 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 0.5 * clockRate, 5);
        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 0.5* clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleWC, 0.5 * clockRate, 5);
        double partial1200Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 0.5 * clockRate,5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, alleleWC, 0.5 * clockRate, 5);


        double partial0000Internal2 = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 0.5 * clockRate,5) * partial0000Internal1+
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 0.5 * clockRate,5) * partial1000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele12, 0.5 * clockRate,5) * partial1200Internal1)
                * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1.0 * clockRate,5);

        double partial1000Internal2 = (substitutionModel.getSequenceTransitionProbability(allele1, allele0, 0.5 * clockRate,5) * partial0000Internal1+
                substitutionModel.getSequenceTransitionProbability(allele1, allele1, 0.5 * clockRate,5) * partial1000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele1, allele12, 0.5 * clockRate,5) * partial1200Internal1)
                * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1.0 * clockRate,5);

        double partial1200Internal2 = (substitutionModel.getSequenceTransitionProbability(allele12, allele0, 0.5 * clockRate,5) * partial0000Internal1+
                substitutionModel.getSequenceTransitionProbability(allele12, allele1, 0.5 * clockRate,5) * partial1000Internal1 +
                substitutionModel.getSequenceTransitionProbability(allele12, allele12, 0.5 * clockRate,5) * partial1200Internal1)
                * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele12, allele12, 1.0 * clockRate,5);


        double partial0000Root = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1.0 * clockRate,5) * partial0000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1.0 * clockRate,5) * partial1000Internal2 +
                substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1.0 * clockRate,5) * partial1200Internal2)
                * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele21, 2.0 * clockRate,5);

       double partial0000Origin = (substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2.0 * clockRate,5) * partial0000Root);





        //loglikelihood
        double LogPExpected = Math.log(partial0000Origin);


        double LogPCalc = likelihood.calculateLogP();
        Log.info.println("Log_lik:" + LogPCalc);

        assertEquals(LogPCalc, LogPExpected);


    }

//    @Test
//    public void testLikelihood3LeavesSameLengthInsertsMissing2_6tips() {
//
//
//        //Testing the ancestral state reconstruction at internal nodes
//        //tree with 3 tips
//        String newick = "(((CHILD24:0.5,CHILD17:0.5)INTERNAL1:0.5,(CHILD30:0.5,CHILD26:0.5)INT3:0.5)INTERNAL2:1,(CHILD16:1.0,CHILD3:1.0):1.0)";
//
//        Sequence a = new Sequence("CHILD3", "?,?,?,?,?");
//        Sequence b = new Sequence("CHILD16", "?,?,?,?,?");
//        Sequence c = new Sequence("CHILD26", "0,0,0,0,0");
//        Sequence d = new Sequence("CHILD30", "0,0,0,0,0");
//        Sequence e = new Sequence("CHILD17", "0,0,0,0,0");
//        Sequence f = new Sequence("CHILD24", "0,0,0,0,0");
//
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d,"sequence", e, "sequence", f, "dataType", "integer");
//
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();
//
//        //create a sub model with values
//        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
//        RealParameter editprobs = new RealParameter("0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.88");
//        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
//        //missingness parameters
//        RealParameter missRate = new RealParameter("0.0");
//        RealParameter missProbability = new RealParameter("0.5");
//
//        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,
//                "missingRate",missRate,"missingProbability",missProbability);        //site model
//        SiteModel siteM = new SiteModel();
//        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);
//
//        //likelihood class
//        RealParameter meanRate = new RealParameter("0.5");
//        StrictClockModel clockModel = new StrictClockModel();
//        clockModel.initByName("clock.rate", meanRate);
//        RealParameter origin = new RealParameter("10");
//        IntegerParameter arrayLength = new IntegerParameter("5");
//
//        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);
//
//
//        //initialise partialLikelihoods
//        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];
//
//        //Manually calc the likelihood for that tree:
//
//        //internal node partials:
//        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
//        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
//        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
//        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
//        double clockRate = 0.5;
//
//
//        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5);
//        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5);
//
//        double partial1200Internal1 = substitutionModel.getSequenceTransitionProbability(allele12, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele12, allele12, 1 * clockRate, 5);
//
//        double partial0000Internal2 = (partial0000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5) + partial1200Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5)) * (substitutionModel.getSequenceTransitionProbability(allele0, allele21, 2 * clockRate, 5));
//
//        //root node
//        double partialOrigin = partial0000Internal2 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2 * clockRate, 5);
//
//        //loglikelihood
//        double LogPExpected = Math.log(partialOrigin);
//
//        double LogPCalc = likelihood.calculateLogP();
//        Log.info.println("Log_lik:" + LogPCalc);
//
//        assertEquals(LogPExpected, LogPCalc);
//
//
//    }




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

    @Test
    public void testLikelihood3LeavesSameLengthInsertsAllDifferentMissing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "(((MISS:0.5,CHILD1:0.5)INT:0.5,CHILD3:1)INTERNAL:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,1,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Sequence d = new Sequence("MISS", "?,?,?,?,?");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d,"dataType", "integer");


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
        RealParameter missRate = new RealParameter("0.5");
        RealParameter missProb = new RealParameter("0.5");
        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies,"missingRate",missRate,"missingProbability",missProb);

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
        List<Integer> alleleWC = Arrays.asList(-1,-1, -1, -1, -1);


        double clockRate = 0.5;
        double missingRate = 0.5;
        double missingProb = 0.5;



        double partial0000Internal0 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele11, 0.5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 0.5 * clockRate, 5)  ;
        double partial1000Internal0 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele11, 0.5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, alleleWC, 0.5 * clockRate, 5);
        double partial1100Internal0 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele11, allele11, 0.5 * clockRate, 5) * substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, alleleWC, 0.5 * clockRate, 5);



        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele12, 1 * clockRate, 5) * (
                partial0000Internal0 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 0.5 * clockRate, 5) +
                partial1000Internal0 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 0.5 * clockRate, 5) +
                partial1100Internal0* substitutionModel.getSequenceTransitionProbability(allele0, allele11, 0.5 * clockRate, 5)

        );
        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbabilityTipEdge(allele1, allele12, 1 * clockRate, 5) * (
                partial1000Internal0 * substitutionModel.getSequenceTransitionProbability(allele1, allele1, 0.5 * clockRate, 5) +
                partial1100Internal0* substitutionModel.getSequenceTransitionProbability(allele1, allele11, 0.5 * clockRate, 5)) ;

        double partial0000Internal2 = (partial0000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) +
                partial1000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)) *
                (substitutionModel.getSequenceTransitionProbabilityTipEdge(allele0, allele21, 2 * clockRate, 5));


        //root node
        double partialOrigin = partial0000Internal2 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2 * clockRate, 5);

        //loglikelihood
        double LogPExpected = Math.log(partialOrigin);

        double LogPCalc = likelihood.calculateLogP();
        Log.info.println("Log_lik:" + LogPCalc);

        assertEquals(LogPExpected, LogPCalc);

    }

//    @Test
//    public void testLikelihood3LeavesSameLengthInsertsAllDifferentERRORPROPAGATION() {
//
//
//        //Testing the ancestral state reconstruction at internal nodes
//        //tree with 3 tips
//        String newick = "(((MISS:0.5,CHILD1:0.5)INT:0.5,CHILD3:1)INTERNAL:1,CHILD2:2.0)";
//
//        Sequence a = new Sequence("CHILD1", "1,1,0,0,0");
//        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
//        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
//        Sequence d = new Sequence("MISS", "?,?,?,?,?");
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d,"dataType", "integer");
//
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        SciPhyTreeLikelihood likelihood = new SciPhyTreeLikelihood();
//
//        //create a sub model with values
//        SciPhySubstitutionModel substitutionModel = new SciPhySubstitutionModel();
//        RealParameter editprobs = new RealParameter("0.8 0.2");
//        RealParameter stateFrequencies = new RealParameter("1.0 0 0");
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
//        substitutionModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);
//
//        //site model
//        SiteModel siteM = new SiteModel();
//        siteM.initByName("gammaCategoryCount", 0, "substModel", substitutionModel);
//
//        //likelihood class
//
//        RealParameter meanRate = new RealParameter("0.5");
//        StrictClockModel clockModel = new StrictClockModel();
//        clockModel.initByName("clock.rate", meanRate);
//        RealParameter origin = new RealParameter("4");
//        IntegerParameter arrayLength = new IntegerParameter("5");
//
//        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "arrayLength", arrayLength);
//
//
//        //initialise partialLikelihoods
//        likelihood.partialLikelihoods = new double[2][tree1.getNodeCount()][];
//
//        //Manually calc the likelihood for that tree:
//
//        //internal node partials:
//        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);
//        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
//        List<Integer> allele11 = Arrays.asList(1, 1, 0, 0, 0);
//        List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
//        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
//        double clockRate = 0.5;
//
//        double partial0000Internal1 = substitutionModel.getSequenceTransitionProbability(allele0, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele0, allele11, 1 * clockRate, 5);
//        double partial1000Internal1 = substitutionModel.getSequenceTransitionProbability(allele1, allele12, 1 * clockRate, 5) * substitutionModel.getSequenceTransitionProbability(allele1, allele11, 1 * clockRate, 5);
//
//        double partial0000Internal2 = (partial0000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 1 * clockRate, 5) + partial1000Internal1 * substitutionModel.getSequenceTransitionProbability(allele0, allele1, 1 * clockRate, 5)) * (substitutionModel.getSequenceTransitionProbability(allele0, allele21, 2 * clockRate, 5));
//
//
//        //root node
//        double partialOrigin = partial0000Internal2 * substitutionModel.getSequenceTransitionProbability(allele0, allele0, 2 * clockRate, 5);
//
//        //loglikelihood
//        double LogPExpected = Math.log(partialOrigin);
//
//        double LogPCalc = likelihood.calculateLogP();
//        assertEquals(LogPExpected, LogPCalc);
//
//    }

}