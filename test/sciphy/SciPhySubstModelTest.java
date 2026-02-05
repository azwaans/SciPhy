package sciphy;

import beast.base.core.Description;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.substitutionmodel.Frequencies;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;


@Description("Test substitution model")
public class SciPhySubstModelTest {

    // test that inputs are handled as expected

    SciPhySubstitutionModel substModel;

    @Before
    public void setUp(){
        substModel = new SciPhySubstitutionModel();
    }

    //----------------------------------------------------------------------------------//
    // test valid inputs
    @Test (expected = RuntimeException.class)
    public void testExceptionForNegEditFreqInput(){

        RealParameter editprobs = new RealParameter("-1 2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", editprobs, "estimate", false);

        // this should fail
        substModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);
    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForEditFreqSumNot1(){

        RealParameter editprobs = new RealParameter("0.6 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", editprobs, "estimate", false);

        // this should fail
        substModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);
    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForEmptyEditFreqInput(){

        RealParameter editprobs = new RealParameter();
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", editprobs, "estimate", false);

        // this should fail
        substModel.initByName("editProbabilities", editprobs, "frequencies", frequencies);
    }

    //----------------------------------------------------------------------------------//
    // test transition probability calculations
    @Test
    public void testTransitionProbabilities(){


        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();
        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);
        int targetBClength = 5;


        Sequence a = new Sequence("cell1", "1,2,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,2,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3
        // P(1)*0.2

        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);
        Double expectedProbability = dist.probability(1) * 0.2 ;
        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5,targetBClength);


        // Assert
        assertEquals(expectedProbability, calculatedProbability,0.00001);

    }


    @Test
    public void testTransitionProbabilitiesNoEdit(){

        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);
        int targetBClength = 5;



        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "2,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5, targetBClength);
        Double expectedProbability = 0.6065306597126334;

        assertEquals(expectedProbability, calculatedProbability, 1e-10);

    }

    @Test
    public void testTransitionProbabilitiesSingleEdit(){

        // Arrange
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();

        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);
        int targetBClength = 5;



        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "2,1,1,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5, targetBClength);
        // expectedProbability : draw 1 event on a Poisson process, with event frequency 0.8
        // P(1) * 0.8
        Double expectedProbability = 0.2426122638850534;


        assertEquals(expectedProbability, calculatedProbability, 1e-10);


    }


    @Test
    public void testTransitionProbabilities2Edits(){

        // Arrange
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();

        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);

        int targetBClength = 5;

        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "2,1,1,2,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5, targetBClength);

        // expectedProbability : draw 2 events on a Poisson process, with event probabilities 0.8 and 0.2 resp
        // P(2) * 0.8 * 0.2
        Double expectedProbability = 0.012130613194252673;

        assertEquals(expectedProbability, calculatedProbability, 1e-10);


    }

    @Test
    public void testTransitionProbabilitiesEditsAndSaturation(){

        // Arrange
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        Frequencies frequencies = new Frequencies();
        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);

        int targetBClength = 5;

        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "2,1,1,2,2");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5, targetBClength);

        // expectedProbability : draw 3 events on a Poisson process bounded to 3:
        // with event frequency 0.8, 0.2
        // (1 - P(0) + P(1) + P(2)) * (0.2*0.2*0.8)
        Double expectedProbability = 4.6040569494306294E-4;

        assertEquals(expectedProbability, calculatedProbability, 1e-10);


    }

    @Test
    public void test1ProbabilityForFullySaturated(){

        // Arrange
        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
        RealParameter editProbabilities = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);

        int targetBClength = 5;

        Sequence a = new Sequence("cell1", "2,1,1,2,2");
        Sequence b = new Sequence("cell2", "2,1,1,2,2");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        Double calculatedProbability = substModel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5, targetBClength);

        Double expectedProbability = 1.0;

        assertEquals(expectedProbability, calculatedProbability, 1e-10);

    }

    // TODO test this elsewhere: the test below fails, but is irrelevant assuming that only correct sequences of states will ever be proposed, as tested in LikelihoodTests
//    @Test
//    public void test0ProbabilityForForbiddenTransition(){
//
//        RealParameter stateFrequencies = new RealParameter("1.0 0 0 ");
//        Frequencies frequencies = new Frequencies();
//        RealParameter editProbabilities = new RealParameter("0.8 0.2");
//
//        frequencies.initByName("frequencies", stateFrequencies, "estimate", false);
//        substModel.initByName( "editProbabilities", editProbabilities, "frequencies" ,frequencies);
//
//        int targetBClength = 5;
//
//        Sequence a = new Sequence("cell1", "2,0,0,0,0");
//        Sequence b = new Sequence("cell2", "1,2,0,0,0");
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "dataType", "integer");
//        alignment.initByName("sequence", b, "dataType", "integer");
//
//        //internal representation of the sequences for the package:
//        List<Integer> sequence_a = alignment.getCounts().get(0);
//        List<Integer> sequence_b = alignment.getCounts().get(1);
//
//        Double calculatedProbability = substModel.getSequenceTransitionProbability(sequence_a, sequence_b,0.5, targetBClength);
//
//        Double expectedProbability = 0.0;
//
//        assertEquals(expectedProbability, calculatedProbability, 1-10);
//    }

    //-----------------------------------------------------------------------------------//
    // Validate that the PoissonDistributionImpl we use to calculate the poisson probabilities
    // behaves as expected. We test this by comparing against the values by Rpois.

    @Test
    public void testPoissonDistAgainstR(){

        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);
        double p0 = dist.probability(0);
        assertEquals( 0.6065307, p0, 0.00001);
        double p1 = dist.probability(1);
        assertEquals( 0.3032653, p1, 0.00001);
        double p2 = dist.probability(2);
        assertEquals(  0.07581633, p2, 0.00001);
        double p3 = dist.probability(3);
        assertEquals(  0.01263606, p3, 0.00001);
        double p4 = dist.probability(4);
        assertEquals(  0.001579507, p4, 0.00001);


         dist = new PoissonDistributionImpl(0.001);
         p0 = dist.probability(0);
        assertEquals( 0.9990005, p0, 0.00001);
         p1 = dist.probability(1);
        assertEquals( 0.0009990005, p1, 0.00001);
         p2 = dist.probability(2);
        assertEquals(  4.995002e-07, p2, 0.00001);
         p3 = dist.probability(3);
        assertEquals(  1.665001e-10, p3, 0.00001);
         p4 = dist.probability(4);
        assertEquals(  4.162502e-14, p4, 0.00001);
    }
}