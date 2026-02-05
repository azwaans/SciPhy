package sciphy;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.AlignmentImporter;
import org.junit.Test;
import sciphy.evolution.simulation.SimulatedSciPhyAlignment;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.io.File;

import static junit.framework.Assert.assertEquals;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry() {

        Integer sequenceLength = 1;
        String outputFileName = "test/simAl.sciphy";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        RealParameter insertrates = new RealParameter("0.8 0.2");
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(2);
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("test/sciphy/expectedAlignment.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }


}
