package java.sciphy;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.Real;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.AlignmentImporter;
import org.junit.Test;
import sciphy.evolution.simulation.SimulatedSciPhyAlignment;
import sciphy.evolution.simulation.SimulatedSciPhyAlignmentHeritableMissingBcodes;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.io.File;

import static junit.framework.Assert.assertEquals;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry() {

        //Integer sequenceLength = 1;
        String outputFileName = "src/test/sciphy/simAl.sciphy";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        SimplexParam insertrates = new SimplexParam(new double[]{0.8, 0.2});
        SimplexParam freqs = new SimplexParam(new double[] {1.0, 0, 0});
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealScalarParam mutationRate = new RealScalarParam(10.0, Real.INSTANCE);
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(2);
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
               // "sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/sciphy/expectedAlignment.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryHeritableMissingOnly() {

        //Integer sequenceLength = 1;
        String outputFileName = "test/simAlWithHeritableMissingOnly.sciphy";
        String newick = "(CHILD1:5,CHILD2:5):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        SimplexParam insertrates = new SimplexParam(new double[]{0.8, 0.2});
        SimplexParam freqs = new SimplexParam(new double[] {1.0, 0, 0});
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealScalarParam mutationRate = new RealScalarParam(10.0,Real.INSTANCE);


        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);


        RealScalarParam missingRate = new RealScalarParam(0.03, Real.INSTANCE);
        RealScalarParam origin = new RealScalarParam(6.0,Real.INSTANCE);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(2);
        SimulatedSciPhyAlignmentHeritableMissingBcodes simAlignment = new SimulatedSciPhyAlignmentHeritableMissingBcodes();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
               // "sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingRate",missingRate,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/sciphy/expectedAlignmentMissingRates.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryTipMissingOnly() {

        //Integer sequenceLength = 1;
        String outputFileName = "src/test/sciphy/simAlWithTipMissingOnly.sciphy";
        String newick = "(CHILD1:5,CHILD2:5):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        SimplexParam insertrates = new SimplexParam(new double[]{0.8, 0.2});
        SimplexParam freqs = new SimplexParam(new double[] {1.0, 0, 0});
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealScalarParam mutationRate = new RealScalarParam(10.0,Real.INSTANCE);


        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);


        RealScalarParam missingProb = new RealScalarParam(0.5,UnitInterval.INSTANCE);
        RealScalarParam origin = new RealScalarParam(6.0,Real.INSTANCE);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(12345);
        SimulatedSciPhyAlignmentHeritableMissingBcodes simAlignment = new SimulatedSciPhyAlignmentHeritableMissingBcodes();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
               // "sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingProbability",missingProb,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/sciphy/expectedAlignmentMissingRates.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryBothMissingTypes() {

        //Integer sequenceLength = 1;
        String outputFileName = "src/test/sciphy/simAlBothMissingTypes.sciphy";
        String newick = "(CHILD1:5,CHILD2:5):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        SimplexParam insertrates = new SimplexParam(new double[]{0.8, 0.2});
        SimplexParam freqs = new SimplexParam(new double[] {1.0, 0, 0});
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealScalarParam mutationRate = new RealScalarParam(10.0,Real.INSTANCE);


        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);


        RealScalarParam missingRate = new RealScalarParam(0.03, Real.INSTANCE);
        RealScalarParam missingProb = new RealScalarParam(0.5, UnitInterval.INSTANCE);
        RealScalarParam origin = new RealScalarParam(6.0,Real.INSTANCE);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(2);
        SimulatedSciPhyAlignmentHeritableMissingBcodes simAlignment = new SimulatedSciPhyAlignmentHeritableMissingBcodes();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                //"sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingProbability",missingProb,
                "missingRate",missingRate,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/sciphy/expectedAlignmentMissingRates.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }


}
