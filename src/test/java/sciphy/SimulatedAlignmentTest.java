package sciphy;

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
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;
import java.io.File;
import static junit.framework.Assert.assertEquals;


public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry() {
        Randomizer.setSeed(2);
        String outputFileName = "src/test/java/sciphy/simAl.sciphy";
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
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "arrayLength", 100,
                "outputFileName", outputFileName,
                "userDataType", integerData
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/java/sciphy/expectedAlignment.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryHeritableMissingDataOnly() {
        Randomizer.setSeed(2);

        String outputFileName = "src/test/java/sciphy/simAlHeritableMissingData.sciphy";
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
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "arrayLength", 100,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingRate",missingRate,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/java/sciphy/expectedAlHeritableMissingData.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryDropoutMissingDataOnly() {
        Randomizer.setSeed(12345);

        String outputFileName = "src/test/java/sciphy/simAlDropoutMissingData.sciphy";
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
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "arrayLength", 100,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingProbability",missingProb,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/java/sciphy/expectedAlDropoutMissingData.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

    @Test
    public void testSimulationOnCherryBothMissingDataTypes() {
        Randomizer.setSeed(2);

        String outputFileName = "src/test/java/sciphy/simAlBothMissingTypes.sciphy";
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
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "arrayLength", 100,
                "outputFileName", outputFileName,
                "userDataType", integerData,
                "missingProbability",missingProb,
                "missingRate",missingRate,
                "origin",origin
        );

        AlignmentImporter expectedAlignment = new sciphy.util.NexusImporter();
        beast.base.evolution.alignment.Alignment alignment = (Alignment) expectedAlignment.loadFile(new File("src/test/java/sciphy/expectedAlBothMissingTypes.sciphy")).get(0);

        assertEquals(alignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(alignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }


}
