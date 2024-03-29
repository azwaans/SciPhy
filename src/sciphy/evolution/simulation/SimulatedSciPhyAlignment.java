package sciphy.evolution.simulation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.pkgmgmt.BEASTClassLoader;
import beast.pkgmgmt.PackageManager;
import beast.base.util.Randomizer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@Description("A more flexible alignment simulator adapted from Tim Vaughan's feast implementation")
public class SimulatedSciPhyAlignment extends Alignment {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree down which to simulate sequence evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "Site model to use in simulation.",
            Input.Validate.REQUIRED);

    //TODO rename number of targets
    public Input<Integer> nrOfTargetsInput = new Input<>(
            "numberOfTargets",
            "Number of targets to simulate.",
            Input.Validate.REQUIRED);

    public Input<Integer> arrayLengthInput = new Input<>(
            "arrayLength",
            "Number of insertions to add per site",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Start of the process, usually the experiment",
            Input.Validate.OPTIONAL);


    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "Name of file (if any) simulated alignment should be saved to.");

    private Tree tree;
    private SiteModel siteModel;
    private int numberOfTargets;
    private int arrayLength;
    private DataType dataType;
    private double originHeight;

    private String ancestralSeqStr;

    public SimulatedSciPhyAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        siteModel = siteModelInput.get();
        Log.info.println("category rates: " + Arrays.toString(siteModel.getCategoryRates(tree.getRoot())));

        numberOfTargets = 1; //TODO rewrite for arbitrary #targets
        arrayLength = arrayLengthInput.get();
        sequences.clear();

        if (originInput.get() != null) {
            originHeight = originInput.get().getValue();
        }

        grabDataType();

        simulate();

        super.initAndValidate();

        // Write simulated alignment to disk if required
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(new TaxonSet(this)));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + outputFileNameInput.get() + ".");
            }
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        int nTaxa = tree.getLeafNodeCount();

        SciPhySubstitutionModel substModel = (SciPhySubstitutionModel) siteModel.getSubstitutionModel();

        double[] transitionProbs = substModel.getInsertProbabilities();

        int[][] alignment = new int[nTaxa][arrayLength];

        Node root = tree.getRoot();

        int[] rootSequence = new int[arrayLength];

        if (originHeight != 0) {
            // then parent sequence is sequence at origin and we evolve sequence first down to the root
            double deltaT = originHeight - root.getHeight();
            double clockRate = siteModel.getRateForCategory(0, root);
            Log.info.println("clock rates: " + clockRate);

            int nPossibleInserts = arrayLength;
            long nPotentialInserts = Randomizer.nextPoisson(deltaT * clockRate);

            int insertionIndex = 0;

            // Add potential inserts while there are still possible insertion positions
            while (nPossibleInserts > 0 && nPotentialInserts > 0) {

                int newInsertion = Randomizer.randomChoicePDF(transitionProbs) + 1;
                rootSequence[insertionIndex] = newInsertion;

                insertionIndex++;
                nPossibleInserts--;
                nPotentialInserts--;
            }
        }

        //ancestralSeqStr = dataType.encodingToString(parentSequence);

        traverse(root, rootSequence,
                transitionProbs,
                alignment);

        for (int leafIdx = 0; leafIdx < nTaxa; leafIdx++) {
            String seqString = dataType.encodingToString(alignment[leafIdx]);

            String taxonName;
            if (tree.getNode(leafIdx).getID() != null)
                taxonName = tree.getNode(leafIdx).getID();
            else
                taxonName = "t" + leafIdx;

            sequenceInput.setValue(new Sequence(taxonName, seqString), this);
        }
    }

    /**
     * Traverse a tree, simulating a sequence alignment down it.
     *
     * @param node            Node of the tree
     * @param parentSequence  Sequence at the parent node in the tree
     * @param transitionProbs transition probabilities
     * @param regionAlignment alignment for particular region
     */
    private void traverse(Node node,
                          int[] parentSequence,
                          double[] transitionProbs,
                          int[][] regionAlignment) {


        // ignore categories so far

        for (Node child : node.getChildren()) {

            double deltaT = node.getHeight() - child.getHeight();
            double clockRate = siteModel.getRateForCategory(0, child);

            // Draw characters on child sequence
            int[] childSequence = parentSequence.clone();

            // find site where next insertion could happen, i.e. the next unedited state '0'
            int insertionIndex = 0;
            while ((insertionIndex < arrayLength) && (childSequence[insertionIndex] != 0)) {
                insertionIndex++;
            }

            if (insertionIndex == arrayLength) {
                // then sequence is fully edited, no further simulation necessary
                ;
            } else {
                // sample number of new insertions
                int nPossibleInserts = arrayLength - insertionIndex;
                long nPotentialInserts = Randomizer.nextPoisson(deltaT * clockRate);

                // Add potential inserts while there are still possible insertion positions
                while (nPossibleInserts > 0 && nPotentialInserts > 0) {

                    int newInsertion = Randomizer.randomChoicePDF(transitionProbs) + 1;
                    childSequence[insertionIndex] = newInsertion;

                    insertionIndex++;
                    nPossibleInserts--;
                    nPotentialInserts--;
                }
            }

            if (child.isLeaf()) {
                System.arraycopy(childSequence, 0,
                        regionAlignment[child.getNr()], 0, childSequence.length);
            } else {
                traverse(child, childSequence,
                        transitionProbs,
                        regionAlignment);
            }
        }
    }

    /**
     * HORRIBLE function to identify data type from given description.
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = PackageManager.find(beast.base.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                        | InstantiationException
                        | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }
}
