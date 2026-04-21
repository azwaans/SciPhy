package sciphy.evolution.likelihood;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.util.*;

import static sciphy.util.LogSum.logSum;

@Description("tree likelihood for a SciPhy alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a version of the SciPhy likelihood using caching without a likelihoodCore implementation ")

public class SciPhyParsimony extends GenericTreeLikelihood {


    final public Input<IntegerParameter> arrayLengthInput = new Input<>("arrayLength", "Number of positions in the target BC", Validate.REQUIRED);
    final public Input<RealParameter> originTimeInput = new Input<>("origin", "Duration of the experiment");


    protected int nodeCount;
    protected int arrayLength;
    protected double parsimonyScore;
    protected double originTime;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;



    //to be able to have current/stored states in an analog way to the partials array, ancestral states are accessed/added
    //states with key : (NodeNr + 1) + (current ? 0:1) * (NodeNr+1)
    public Hashtable<Integer, List<List<Integer>>> ancestralStates;
    public double[][][] partialParsimonies;



    protected int[] currentParsimoniesIndex;
    protected int[] storedParsimoniesIndex;

    protected int[] currentStatesIndex;
    protected int[] storedStatesIndex;


    @Override
    public void initAndValidate() {

        arrayLength = arrayLengthInput.get().getValue();
        if (arrayLength < 1 || (dataInput.get().getSiteCount() != arrayLength)) {
            throw new IllegalArgumentException(String.format(
                    "Invalid array length: Ensure that length >= 1 and matches alignment "));
        }
        nodeCount = treeInput.get().getNodeCount();
        if (nodeCount <= 2) {
            throw new IllegalArgumentException(String.format(
                    "Invalid tree input: single node/branch. Ensure that #nodes>2 "));
        }

        originTime = 0.0;
        if (originTimeInput.get() != null) {
            originTime = originTimeInput.get().getValue();
            if (originTime < 0.0) {
                throw new IllegalArgumentException(String.format(
                        "Invalid origin time input: ensure that origin>0"));
            }
        }


        //TODO check that state count from alignment (i.e. data type) and substitution model are the same
        ancestralStates = new Hashtable<>();
        partialParsimonies = new double[2][nodeCount][];

        currentParsimoniesIndex = new int[nodeCount];
        storedParsimoniesIndex = new int[nodeCount];

        currentStatesIndex = new int[nodeCount];
        storedStatesIndex = new int[nodeCount];
        parsimonyScore = 0;




        hasDirt = Tree.IS_FILTHY;

        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
            initLeafAncestors(i);
        }

        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
            initLeafPartials(i);
        }

    }


    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    public double calculateLogP() {


        logP = 0.0;
        return logP;
    }

    public double calculateParsimony() {
        final TreeInterface tree = treeInput.get();



        //adjust clock rate for the given category
        traverse(tree.getRoot(), 0);

        if (originTime == 0.0) {
            //sum of all partial likelihoods at the root
            int rootNr = tree.getRoot().getNr();
            parsimonyScore = Arrays.stream(partialParsimonies[currentParsimoniesIndex[rootNr]][rootNr]).sum();
        } else {
            //the tree log likelihood is the log(p) of unedited state at the origin
            parsimonyScore = calculateOriginPartial(tree.getRoot(), 0);

        }


        return parsimonyScore;

    }




    /**
     * Calculate partial likelihoods for a given leaf node, and fill the corresponding partialLikelihood array
     */
    protected void initLeafPartials(int nodeNr) {

        double[] leafpartialParsimonies = initpartialParsimoniesLeaf(ancestralStates.get(makeCachingIndexStates(nodeNr)).size());
        this.partialParsimonies[0][nodeNr] = new double[leafpartialParsimonies.length];
        this.partialParsimonies[1][nodeNr] = new double[leafpartialParsimonies.length];
        System.arraycopy(leafpartialParsimonies, 0, this.partialParsimonies[0][nodeNr], 0, leafpartialParsimonies.length);

    }

    /**
     * Calculate the set of ancestral states for a given leaf node, and fill the corresponding AncestralStates hashmap
     */
    protected void initLeafAncestors(int nodeNr) {

        List<List<Integer>> possibleLeafAncestors = getPossibleAncestors(dataInput.get().getCounts().get(nodeNr));
        ancestralStates.put(makeCachingIndexStates(nodeNr), possibleLeafAncestors);

    }


    /**
     * This implements a postorder traversal of the tree to fill the ancestralStates hashmap and corresponding partialLikelihood array.
     */
    protected int traverse(Node node, int categoryId) {

        int update = (node.isDirty() | hasDirt);
        int nodeIndex = node.getNr();


        if (!node.isRoot() && (update != Tree.IS_CLEAN )) {

            update |= Tree.IS_DIRTY;
        }

        if (!node.isLeaf()) {

            final Node child1 = node.getLeft();
            final int update1 = traverse(child1, categoryId);
            final Node child2 = node.getRight();
            final int update2 = traverse(child2, categoryId);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                update |= (update1 | update2);

                if (update >= Tree.IS_FILTHY) {
                    setNodeStatesForUpdate(nodeIndex);
                    calculateStates(nodeIndex, child1.getNr(), child2.getNr());
                }

                setNodePartialsForUpdate(nodeIndex);
                calculatePartialsParsimony(nodeIndex, child1, child2, categoryId);



            }
        }

        return update;
    }

    /**
     * Construct a set of possible ancestral states at an internal node by intersection of children sets, updates the
     * AncestralStates hashmap with the resulting set.
     */
    public void calculateStates(int nodeNr, int child1Nr, int child2Nr) {

        List<List<Integer>> ancSetChild1 = ancestralStates.get(makeCachingIndexStates(child1Nr));
        List<List<Integer>> ancSetChild2 = ancestralStates.get(makeCachingIndexStates(child2Nr));

        List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);

        // intersection of children ancestral states
        ancSetNode.retainAll(ancSetChild2);

        if(ancSetNode.size() ==0 ) {
            // in this case the state is the unedited
            List<Integer> startState = Arrays.asList(0, 0, 0, 0, 0);
            List<List<Integer>> parsimonySetNode = new ArrayList<>();
            parsimonySetNode.add(startState);
            ancestralStates.put(makeCachingIndexStates(nodeNr), parsimonySetNode);
        }
        else {
            int indexLongest = 0;

            int lengthLongest = 0;

            for (int index = 0; index < ancSetNode.size(); ++index) {
                //find the most edited state

                //to do so, we need to have access to the number of edited positions
                //create an unedited state to subtract from sequences to get only edited sites
                List<Integer> zero = Arrays.asList(0);

                //removing all unedited sites from each sequence
                List<Integer> state = new ArrayList<>(ancSetNode.get(index));
                state.removeAll(zero);

                if (state.size() > lengthLongest) {
                    indexLongest = index;
                    lengthLongest = state.size();
                }

            }

            List<List<Integer>> parsimonySetNode = new ArrayList<>();
            parsimonySetNode.add(ancSetNode.get(indexLongest));


            ancestralStates.put(makeCachingIndexStates(nodeNr), parsimonySetNode);
        }
    }

    public void setNodePartialsForUpdate(int nodeIndex) {
        currentParsimoniesIndex[nodeIndex] = 1 - currentParsimoniesIndex[nodeIndex];
    }

    public void setNodeStatesForUpdate(int nodeIndex) {
        currentStatesIndex[nodeIndex] = 1 - currentStatesIndex[nodeIndex];
    }

    public int makeCachingIndexStates(int nodeIndex) {
        int node = nodeIndex + 1;
        String forHashing = node + "" +  currentStatesIndex[nodeIndex] + ""+ node;
        return forHashing.hashCode();

    }


    /**
     * This function calculates partial likelihoods for all possible states at a node given its children partials
     * and sets the corresponding partial likelihoods, for all possible states at node nodeNr
     */
    public void calculatePartialsParsimony(int nodeNr, Node child1, Node child2, int categoryId) {

        //initialize an array for the partials
        double[] partials = new double[ancestralStates.get(makeCachingIndexStates(nodeNr)).size()];

        for (int stateIndex = 0; stateIndex < ancestralStates.get(makeCachingIndexStates(nodeNr)).size(); ++stateIndex) {

            List<Integer> startState = ancestralStates.get(makeCachingIndexStates(nodeNr)).get(stateIndex);

            double child1partialParsimoniestate = calculatepartialParsimoniestate(startState, child1, categoryId);
            double child2partialParsimoniestate = calculatepartialParsimoniestate(startState, child2, categoryId);

            partials[stateIndex] = child1partialParsimoniestate + child2partialParsimoniestate;
        }

        partialParsimonies[currentParsimoniesIndex[nodeNr]][nodeNr] = partials;

    }

    /**
     * This function calculates the likelihood of the unedited state at the origin given partial likelihoods at the root
     * node
     *
     * @return likelihood of the unedited barcode at t = origin
     */

    public double calculateOriginPartial(Node rootNode, int categoryId) {

        //the start state is the unedited sciphy barcode
        List<Integer> startState = Arrays.asList(0, 0, 0, 0, 0);
        double partialAtOrigin = calculatepartialParsimoniestate(startState, rootNode, categoryId);

        return partialAtOrigin;

    }

    /**
     * This function calculates the partial likelihood term of a specific state at a node derived on a branch leading to
     * a child node
     *
     * @return partial likelihood for a state at a node given partials at a node childNode
     */
    public double calculatepartialParsimoniestate(List<Integer> startState, Node childNode, int categoryId) {

        double statePartialLikelihood = 0;
        // calculate partials
        if (childNode.isLeaf()) {

            List<Integer> endState = ancestralStates.get(makeCachingIndexStates(childNode.getNr())).get(0);
            statePartialLikelihood += getSequenceParsimonyTransition(startState, endState, this.arrayLength);
        } else {

            for (int endStateIndex = 0; endStateIndex < ancestralStates.get(makeCachingIndexStates(childNode.getNr())).size(); ++endStateIndex) {

                List<Integer> endState = ancestralStates.get(makeCachingIndexStates(childNode.getNr())).get(endStateIndex);

                // if the end state has non-null partial likelihood
                //if (partialParsimonies[currentParsimoniesIndex[childNode.getNr()]][childNode.getNr()][endStateIndex] != 0.0) {

                statePartialLikelihood = partialParsimonies[currentParsimoniesIndex[childNode.getNr()]][childNode.getNr()][endStateIndex] + getSequenceParsimonyTransition(startState, endState, this.arrayLength) ;

                // }
            }
        }
        return statePartialLikelihood;
    }

    public double[] getSafePartials(int nodeNr) {
        int currentIndex = currentParsimoniesIndex[nodeNr];
        double[] buffer = partialParsimonies[currentIndex][nodeNr];

        // If current is null or empty, use the stored one
        if (buffer == null || buffer.length == 0) {
            int storedIndex = 1 - currentIndex;
            buffer = partialParsimonies[storedIndex][nodeNr];
        }

        // Final fallback: if both are empty, return a tiny dummy array or null
        return (buffer != null && buffer.length > 0) ? buffer : null;
    }

    public double getSequenceParsimonyTransition(final List<Integer> startSequence, final List<Integer> endSequence, int arrayLength) {

        List<Integer> startState = new ArrayList(startSequence);
        List<Integer> endState = new ArrayList(endSequence);

        //create an unedited state to subtract from sequences to get only edited sites
        List<Integer> zero = Arrays.asList(0);

        //removing all unedited sites from each sequence
        startState.removeAll(zero);
        endState.removeAll(zero);

        //if endState is less edited than the start state, violates ordering
        if(startState.size() > endState.size() ){
            return 0.0;
        }

        //subtracting start sequence from end sequence: edits introduced
        // if start state has identical elements to end state remove
        startState.forEach(endState::remove);
        List<Integer> newInserts = endState;

        //available positions are targetBClength length - number of edited positions
        return newInserts.size();


    }

    /**
     * This function initialises an array of partial likelihoods for a leaf node, the partial likelihood is 1 for
     * the observed sequence and 0 for everything else. The size corresponds to the total number of possible ancestral states.
     *
     * @return array of partial likelihoods at leaf node
     */
    public double[] initpartialParsimoniesLeaf(int size) {

        double[] leafPartials = new double[size];
        leafPartials[0] = 0;
        return leafPartials;
    }

    /**
     * This function returns all possible ancestral states given a sequence.
     * Because sciphy sequences record ordered edits, ancestral states are obtained by sequentially removing edits
     * along the sequence (from any insert (1 to N)  to 0)
     *
     * @return a list of possible ancestral sciphy barcode states
     */
    public static List<List<Integer>> getPossibleAncestors(List<Integer> sequence) {

        List<List<Integer>> ancestors = new ArrayList();
        ancestors.add(sequence);

        List<Integer> ancestor = new ArrayList<>(sequence);
        for (int i = sequence.size() - 1; i >= 0; --i) {
            if (sequence.get(i) != 0) {
                ancestor.set(i, 0);
                ancestors.add(new ArrayList<>(ancestor));
            }
        }
        return ancestors;
    }



    /**
     * check state for changed variables and update temp results if necessary *
     */
    //requires recalculation if the data has changed, if the sitemodel (or sub model has changed), the branch rate model, or the tree has changed.
    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }


        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {

        super.store();
        System.arraycopy(currentParsimoniesIndex, 0, storedParsimoniesIndex, 0, nodeCount);
        System.arraycopy(currentStatesIndex, 0, storedStatesIndex, 0, nodeCount);
    }

    //TODO do we need unstore??? We think we don't because when scaling is active, it is for the entire likelihood

    @Override
    public void restore() {

        super.restore();

        int[] tmp2 = currentParsimoniesIndex;
        currentParsimoniesIndex = storedParsimoniesIndex;
        storedParsimoniesIndex = tmp2;

        int[] tmp3 = currentStatesIndex;
        currentStatesIndex = storedStatesIndex;
        storedStatesIndex = tmp3;
    }


}