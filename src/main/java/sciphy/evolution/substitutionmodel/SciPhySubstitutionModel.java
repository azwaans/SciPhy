package sciphy.evolution.substitutionmodel;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.evolution.substitutionmodel.Base;
import beast.base.evolution.tree.Node;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.Simplex;
import org.apache.commons.math3.distribution.PoissonDistribution;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



@Description("Calculates transition probabilities sequence transitions from an ancestral node to a child node. " +
        "Assume that the sequences proposed represent valid transition pairs.")

public class SciPhySubstitutionModel extends Base {
    final public Input<Simplex> editProbabilitiesInput = new Input<>("editProbabilities",
            "Edit probabilities for the typewriter process", Input.Validate.REQUIRED);

    final public Input<RealScalar<NonNegativeReal>>  missingRateInput = new Input<>(
            "missingRate",
            "Rate at which the barcode goes missing heritably. Default 0.0.",
            new RealScalarParam(0.0, NonNegativeReal.INSTANCE));

    final public Input<RealScalar<NonNegativeReal>>  missingProbInput = new Input<>(
            "missingProbability",
            "Probability that a barcode goes missing at the tips. Default 0.0",
            new RealScalarParam(0.0, NonNegativeReal.INSTANCE));

    /**
     * edit insertion rate  *
     */
    protected Simplex editProbabilities;
    protected RealScalar<NonNegativeReal> missRate;
    protected RealScalar<NonNegativeReal> missProbability;
    double[] editProbs;
    double missingRate;
    double missingProbability;


    public List<Integer> missingState;
    public List<Integer> lostState;

    @Override
    public void initAndValidate() {

        super.initAndValidate();
        // Get edit probabilities and check correct input
        editProbabilities = editProbabilitiesInput.get();
        editProbs = new double[editProbabilities.size()];
        for (int i = 0; i < editProbabilities.size(); i++) {
            editProbs[i] = editProbabilities.get(i);
        }


        double insertProbabilitiesSum = Arrays.stream(editProbs).sum();
        if (Math.abs(insertProbabilitiesSum - 1.0) > 1e-6) {
            throw new IllegalArgumentException(String.format(
                    "sum of insert probabilities is not 1"));

        }

        for( double insertProbability : editProbs) {
            if(insertProbability>1 || insertProbability<0) {
                throw new IllegalArgumentException(String.format(
                        "insert probabilities is not between 0 and 1"));
            }
        }


        missingRate = 0.0;
        if (missingRateInput.get() != null) {
            missRate = missingRateInput.get();
            missingRate = missRate.get();
        }

        missingProbability = 0.0;
        if (missingProbInput.get() != null) {
            missProbability = missingProbInput.get();
            missingProbability = missProbability.get();
        }

    }

    @Override
    public double[] getRateMatrix(Node node) {
        return new double[0];
    }


    /**
     * This function calculates the probability of transitioning between 2 sequences states in given evolutionary time (distance)
     * (with potentially multiple edits having happened), on internal branches of a tree
     *
     * @param startSequence  is a sequence state at a parent node
     * @param endSequence is a sequence state at a child node
     */
    public double getSequenceTransitionProbability(final List<Integer> startSequence, final List<Integer> endSequence, double distance, double lossDistance, int arrayLength) {

        List<Integer> startState = new ArrayList(startSequence);
        List<Integer> endState = new ArrayList(endSequence);


      if(startState.equals(lostState)) {
            //the lost state is an absorbing state P(lost -> lost) is 1.0
            if(endState.equals(lostState)) {
                return 1.0;
            }
            //see above P(lost -> WC) and P(lost -> normal state) is 0.0
            else {
                return 0.0;
            }

      }
      else  if(startState.equals(missingState)) {
          //P(WC state -> lost)
          if(endState.equals(lostState)) {
              return 1 - Math.exp(-missingRate * lossDistance);
          }
          //P(WC state -> WC state)
          else {
              return Math.exp(-missingRate * lossDistance);
          }
          //(WC state -> normal doesn't arise through ancestral state reconstruction)
      }

        else {
            //normal state (neither lost nor missing):
            if(endState.equals(missingState)) {
                //P(normal -> WC)
                return Math.exp(-missingRate * lossDistance);
            }

            else if(endState.equals(lostState)) {
                //P(normal -> lost)
                return 1.0 - Math.exp(- missingRate * lossDistance );
            }

            else {
                //create an unedited state to subtract from sequences to get only edited sites
                List<Integer> zero = Arrays.asList(0);

                //removing all unedited sites from each sequence
                startState.removeAll(zero);
                endState.removeAll(zero);

                //if endState is less edited than the start state, violates ordering
                if (startState.size() > endState.size()) {
                    return 0.0;
                }

                //subtracting start sequence from end sequence: edits introduced
                // if start state has identical elements to end state remove
                startState.forEach(endState::remove);
                List<Integer> newInserts = endState;

                //available positions are targetBClength length - number of edited positions
                int nrOfPossibleInserts = arrayLength - startState.size();

                //initialise the poisson distribution with mean rate * distance
                org.apache.commons.math3.distribution.PoissonDistribution poissonDistribution = new PoissonDistribution(distance);

                //calculate the transition probability for the case where all available positions are edited in
                // This is the absorbing state in the poisson process
                // P(max) = 1- sum(P(n)) * probability of this insert combination

                //The probability of going to the absorbing state is P(barcode not going missing AND the absorbing state being reached)
                if (newInserts.size() == nrOfPossibleInserts) {

                    return (Math.exp(-missingRate * lossDistance)) * calculateAbsorbingStateProbability(poissonDistribution, nrOfPossibleInserts) * combinedInsertProbabilities(newInserts);
                }
                //calculate the transition probability for the case where a #edits < available positions
                //this is a regular draw from the poisson process * probability of this insert combination

                //The probability of going to this edited state is P(barcode not going missing AND this state being reached)
                else if (newInserts.size() < nrOfPossibleInserts) {

                    return (Math.exp(-missingRate * lossDistance)) * poissonDistribution.probability(newInserts.size()) * combinedInsertProbabilities(newInserts);

                } else {

                    throw new RuntimeException("Error! Number of new inserts is larger than nr of possible inserts!");
                }
            }

        }
    }


    /**
     * This function calculates the probability of transitioning between 2 sequences states in given evolutionary time (distance)
     * (with potentially multiple edits having happened),on on tip edges of a tree
     *
     * @param startSequence  is a sequence state at a parent node
     * @param endSequence is a sequence state at a child node
     */
    public double getSequenceTransitionProbabilityTipEdge(final List<Integer> startSequence, final List<Integer> endSequence, double distance, double lossDistance,int arrayLength) {

        List<Integer> startState = new ArrayList(startSequence);
        List<Integer> endState = new ArrayList(endSequence);


        if(startState.equals(lostState)) {

            if(endState.equals(missingState)) {
                //P(lost -> dropout at tips)
                return 1.0;
            }
             if (endState.equals(lostState)) {
                 //P(lost -> lost at tips), this only arises as an artifact of assigning ancestral states at tips (if missing only missingState is valid at tips)
                 //a lost state at tips is undefined , only dropout
                return 0.0;
            }

        }

        //start state in a wildcard at the internal node
        else if(startState.equals(missingState)) {

            if(endState.equals(missingState)) {
                //P(WC -> dropout)
                // a WC state will be missing at tips by getting lost through heritable loss or by dropping out at sequencing
                return  (1 - Math.exp(- missingRate * lossDistance)) + ( Math.exp(- missingRate * lossDistance ) * missingProbability);
                //used to be: (1 - Math.exp(- missingRate * distance )) + above
            }

            else if(endState.equals(lostState)) {
                //P(WC -> lost at tips), this only arises as an artifact of assigning ancestral states at tips (if missing only missingState is valid at tips)
                //a lost state at tips is undefined , only dropout
                return 0.0 ;

            }


        }

        else {
            //normal start state
            if(endState.equals(missingState)) {
                //P(normal -> dropout)
                // a normal state will be missing at tips by getting lost through heritable loss or by dropping out at sequencing
                return  (1 - Math.exp(- missingRate * lossDistance)) + ( Math.exp(- missingRate * lossDistance ) * missingProbability);

            }

             if(endState.equals(lostState)) {
                 //P(normal -> lost at tips), this only arises as an artifact of assigning ancestral states at tips (if missing only missingState is valid at tips)
                 //a lost state at tips is undefined , only dropout
                 return 0.0;

             }

            //create an unedited state to subtract from sequences to get only edited sites
            List<Integer> zero = Arrays.asList(0);

            //removing all unedited sites from each sequence
            startState.removeAll(zero);
            endState.removeAll(zero);

            //if endState is less edited than the start state, violates ordering
            if (startState.size() > endState.size()) {
                return 0.0;
            }

            //subtracting start sequence from end sequence: edits introduced
            // if start state has identical elements to end state remove
            startState.forEach(endState::remove);
            List<Integer> newInserts = endState;

            //available positions are targetBClength length - number of edited positions
            int nrOfPossibleInserts = arrayLength - startState.size();

            //initialise the poisson distribution with mean rate * distance
            org.apache.commons.math3.distribution.PoissonDistribution poissonDistribution = new PoissonDistribution(distance);

            //calculate the transition probability for the case where all available positions are edited in
            // This is the absorbing state in the poisson process
            // P(max) = 1- sum(P(n)) * probability of this insert combination

            //The probability of going to the absorbing state is P(barcode not going missing AND the absorbing state being reached)
            if (newInserts.size() == nrOfPossibleInserts) {
                //used to have the proba
                return (1 - missingProbability) * ( Math.exp(- missingRate * lossDistance)) * calculateAbsorbingStateProbability(poissonDistribution, nrOfPossibleInserts) * combinedInsertProbabilities(newInserts);
            }
            //calculate the transition probability for the case where a #edits < available positions
            //this is a regular draw from the poisson process * probability of this insert combination

            //The probability of going to this edited state is P(barcode not going missing AND this state being reached)
            else if (newInserts.size() < nrOfPossibleInserts) {

                return (1 - missingProbability) * ( Math.exp(- missingRate * lossDistance)) * poissonDistribution.probability(newInserts.size()) * combinedInsertProbabilities(newInserts);

            } else {

                throw new RuntimeException("Error! Number of new inserts is larger than nr of possible inserts!");
            }

        }
        throw new RuntimeException("Error! a condition is missing somehow");

    }



    /**
     * This function calculates the probability of reaching/editing the last unedited position in the barcode with nbrOfPossibleInserts
     * available positions
     *
     * @param nbrOfPossibleInserts  is the number of available positions until the absorbing state is reached
     * @param dist is the poisson
     */
    public double calculateAbsorbingStateProbability(org.apache.commons.math3.distribution.PoissonDistribution dist,int nbrOfPossibleInserts) {

        double absorbingStateProbability = 1.0;

        for(int i = 0;  i < nbrOfPossibleInserts ; i++) {
            absorbingStateProbability -= dist.probability(i);
        }

        return absorbingStateProbability;

    }


    /**
     * Function to obtain the probability factor induced by insert frequencies
     * combineInsertProbabilities
     */
    public double combinedInsertProbabilities(List<Integer> inserts) {

        double factor = 1.0;
        double[] editProbabilitiesValue = editProbs;

        for(Integer i : inserts){
            //inserts are in {1, ..., nInserts}; insertProbabilities are in {0, ..., nInserts - 1}
            factor = factor * editProbabilitiesValue[i-1];
        }
        return factor;

    }

    /**
     * Function to obtain the array of insert probabilities
     *
     */
    public double[] getInsertProbabilities() {
        return this.editProbs;
    }

    /**
     * Function to obtain the dropout probability
     *
     */
    public double getMissingProbability() {
        return missingProbability;
    }

    /**
     * Function to obtain heritable loss rate
     *
     */
    public double getMissingRate() {
        return missingRate;
    }



    /**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        for (int i = 0; i < editProbabilities.size(); i++) {
            editProbs[i] = editProbabilities.get(i);
        }
        missingProbability = missProbability.get();
        missingRate = missRate.get();
        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {
        for (int i = 0; i < editProbabilities.size(); i++) {
            editProbs[i] = editProbabilities.get(i);
        }
        missingProbability = missProbability.get();
        missingRate = missRate.get();
        super.restore();

    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        //todo check better ways to do so
        for (int i = 0; i < editProbabilities.size(); i++) {
            editProbs[i] = editProbabilities.get(i);
        }
        missingProbability = missProbability.get();
        missingRate = missRate.get();
        return true;
    }

    /**
     * TODO: remove this (may require inheriting from another superclass)
     * methods that we do not need.
     * this not used in the current implementation because we do not use a matrix format
     */
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        if (dataType instanceof StandardData || dataType instanceof Binary || dataType instanceof IntegerData) {
            return true;
        }
        return false;
    }

}
