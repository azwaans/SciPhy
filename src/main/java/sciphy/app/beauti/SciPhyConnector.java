package sciphy.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.Alignment;
import beast.base.spec.inference.distribution.Dirichlet;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Simplex;
import beastfx.app.inputeditor.BeautiDoc;
import sciphy.evolution.likelihood.SciPhyTreeLikelihood;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

//Connector for SciPhy data:
// SciPhy data doesn't have a fixed/a priori-known character state count. Each alignment
// may contain different state counts, which must be extracted to set correct dimensions and starting values for
// state probabilities (editing probabilities in SciPhy terminology) and the corresponding priors (in the current template,
// dirichlet prior).
// The SciPhy likelihood requires specifying an alignment length per partition which is extracted from
// the SiteCount.
// The SciPhy tree likelihood usually is conditioned on an unedited sequence state at the origin of the process, which
// corresponds to the experimental length (in single-cell analyses, this is known).
public class SciPhyConnector {

    public static void customConnector(BeautiDoc doc) {
        for (BEASTInterface plugin : doc.pluginmap.values()) {
            if (plugin instanceof SciPhySubstitutionModel substModel) {

                //  Fetch the specific Alignment belonging to this partition
                Alignment alignment = getAlignmentForSubstModel(doc, substModel);
                if (alignment == null) continue;

                int targetDimension = getDimensionFromAlignment(alignment); // maxState
                int siteCount = alignment.getSiteCount();                   // sequence length for THIS partition

                if (targetDimension <= 0 || siteCount <= 0) continue;

                Simplex currentSimplex = substModel.editProbabilitiesInput.get();

                if (currentSimplex instanceof SimplexParam editProbsParam) {

                    // Retrieve related partition parameters
                    RealVectorParam alphaParam = findAlphaParam(doc, editProbsParam);
                    IntScalarParam arrayLengthParam = findArrayLengthParam(doc, editProbsParam);

                    // Skip if ALL parameters for this partition are already updated
                    boolean editProbsDone = (editProbsParam.size() == targetDimension);
                    boolean alphaDone = (alphaParam == null || alphaParam.size() == targetDimension);
                    boolean arrayLengthDone = (arrayLengthParam == null || isArrayLengthUpdated(arrayLengthParam, siteCount));

                    if (editProbsDone && alphaDone && arrayLengthDone) {
                        continue; // Prevents infinite refresh loops in BEAUti
                    }

                    // Update editProbabilities [1/maxState, ..., 1/maxState]
                    if (!editProbsDone) {
                        Double[] uniformValues = new Double[targetDimension];
                        Arrays.fill(uniformValues, 1.0 / targetDimension);

                        List<Double> valList = editProbsParam.valuesInput.get();
                        if (valList != null) {
                            valList.clear();
                            valList.addAll(Arrays.asList(uniformValues));
                        } else {
                            editProbsParam.valuesInput.setValue(new ArrayList<>(Arrays.asList(uniformValues)), editProbsParam);
                        }
                        editProbsParam.isEstimatedInput.setValue(true, editProbsParam);
                        editProbsParam.initAndValidate();
                    }

                    // Update Dirichlet alpha RealVectorParam [1.5, ..., 1.5]
                    if (alphaParam != null && !alphaDone) {
                        Double[] alphaValues = new Double[targetDimension];
                        Arrays.fill(alphaValues, 1.5);

                        Object currentVal = alphaParam.valuesInput.get();
                        if (currentVal instanceof List<?>) {
                            @SuppressWarnings("unchecked")
                            List<Double> alphaValList = (List<Double>) currentVal;
                            alphaValList.clear();
                            alphaValList.addAll(Arrays.asList(alphaValues));
                        } else {
                            alphaParam.valuesInput.setValue(new ArrayList<>(Arrays.asList(alphaValues)), alphaParam);
                        }
                        alphaParam.initAndValidate();
                    }

                    // Update arrayLength IntScalarParam to match this partition's site count
                    if (arrayLengthParam != null && !arrayLengthDone) {
                        arrayLengthParam.valuesInput.setValue(String.valueOf(siteCount), arrayLengthParam);
                        arrayLengthParam.initAndValidate();
                    }
                }
            }
        }
    }

    /**
     * Partition-aware lookup: locates the specific Alignment associated with the given substitution model.
     */
    private static Alignment getAlignmentForSubstModel(BeautiDoc doc, SciPhySubstitutionModel substModel) {
        if (substModel.getID() != null) {
            // Extracts "partitionName" from "SciPhySubstmodel.s:partitionName"
            String id = substModel.getID();
            String partitionName = id.contains(":") ? id.substring(id.indexOf(":") + 1) : id;

            for (Alignment alignment : doc.alignments) {
                if (alignment.getID() != null && alignment.getID().equals(partitionName)) {
                    return alignment;
                }
            }
        }

        // Fallback: search SciPhyTreeLikelihood structure for matching substitution model linkage
        for (BEASTInterface plugin : doc.pluginmap.values()) {
            if (plugin instanceof SciPhyTreeLikelihood likelihood) {
                if (likelihood.siteModelInput.get() instanceof BEASTInterface siteModel) {
                    for (var input : siteModel.listInputs()) {
                        if (input.get() == substModel) {
                            return likelihood.dataInput.get();
                        }
                    }
                }
            }
        }

        // Fallback for single-alignment setups
        if (doc.alignments != null && doc.alignments.size() == 1) {
            return doc.alignments.get(0);
        }

        return null;
    }

    /**
     * Checks if arrayLength parameter already equals the partition's site count.
     */
    private static boolean isArrayLengthUpdated(IntScalarParam arrayLengthParam, int siteCount) {
        if (arrayLengthParam == null) return true;

        try {
            if (arrayLengthParam != null && arrayLengthParam.get() == siteCount) {
                return true;
            }
        } catch (Exception ignored) {}

        Object val = arrayLengthParam.valuesInput.get();
        if (val != null) {
            String valStr = val.toString().trim();
            return valStr.equals(String.valueOf(siteCount));
        }
        return false;
    }

    /**
     * Finds IntScalarParam for arrayLength corresponding strictly to the current partition.
     */
    private static IntScalarParam findArrayLengthParam(BeautiDoc doc, SimplexParam editProbsParam) {
        if (editProbsParam.getID() != null) {
            String partitionName = editProbsParam.getID()
                    .replace("editProbabilities.s:", "")
                    .replace("editProbabilities.", "");

            String arrayLengthID = "arrayLength." + partitionName;
            BEASTInterface plugin = doc.pluginmap.get(arrayLengthID);
            if (plugin instanceof IntScalarParam param) {
                return param;
            }
        }
        return null;
    }

    /**
     * Scans alignment matrix for maximum state code (maxState).
     */
    public static int getDimensionFromAlignment(Alignment alignment) {
        int maxState = 0;
        for (int i = 0; i < alignment.getPatternCount(); i++) {
            for (int code : alignment.getPattern(i)) {
                if (code > maxState) {
                    maxState = code;
                }
            }
        }
        return maxState;
    }

    /**
     * Finds matching RealVectorParam for alpha strictly within the partition.
     */
    private static RealVectorParam findAlphaParam(BeautiDoc doc, SimplexParam editProbsParam) {
        if (editProbsParam.getID() != null) {
            String alphaID = editProbsParam.getID().replace("editProbabilities.", "alpha.");
            BEASTInterface plugin = doc.pluginmap.get(alphaID);
            if (plugin instanceof RealVectorParam param) {
                return param;
            }
        }

        for (BEASTInterface plugin : doc.pluginmap.values()) {
            if (plugin instanceof Dirichlet dirichlet) {
                if (dirichlet.paramInput.get() == editProbsParam) {
                    RealVector<?> alphaVector = dirichlet.alphaInput.get();
                    if (alphaVector instanceof RealVectorParam param) {
                        return param;
                    }
                }
            }
        }
        return null;
    }
}