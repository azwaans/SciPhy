open module sciphy {
    requires beast.pkgmgmt;
    requires beast.base;
    requires feast;
    requires beast.fx;
    requires commons.math3;

    exports sciphy.evolution.datatype;
    exports sciphy.evolution.likelihood;
    exports sciphy.evolution.simulation;
    exports sciphy.evolution.substitutionmodel;
    exports sciphy.evolution.tree;
    exports sciphy.util;
    exports sciphy.app.beauti;

    provides beast.base.core.BEASTInterface with
            sciphy.evolution.datatype.SciPhyData,
            sciphy.evolution.likelihood.SciPhyParsimony,
            sciphy.evolution.likelihood.SciPhyTreeLikelihood,
            sciphy.evolution.simulation.SimulatedSciPhyAlignment,
            sciphy.evolution.simulation.SimulatedSciPhyAlignmentHeritableMissingBcodes,
            sciphy.evolution.substitutionmodel.SciPhySubstitutionModel,
            sciphy.evolution.tree.startingTree,
            sciphy.util.LogSum,
            sciphy.util.NexusImporter,
            sciphy.util.NexusParser,
            sciphy.util.SiteRatesLogger;
}