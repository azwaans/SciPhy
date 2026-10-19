/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package sciphy.util;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.Logger;
import beast.base.inference.util.ESS;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.util.DiscreteStatistics;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Special logger for logging the Gamma Site Category rates
 *
 * @author azwaans 
 * based on 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SiteRatesLogger extends BEASTObject implements Loggable {

    public Input<SiteModel> siteModelInput = new Input<SiteModel>(
            "siteModel",
            "site tree whose rates to log.",
            Validate.REQUIRED);

    public Input<Tree> treeInput = new Input<Tree>(
            "tree",
            "tree is needed for getCategoryRates.",
            Validate.REQUIRED);



    private SiteModel siteModel;
    private Tree tree;



    @Override
    public void initAndValidate() {
        siteModel = siteModelInput.get();
        tree = treeInput.get();

    };

    @Override
    public void init(PrintStream out) {
        String outName;
        outName = "siteModel";
        for (int i=0; i<siteModel.getCategoryCount(); i++) {
            out.print(outName + ".rateCategory_" + i + "\t");
        }

    }

    @Override
    public void log(long nSample, PrintStream out) {

        for (int i = 0; i < siteModel.getCategoryRates(tree.getRoot()).length; i++) {
            out.print(siteModel.getCategoryRates(tree.getRoot())[i] + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
    }

//    /**
//     * Compute statistics from completed traces.
//     */
//    public void computeStatistics() {
//
//        // Truncate burnin
//        heights = heights.subList((int)(burninFrac*heights.size()),
//                heights.size()-1);
//
//        // Transfer to array for DiscreteStatistics methods:
//        heightsArray = new double[heights.size()];
//        for (int i=0; i<heights.size(); i++)
//            heightsArray[i] = heights.get(i);
//
//        // Compute height statistics:
//        heightMean = DiscreteStatistics.mean(heightsArray);
//        heightVar = DiscreteStatistics.variance(heightsArray);
//        heightESS = ESS.calcESS(heights);
//    }
//
//    public double getRateMean() {
//        return heightMean;
//    }
//
//    public double getHeightVar() {
//        return heightVar;
//    }
//
//    public double getHeightESS() {
//        return heightESS;
//    }
    
}
