package sciphy.util;

import static java.lang.Math.log1p;

public class LogSum {

    /**
     * Calculates the log of the sum of a collection from the collection of log transformed values
     * without having to exponentiate all elements
     *
     * @param la array of log
     * @param la size
     */
    public static double logSum(double la[], int numElements) {
        // Assume index_of_max() finds the maximum element
        // in the array and returns its index
        double max = la[0];
        int index = 0;
        for (int i = 0; i < la.length; i++) {
            if (max < la[i]) {
                max = la[i];
                index = i;
            }
        }

        double sum_exp = 0;
        for (int i = 0; i < numElements; i++) {
            if (i == index) {
                continue;
            }
            sum_exp += Math.exp(la[i] - la[index]);
        }

        return la[index] + log1p(sum_exp);
    }
}
