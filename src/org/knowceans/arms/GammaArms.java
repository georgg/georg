/*
 * Created on Jul 31, 2005
 */
/*
 * (C) Copyright 2005, Gregor Heinrich (gregor :: arbylon : net) (This file is
 * part of the knowceans-arms (org.knowceans.arms.*) experimental software package.)
 */
/*
 * knowceans-arms is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 */
/*
 * knowceans-arms is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 */
/*
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.knowceans.arms;

import org.knowceans.util.Densities;
import org.knowceans.util.Histogram;

public class GammaArms extends ArmSampler {

    /**
     * implements the gamma density.
     * 
     * @param x
     * @param params double[]{alpha, beta}
     * @return
     */
    public double logpdf(double x, Object params) {
        double[] pars = (double[]) params;
        return Math.log(Densities.pdfGamma(x, pars[0], pars[1]));
    }

    public static void main(String[] args) {
        GammaArms gars = new GammaArms();
        double[] xprev = new double[] {0.3};
        double[] params = new double[] {40.0, 30.0};

        int nsamp = 1000;
        double[] samples = new double[nsamp];
        for (int i = 0; i < nsamp; i++) {
            try {
                samples[i] = gars.armsSimple(params, 4, new double[] {0.2},
                    new double[] {8000.}, true, xprev);
                System.out.println(samples[i]);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        Histogram.hist(System.out, samples, 100);
        
    }
}
