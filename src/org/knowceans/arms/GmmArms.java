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

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import org.knowceans.util.Densities;
import org.knowceans.util.Histogram;

/**
 * GmmArms simulates a normal mixture. This reimplements the ARMS example 02.
 * 
 * @author heinrich
 */
public class GmmArms extends ArmSampler {

    public class GmmParams {
        double[] prob;

        double[] mean;

        double[] sigma;

        int k;
    }

    /**
     * log of Gaussian mixture pdf.
     * 
     * @param x
     * @param double[][]{{mean}, {sigma}, {pi}}
     */
    @Override
    public double logpdf(double x, Object params) {

        double[][] a = (double[][]) params;
        double[] mean = a[0];
        double[] sigma = a[1];
        double[] prob = a[2];
        return Math.log(Densities.pdfGmm(x, mean.length, prob, mean, sigma));
    }

    public static void main(String[] args) {
        int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4, i;
        int[] neval = new int[1];
        double[] xinit = {-3.0, 0.0, 12.0, 20.0}, xl = {-100.0}, xr = {100.0};
        double[] xsamp = new double[100];
        double[] xcent = new double[10];
        double[] qcent = new double[] {5., 30., 70., 95.};
        long seed = 44;
        double[] convex = {1.0};
        boolean dometrop = true;
        double[] xprev = {0.0};

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("arms.out02"));

            GmmArms ga = new GmmArms();

            /*
             * initialise data needed by normal density function
             */
            double[] mean = new double[] {5., 10.};
            double[] sigma = new double[] {1., 2.5};
            double[] prob = new double[] {.3, .7};
            int totsamp = 100000;
            double[] samples = new double[totsamp];
            double[] sample = new double[1]; 

            for (i = 0; i < totsamp; i++) {
                sample = ga.arms(new double[][] {mean, sigma, prob}, xinit, ninit, xl,
                    xr, convex, npoint, dometrop, xprev, xsamp, nsamp, qcent,
                    xcent, ncent, neval);
                samples[i] = sample[0];

                bw.write(i + " " + xsamp[0] + "  " + neval[0] + "\n");

                /*
                 * update xprev to get a Markov chain
                 */
                xprev[0] = xsamp[0];
            }
            bw.close();
            Histogram.hist(System.out, samples, 100);            
            
            
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
