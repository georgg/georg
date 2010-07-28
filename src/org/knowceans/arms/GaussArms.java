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

/**
 * GaussArms simulates a normal distribution with ARMS. This reimplements the
 * ARMS example 01.
 * 
 * @author heinrich
 */
public class GaussArms extends ArmSampler {

    public class GmmParams {
        double[] prob;

        double[] mean;

        double[] sigma;

        int k;
    }

    /**
     * log of Gaussian pdf.
     * 
     * @param x
     * @param params double[]{mean, sdev}
     */
    @Override
    public double logpdf(double x, Object params) {

        double[] a = (double[]) params;
        double mu = a[0];
        double sigma = a[1];
        return Math.log(Densities.pdfNorm(x, mu, sigma));
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
            BufferedWriter bw = new BufferedWriter(new FileWriter("arms.out01"));

            GaussArms ga = new GaussArms();

            /*
             * initialise data needed by normal density function
             */
            double mean = 10.;
            double sdev = 5.;

            for (i = 0; i < 10000; i++) {
                ga.arms(new double[] {mean, sdev}, xinit, ninit, xl, xr,
                    convex, npoint, dometrop, xprev, xsamp, nsamp, qcent,
                    xcent, ncent, neval);

                bw.write(i + " " + xsamp[0] + "  " + neval[0] + "\n");

                /*
                 * update xprev to get a Markov chain
                 */
                xprev[0] = xsamp[0];
            }

            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
