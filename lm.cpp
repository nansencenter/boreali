#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include Armadillo: C++ linear algebra library 
//http://arma.sourceforge.net/
//#include <armadillo>

//include CMINPACK for Levenberg-Marquardt optimization
//http://devernay.free.fr/hacks/cminpack/index.html
#include <cminpack.h>
//#define real __cminpack_real__
#define real long double

using namespace std;
//using namespace arma;

#include "lm.h"

Hydrooptics :: Hydrooptics(int inBands, double * inModel){
    // set number of bands and model.
    // number of bands
    bands = inBands;
    // hydro-optical model
    set_model(inModel);
    // measured reflectance
    s = new double [bands];
};

Hydrooptics :: ~Hydrooptics(){
    delete [] aaw;
    delete [] bbw;
    delete [] bw;
    delete [] aam;
    delete [] bbm;
    delete [] bm;
    delete [] s;

}

int Hydrooptics :: set_model(double * model){
    // Copy absorption and backscattering coefficients from input model
    // into the object

    int bn;
    
    aaw = new double [bands];
    bbw = new double [bands];
    bw  = new double [bands];
    aam = new double [bands * 3];
    bbm = new double [bands * 3];
    bm  = new double [bands * 3];
    
    for (bn = 0; bn < bands; bn ++){
        aaw[bn] = model[bn + 0 * bands];
        aam[bn + 0 * bands] = model[bn + 1 * bands];
        aam[bn + 1 * bands] = model[bn + 2 * bands];
        aam[bn + 2 * bands] = model[bn + 3 * bands];
        bbw[bn] = model[bn + 4 * bands];
        bw[bn]  = model[bn + 4 * bands] / B_WAT;
        bbm[bn + 0 * bands] = model[bn + 5 * bands];
        bbm[bn + 1 * bands] = model[bn + 6 * bands];
        bbm[bn + 2 * bands] = model[bn + 7 * bands];
        // TODO:
        // That has to be fixed later with variable backscattering
        // efficiency
        bm[bn + 0 * bands] = model[bn + 5 * bands] / B_CHL;
        bm[bn + 1 * bands] = model[bn + 6 * bands] / B_TSM;
        bm[bn + 1 * bands] = 0;
    }
        
    return 0;
}

int Hydrooptics :: set_params(double * inS, double inTheta){
    //set hydro-optical conditions
    int bn;
    
    for (bn = 0; bn < bands; bn ++){
        //measured reflectance
        s[bn] = inS[bn];
    }

    set_params(inTheta);
}

int Hydrooptics :: set_params(double inTheta){
    //set hydro-optical conditions

    // sun zenith
    theta = inTheta * PI / 180.;
    //inwater sun zenith
    double theta1 = asin(sin(theta) / NW);
    //inwater cos(sun zenith)
    double mu01 = cos(theta1);
    //internal surface reflectivity
    double rhoInt = 0.271 + 0.249 * mu01;
    //Calculate Q/f factor (from theta, from the object)
    //Q/f - conversion factor: Rrsw = R / Q/f
	qf = PI / (1 - rhoInt);
    //cos(sun zenith)
    mu0 = cos(theta);

    return 0;
}

double Hydrooptics :: rrsw (const double * c, int bn) const{
    // calculate Rrsw from given C, for a given band
    double a, bb, b, kd, r;
    
    // deep
    a  = aaw[bn] + aam[bn + 0 * bands] * c[0] + aam[bn + 1 * bands] * c[1] + aam[bn + 2 * bands] * c[2];
    bb = bbw[bn] + bbm[bn + 0 * bands] * c[0] + bbm[bn + 1 * bands] * c[1] + bbm[bn + 2 * bands] * c[2];
    r  = RW0 + RW1 * bb / a + RW2 * (bb * bb) / (a * a);
    
    return r;
};

double Hydrooptics :: rs (const double * c, int bn) const{
    // calcucate cost function per band
    double rs;

    rs = rrsw(c, bn) - s[bn];
    //maybe rs ^ 2 ?
    
    if (c[0] < 0 || c[1] < 0 || c[2] < 0)
        rs = 100;

    return rs;
};

double  Hydrooptics :: sse (const double * c) const{
    // calculate sum square error of all bands
    double sse = 0;
    int bn;

    for (bn = 0; bn < bands; bn ++)
        sse += pow(rs(c, bn),2);

    return sse;
};

double Hydrooptics :: j_deep(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT) const{


    return (RW1*bb0)/(aWAT+a0*c0+a1*c1+a2*c2)-RW2*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,3.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)*2.0-RW1*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)+RW2*bb0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)*2.0;
};

double Hydrooptics :: jacobian(const double * c, int bn, int vn) const{
    // calcluate jacobian for a given concentration/band/variable
    
    int v0n[3] = { 0, 1, 2 };
    int v1n[3] = { 1, 0, 0 };
    int v2n[3] = { 2, 2, 1 };

    double c0, c1, c2, a0, a1, a2, bb0, bb1, bb2, b0, b1, b2, aWAT, bbWAT, bWAT;
    double j;
    
    aWAT = aaw[bn];
    bbWAT = bbw[bn];
    bWAT = bw[bn];

    c0 = c[v0n[vn]];
    c1 = c[v1n[vn]];
    c2 = c[v2n[vn]];
    
    // define aw, a, bb, b in constructor
    a0 = aam[bn + v0n[vn]*bands];
    a1 = aam[bn + v1n[vn]*bands];
    a2 = aam[bn + v2n[vn]*bands];

    bb0 = bbm[bn + v0n[vn]*bands];
    bb1 = bbm[bn + v1n[vn]*bands];
    bb2 = bbm[bn + v2n[vn]*bands];

    b0 = bm[bn + v0n[vn]*bands];
    b1 = bm[bn + v1n[vn]*bands];
    b2 = bm[bn + v2n[vn]*bands];

    j = j_deep(s[bn],
        c0, c1, c2,
        a0, a1, a2,
        bb0, bb1, bb2,
        b0, b1, b2,
        aWAT, bbWAT, bWAT);

    return j;
};

int startingCPA(double parameters[6], double * startC){
    int ci0, ci1, ci2;
    double c0, c1, c2, dc0, dc1, dc2;
    double starts = 10;
    int k;
    int fullSize = starts * starts * starts;
    
    double min0 = parameters[0], max0 = parameters[1];
    double min1 = parameters[2], max1 = parameters[3];
    double min2 = parameters[4], max2 = parameters[5];
    
    dc0 = (max0 - min0) / (starts - 1);
    dc1 = (max1 - min1) / (starts - 1);
    dc2 = (max2 - min2) / (starts - 1);
    
    k = 0;
    c0 = min0;
    for (ci0 = 0; ci0 < starts; ci0 ++){
        c1 = min1;
        for (ci1 = 0; ci1 < starts; ci1 ++){
            c2 = min2;
            for (ci2 = 0; ci2 < starts; ci2 ++){
                startC[k * 3 + 0] = c0;
                startC[k * 3 + 1] = c1;
                startC[k * 3 + 2] = c2;
                k ++;
                c2 += dc2;
            }
            c1 += dc1;
        }
        c0 += dc0;
    }

    return fullSize;
}

//iterface to python
//calculate Rrsw from given concentrations, solar zenith
extern int get_rrsw(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double theta,
         double *outR, int outR_n0){

    int bn;

    //init HO-object
    Hydrooptics ho(outR_n0, model);
    //set ho-conditions
    ho.set_params(theta);

    for (bn = 0; bn < outR_n0; bn ++){
        //use HO-object for estimation of Rrsw
        outR[bn] = ho.rrsw(inC, bn);
    };

    return 0;
};

int compare (const void * v1, const void * v2)
{
    const double d1 = **(const double **)v1;
    const double d2 = **(const double **)v2;
    
    return d1<d2?-1:(d1>d2);
}

extern int get_c(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *theta, int theta_rows,
         double *outC, int outC_length){

    int bn, i, j, k, kBest, bands = inR_cols, pixels = inR_rows;
    printf("Retrieval from %d bands x %d pixels (W/O ARMA)...\n", bands, pixels);
    
    //result of optimization
    double xBest[4];
    
    //init HO-object
    Hydrooptics ho(bands, model);
    
    //create starting CPA for max 10x10x10 input vectors
    double sse[1000];
    double * ssep[1000];
    double startC[3000];
    // number of all starting vecotrs
    int startCN = startingCPA(parameters, startC);
    // number of best starting vecortrs
    int startBestCN = 10;
    
    //prepare for optimization with CMINPACK
    real x[5], fvec[10], fjac[30], tol, wa[300], fnorm;
    int info, ipvt[3], lwa = 100;
    //set tolerance to square of the machine recision
    tol = sqrt(__cminpack_func__(dpmpar)(1));

    //for all pixels
    for (i = 0; i < pixels; i ++){

        //set measured Rrsw and solar zenith in the HO-object
        ho.set_params(inR + i*bands, theta[i]);
        
        //erase xBest
        for (j = 0; j < 4; j ++)
            xBest[j] = 100;

        
        //estimate SSE for all starting vectors
        for (k = 0; k < startCN; k ++){
            
            // set initial concentrations
            //printf("start [%d]: ", k);
            for (j = 0; j < 3; j ++){
                x[j] = startC[k * 3 + j];
                //printf("%5.2g ", (double)x[j]);
            };
            
            // keep SSE and pointer to SSE
            sse[k] = ho.sse(x);
            ssep[k] = &sse[k];
            //printf("orig: %d %f %u\n", k, sse[k], ssep[k]);
        }
        
        //sort SSE pointers
        qsort(ssep, startCN, sizeof *ssep, compare);
        
        // start optimization from N best starting vectors
        for (k = 0; k < startBestCN; k ++){
            //get index of the k-th best starting vector
            kBest = ssep[k]-sse;
            
            // set initial concentrations
            //printf("start [%d]: ", kBest);
            for (j = 0; j < 3; j ++){
                x[j] = startC[kBest * 3 + j];
                //printf("%5.2g ", (double)x[j]);
            };

            //perform optimization
            info = __cminpack_func__(lmder1)(fcn, &ho, bands, 3, x, fvec, fjac, bands, tol, ipvt, wa, lwa);
            //estimate norm of residuals
            fnorm = __cminpack_func__(enorm)(bands, fvec);
    
            //printf(" ==> ");
            //for (j = 0; j < 3; j ++)
            //    printf("%5.2g ", (double)x[j]);
            //printf("%7.4g\n", (double)fnorm);
            
            if (fnorm < xBest[3]){
                for (j = 0; j < 3; j ++)
                    xBest[j] = (double)x[j];
                xBest[3] = (double)fnorm;
            }
        }
        
        //sent values of concentrtions and residuals back to Python
        //printf("final result: %d ", i);
        for (j = 0; j < 4; j ++){
            //printf("%g ", xBest[j]);
            outC[j + i*4] = xBest[j];
        }
        //printf("\n");

        if (fmod(i, 100.) == 0.)
            printf("%d\n", i);
    }

    printf("OK!\n");

    return 0;
};


int fcn(void *p, int m, int n, const real *x, real *fvec, real *fjac, 
	 int ldfjac, int iflag){

    int bn, i1;
    const Hydrooptics *ho = (Hydrooptics *)p;
    double rsval;
    
    if (iflag == 1){

        // calculate cost function
        for (bn = 0; bn < m; bn ++){
            rsval = ho -> rs(x, bn);
            fvec[bn] = rsval;
            //printf("fcn: fvec:%d %f\n", i0, fvec[i0]);
        };

    } else {

        //printf("fcn: j:\n");
        //cout << j << endl;
        //calculate jacobians
        for (bn = 0; bn < m; bn ++){
            fjac[bn + ldfjac*0] = ho->jacobian(x, bn, 0);
            fjac[bn + ldfjac*1] = ho->jacobian(x, bn, 1);
            fjac[bn + ldfjac*2] = ho->jacobian(x, bn, 2);
            
            //printf("fcn: fjac:%d %f %f %f\n", i0, fjac[i0 + ldfjac*0], fjac[i0 + ldfjac*1], fjac[i0 + ldfjac*2]);
        };
    };
    
    return 0;
};
