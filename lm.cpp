#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include CMINPACK for Levenberg-Marquardt optimization
#include <cminpack.h>

using namespace std;

#include "lm.h"

/// ======================== Hydrooptics ===============================

Hydrooptics :: Hydrooptics(int inBands, double * inModel){
    // set number of bands and model.
    // number of bands
    bands = inBands;
    // hydro-optical model
    aaw = new double [bands];
    bbw = new double [bands];
    bw  = new double [bands];
    aam = new double [bands * 3];
    bbm = new double [bands * 3];
    bm  = new double [bands * 3];
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

int Hydrooptics :: set_theta(double inTheta){
    //set gemetrical-optical conditions
    //used when calculating Rrsw and C
    int bn;
    
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

int Hydrooptics :: set_s(double * inS){
    //set values of measured remote sensing reflectance
    //used when calculating C
    int bn;
    
    for (bn = 0; bn < bands; bn ++){
        //measured reflectance
        s[bn] = inS[bn];
    }
}

double Hydrooptics :: rrsw (const double * c, int bn){
    // calculate Rrsw from given C, for a given band
    double a, bb, r;

    // deep
    a = aaw[bn] + aam[bn + 0 * bands] * c[0] + aam[bn + 1 * bands] * c[1] + aam[bn + 2 * bands] * c[2];
    bb = bbw[bn] + bbm[bn + 0 * bands] * c[0] + bbm[bn + 1 * bands] * c[1] + bbm[bn + 2 * bands] * c[2];
    r = RW0 + RW1 * bb / a + RW2 * (bb * bb) / (a * a);

    return r;
};

double Hydrooptics :: rs (const double * c, int bn){
    // calcucate cost function per band
    double rs;
    rs = rrsw(c, bn) - s[bn];

    if (c[0] < 0 || c[1] < 0 || c[2] < 0)
        rs = 100;

    return rs;
};

double  Hydrooptics :: sse (const double * c){
    // calculate sum square error of all bands
    double sse = 0;
    int bn;

    for (bn = 0; bn < bands; bn ++)
        sse += pow(rs(c, bn), 2);

    return sse;
};

int Hydrooptics :: retrieve(int cpas, int startCN, double * startC, double * xBest){
    // retreive concentrations

    bool log1 = false;
    bool log2 = false;
    
    int startBestCN = 10;
    int j, k, kBest;
    
    //data for best sse
    double ssev[4000];
    double * ssep[4000];

    //data for optimization with CMINPACK
    double x[6], fvec[10], fjac[30], tol, wa[300], fnorm;
    int info, ipvt[3], lwa = 100;
    //set tolerance to square of the machine recision
    tol = sqrt(dpmpar(1));
    
    //erase xBest
    for (j = 0; j < (cpas+1); j ++)
        xBest[j] = 100;

    //estimate SSE for all starting vectors
    for (k = 0; k < startCN; k ++){
        
        // set initial concentrations
        if (log1 and log2) printf("start [%d]: ", k);
        for (j = 0; j < cpas; j ++){
            x[j] = startC[k * cpas + j];
            if (log1 and log2) printf("%5.2g ", (double)x[j]);
        };
        
        // keep SSE value and pointer to SSE
        ssev[k] = sse(x);
        ssep[k] = &ssev[k];
        if (log1 and log2) printf("orig: %f %u\n", ssev[k], ssep[k]);
    }
    
    //sort SSE pointers
    qsort(ssep, startCN, sizeof *ssep, compare);
    
    // start optimization from 5 best starting vectors
    for (k = 0; k < startBestCN; k ++){
        //get index of the k-th best starting vector
        kBest = ssep[k]-ssev;
        
        // set initial concentrations
        if (log1) printf("start [%d / %d / %f]: ", k, kBest, ssev[kBest]);
        for (j = 0; j < cpas; j ++){
            x[j] = startC[kBest * cpas + j];
            if (log1) printf("%5.2g ", (double)x[j]);
        };

        //perform optimization using CMINPACK
        info = lmdif1(fcn, this, bands, cpas, x, fvec, tol, ipvt, wa, lwa);
        
        //estimate norm of residuals
        fnorm = enorm(bands, fvec);

        if (log1) {
            printf(" ==> ");
            for (j = 0; j < cpas; j ++)
                printf("%5.2g ", (double)x[j]);
            printf("%7.4g / %7.4g\n", (double)fnorm, sse(x));
        }
        
        if (fnorm < xBest[cpas]){
            for (j = 0; j < cpas; j ++)
                xBest[j] = (double)x[j];
            xBest[cpas] = (double)fnorm;
        }
    }
    //send values of concentrtions and residuals back to Python
    if (log1) {
        printf("final result:");
        for (j = 0; j < (cpas+1); j ++){
            printf("%g ", xBest[j]);
        }
        printf("\n");
    }
    
    return 0;
}

///======================== HydroopticsShallow ================================

HydroopticsShallow :: HydroopticsShallow(int inBands, double * inModel) : Hydrooptics(inBands, inModel){
    // albedo
    al = new double [bands];
};

HydroopticsShallow :: ~HydroopticsShallow(){
    delete [] al;
}

int HydroopticsShallow :: set_albedo(double * inAL){
    //set bottom spectral albedo
    for (int bn = 0; bn < bands; bn ++) al[bn] = inAL[bn];
}

double HydroopticsShallow :: rrsw (const double * c, int bn){
    // calculate Rrsw from given C, for a given band
    double a, bb, b, kd, r;

    // deep
    a = aaw[bn] + aam[bn + 0 * bands] * c[0] + aam[bn + 1 * bands] * c[1] + aam[bn + 2 * bands] * c[2];
    bb = bbw[bn] + bbm[bn + 0 * bands] * c[0] + bbm[bn + 1 * bands] * c[1] + bbm[bn + 2 * bands] * c[2];
    r  = RW0 + RW1 * bb / a + RW2 * (bb * bb) / (a * a);
    
    // shallow
    b  = bw[bn] + bm[bn + 0 * bands] * c[0] + bbm[bn + 1 * bands] * c[1] + bbm[bn + 2 * bands] * c[2] ;
    kd = sqrt(a * a + a * bb * (KD0 + KD1 * mu0)) / mu0;
    r = r * (1 - exp(-2 * h * kd)) + al[bn] * exp(-2 * h * kd) / qf;

    return r;
};

/// ============================== HydroopticsAlbedo ==========================

HydroopticsAlbedo :: HydroopticsAlbedo(int inBands, double * inModel, double * inLambda) : HydroopticsShallow(inBands, inModel){
    int bn;
    //wavelenthgs
    lambda = new double [bands];
    for (bn = 0; bn < bands; bn ++){
        lambda[bn] = inLambda[bn];
    }
};

HydroopticsAlbedo :: ~HydroopticsAlbedo(){
    delete [] lambda;
}

double HydroopticsAlbedo :: albedo(double al1, double al2, int bn){
    double alval;
    // calculate albedo from parametrization
    // = -0.06 + ASLP * 0.0005 * L
    //= AHMP * 0.25 * exp((555 - L)^2 / (2 * 40^2))
    alval = L0 + lambda[bn] * (al2 * 0.0005)
               + (al1 * 0.25) * exp(- pow(G1 - lambda[bn], 2.0) / 2 / G2 / G2);
    return alval;
}

double HydroopticsAlbedo :: rrsw (const double * c, int bn){
    // set albedo from parametrization
    al[bn] = albedo(c[3], c[4], bn);

    // calculate Rrsw for shallow waters
    return HydroopticsShallow :: rrsw(c, bn);
};


double HydroopticsAlbedo :: rs (const double * c, int bn){
    double rs;
    
    rs = rrsw(c, bn) - s[bn];

    if (c[0] < 0  || c[1] < 0   || c[2] < 0 ||
        c[3] < 0. || c[3] > 1.0 ||
        c[4] < 0. || c[4] > 1.0)
        rs = 100;

    return rs;
};

///============================================================================
int startingCPA(double parameters[6], double * startC){
    int ci0, ci1, ci2;
    double c0, c1, c2, dc0, dc1, dc2;
    double starts = 7.;
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


int startingCPA_al(double parameters[6], double * startC){
    int ci0, ci1, ci2, ci3, ci4;
    double c0, c1, c2, c3, c4, dc0, dc1, dc2, dc3, dc4;
    double starts = 5.;
    int k;
    int fullSize = starts * starts * starts * starts * starts;
    
    double min0 = parameters[0], max0 = parameters[1];
    double min1 = parameters[2], max1 = parameters[3];
    double min2 = parameters[4], max2 = parameters[5];
    double min3 = 0.0, max3 = 1.0;
    double min4 = 0.0, max4 = 1.0;
    
    dc0 = (max0 - min0) / (starts - 1);
    dc1 = (max1 - min1) / (starts - 1);
    dc2 = (max2 - min2) / (starts - 1);
    dc3 = (max3 - min3) / (starts - 1);
    dc4 = (max4 - min4) / (starts - 1);
    
    k = 0;
    c0 = min0;
    for (ci0 = 0; ci0 < starts; ci0 ++){
        c1 = min1;
        for (ci1 = 0; ci1 < starts; ci1 ++){
            c2 = min2;
            for (ci2 = 0; ci2 < starts; ci2 ++){
                c3 = min3;
                for (ci3 = 0; ci3 < starts; ci3 ++){
                    c4 = min4;
                    for (ci4 = 0; ci4 < starts; ci4 ++){
        
                        startC[k * 5 + 0] = c0;
                        startC[k * 5 + 1] = c1;
                        startC[k * 5 + 2] = c2;
                        startC[k * 5 + 3] = c3;
                        startC[k * 5 + 4] = c4;
        
                        k ++;
                        c4 += dc4;
                    }
                    c3 += dc3;
                }
                c2 += dc2;
            }
            c1 += dc1;
        }
        c0 += dc0;
    }

    return fullSize;
}


//iterface to python
//calculate Rrsw from given concentrations, solar zenith in deep waters
extern int get_rrsw_deep(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double theta,
         double *outR, int outR_n0){

    int bn;

    //init HO-object
    Hydrooptics ho(outR_n0, model);
    //set ho-conditions
    ho.set_theta(theta);
    for (bn = 0; bn < outR_n0; bn ++){
        //use HO-object for estimation of Rrsw
        outR[bn] = ho.rrsw(inC, bn);
    };

    return 0;
};


//iterface to python
//calculate Rrsw from given concentrations, solar zenith, depth, albedo
extern int get_rrsw_shal(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double theta,
         double h,
         double *albedo, int albedo_n0,
         double *outR, int outR_n0){

    int bn;

    //init HO-object
    HydroopticsShallow ho(outR_n0, model);
    //set ho-conditions
    ho.set_theta(theta);
    ho.h = h;
    ho.set_albedo(albedo);

    for (bn = 0; bn < outR_n0; bn ++){
        //use HO-object for estimation of Rrsw
        outR[bn] = ho.rrsw(inC, bn);
    };

    return 0;
};

//iterface to python
//calculate Rrsw from given concentrations (including albedo parametrization), solar zenith, depth, wavelengths
extern int get_rrsw_albe(
         double *model, int model_n0, int model_n1,
         double inC[5],
         double theta,
         double h,
         double *lambda, int lambda_n0,
         double *outR, int outR_n0){

    int bn;
    
    //init HO-object
    HydroopticsAlbedo ho(outR_n0, model, lambda);

    //set ho-conditions
    ho.set_theta(theta);
    ho.h = h;

    for (bn = 0; bn < outR_n0; bn ++){
        //use HO-object for estimation of Rrsw
        outR[bn] = ho.rrsw(inC, bn);
    };

    return 0;

};

// compare values of two pointers (used for sorting starting vectors)
int compare (const void * v1, const void * v2)
{
    const double d1 = **(const double **)v1;
    const double d2 = **(const double **)v2;
    
    return d1<d2?-1:(d1>d2);
}

//iterface to python
//calculate C from given Rrsw, solar zenith in deep waters
extern int get_c_deep(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *theta, int theta_rows,
         double *outC, int outC_length){

    int i, j, bands = inR_cols, pixels = inR_rows, cpas=3;

    //init HO-object
    Hydrooptics ho(inR_cols, model);

    printf("Retrieval from %d bands x %d pixels...\n", bands, pixels);
    
    double startC[3000], xBest[6];
    // create starting vecotrs
    int startCN = startingCPA(parameters, startC);

    //for all pixels
    for (i = 0; i < pixels; i ++){

        // set hydro-optical conditions
        ho.set_theta(theta[i]);
        
        // set measured reflectance
        ho.set_s(inR + i*bands);

        // retrieve concentrations
        ho.retrieve(cpas, startCN, startC, xBest);

        //send values of concentrtions and residuals back to Python
        for (j = 0; j < (cpas+1); j ++){
            outC[j + i*(cpas+1)] = xBest[j];
        }
    
        if (fmod(i, 100.) == 0.)
            printf("%d / %f \n", i, (double)(i) / (double)(pixels));

    }

    printf("OK!\n");

    return 0;
};

//iterface to python
//calculate C from given Rrsw, solar zenith, depth, albedo
extern int get_c_shal(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *theta, int theta_rows,
         double *h, int h_rows,
         double *albedo, int albedo_rows, int albedo_cols,
         double *outC, int outC_length){

    int i, j, startCN, bands = inR_cols, pixels = inR_rows, cpas=3;
    double startC[3000], xBest[6];
    
    //init HO-object
    HydroopticsShallow ho(inR_cols, model);

    printf("Retrieval from %d bands x %d pixels...\n", bands, pixels);
    
    // create starting vecotrs
    startCN = startingCPA(parameters, startC);

    //for all pixels
    for (i = 0; i < pixels; i ++){

        // set hydro-optical conditions
        ho.set_theta(theta[i]);
        
        // set measured reflectance, albedo and depth
        ho.set_s(inR + i*bands);
        ho.set_albedo(albedo + i*bands);
        ho.h = h[i];

        // retrieve concentrations
        ho.retrieve(cpas, startCN, startC, xBest);

        //send values of concentrtions and residuals back to Python
        for (j = 0; j < (cpas+1); j ++){
            outC[j + i*(cpas+1)] = xBest[j];
        }
    
        if (fmod(i, 100.) == 0.)
            printf("%d / %f \n", i, (double)(i) / (double)(pixels));

    }

    printf("OK!\n");

    return 0;
}


//iterface to python
//calculate C (including albedo params) from given Rrsw, solar zenith, depth, wavelengths
extern int get_c_albe(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *theta, int theta_rows,
         double *h, int h_rows,
         double *lambda, int lambda_n0,
         double *outC, int outC_length){

    int i, j, startCN, bands = inR_cols, pixels = inR_rows, cpas=5;
    double startC[16000], xBest[6];
    
    //init HO-object
    HydroopticsAlbedo ho(inR_cols, model, lambda);

    printf("Retrieval from %d bands x %d pixels...\n", bands, pixels);
    
    // create starting vecotrs
    startCN = startingCPA_al(parameters, startC);
    printf("Starting CPA - OK!\n");

    //for all pixels
    for (i = 0; i < pixels; i ++){

        // set hydro-optical conditions
        ho.set_theta(theta[i]);
        
        // set measured reflectance and depth
        ho.set_s(inR + i*bands);
        ho.h = h[i];

        // retrieve concentrations
        ho.retrieve(cpas, startCN, startC, xBest);

        //send values of concentrtions and residuals back to Python
        for (j = 0; j < (cpas+1); j ++){
            outC[j + i*(cpas+1)] = xBest[j];
        }
    
        if (fmod(i, 100.) == 0.)
            printf("%d / %f \n", i, (double)(i) / (double)(pixels));
    }

    printf("OK!\n");

    return 0;
}

//interface to CMINPACK
// functions uses finite difference derivatives of Rrsw
int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag){

    int bn, vn;
    Hydrooptics *ho = (Hydrooptics *)p;
    double rsval, rssum;
    
        // calculate cost function
        rssum = 0;
        for (bn = 0; bn < m; bn ++){
            rsval = ho -> rs(x, bn);
            fvec[bn] = rsval;
            rssum += ho -> rs(x, bn);
        };
    
    return 0;
};
