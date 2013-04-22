#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include Armadillo: C++ linear algebra library 
//http://arma.sourceforge.net/
#include <armadillo>

//include CMINPACK for Levenberg-Marquardt optimization
//http://devernay.free.fr/hacks/cminpack/index.html
#include <cminpack.h>
//#define real __cminpack_real__
#define real long double

using namespace std;
using namespace arma;

#include "lm.h"

Hydrooptics :: Hydrooptics(mat inM, mat inAL, double inH, double inTheta){
    m = inM;
    al = inAL;
    h = inH;
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

};

Hydrooptics :: Hydrooptics(mat inM, double inH, double inTheta, mat inLambda){
    m = inM;
    lambda = inLambda;
    h = inH;
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

};

mat Hydrooptics :: cOne(mat c) const{
    mat cc = ones<mat>(1, c.n_cols+1);
    cc.cols(1, c.n_cols) = c;

    return cc;
};
    
mat Hydrooptics :: albedo(double al1, double al2) const{
    mat al = zeros<mat>(1, lambda.n_cols);
    int i;
    for (i = 0; i < lambda.n_cols; i ++)
        al(0, i) = L0+lambda(0,i)*(al2*0.0005)+(al1*0.25)*exp(1.0/(G2*G2)*pow(G1-lambda(0,i),2.0)*(-1.0/2.0));
        //printf("%f %f\n", lambda(0,i), al(0, i));
            
    return al;
}

mat Hydrooptics :: rrsw (mat c) const{
    //[1 c0 c1 c2]
    mat cc = cOne(c);
    //total absorption
    mat a = cc * m.rows(0, 3);
    //total back-scattering
    mat bb = cc * m.rows(4, 7);
    //total subsurface remote sensing reflectance for deep waters
    mat rrsw = RW0 + RW1 * bb / a + RW2 * (bb % bb) / (a % a);
    mat result = rrsw;
    
    if (h > 0){
        // specific scattering (for WATER, CHL and TSM)
        mat b = m.rows(4, 7); //.zeros();
        b.row(0) /= B_WAT;
        b.row(1) /= B_CHL;
        b.row(2) /= B_TSM;

        //total scattering
        b = cc * b;

        //volume attenuation coeff (Kd)
        mat kd = sqrt(a % a + a % b * (KD0 + KD1 * mu0)) / mu0;

        //total subsurface remote sensing reflectance for shallow waters
        result = rrsw % (1 - exp(-2 * h * kd)) + al % exp(-2 * h * kd) / qf;
    }

    return result;
};


mat Hydrooptics :: rrsw (mat c, double al1, double al2) const{
    //[1 c0 c1 c2]
    mat cc = cOne(c);
    //total absorption
    mat a = cc * m.rows(0, 3);
    //total back-scattering
    mat bb = cc * m.rows(4, 7);
    //total subsurface remote sensing reflectance for deep waters
    mat rrsw = RW0 + RW1 * bb / a + RW2 * (bb % bb) / (a % a);

    //albedo
    mat al = albedo(al1, al2);
    //cout << al << endl;
    
    // specific scattering (for WATER, CHL and TSM)
    mat b = m.rows(4, 7); //.zeros();
    b.row(0) /= B_WAT;
    b.row(1) /= B_CHL;
    b.row(2) /= B_TSM;

    //total scattering
    b = cc * b;

    //volume attenuation coeff (Kd)
    mat kd = sqrt(a % a + a % b * (KD0 + KD1 * mu0)) / mu0;

    //total subsurface remote sensing reflectance for shallow waters
    return rrsw % (1 - exp(-2 * h * kd)) + al % exp(-2 * h * kd) / qf;
};

mat Hydrooptics :: rs (mat c) const{
    int i;
    mat rs = rrsw(c) - s;
    if (c(0, 0) < 0 || c(0, 1) < 0 || c(0, 2) < 0)
        for (i = 0; i < rs.n_cols; i ++)
            rs(0, i) = 100;
    return rs;
};

mat Hydrooptics :: rs (mat c, double al1, double al2) const{
    int i;
    mat rs = rrsw(c, al1, al2) - s;
    if (c(0, 0) < 0. ||
        c(0, 1) < 0. ||
        c(0, 2) < 0. ||
        al1 < 0. ||
        al2 < 0.3 ||
        al1 > 1. ||
        al2 > 1.)
        for (i = 0; i < rs.n_cols; i ++)
            rs(0, i) = 100;
    return rs;
};

double Hydrooptics :: j_deep(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT) const{


    return (RW1*bb0)/(aWAT+a0*c0+a1*c1+a2*c2)-RW2*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,3.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)*2.0-RW1*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)+RW2*bb0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)*2.0;
};

double Hydrooptics :: j_shallow(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al) const{
    b0 = b0 == 0 ? 1 : b0;
    b1 = b1 == 0 ? 1 : b1;
    b2 = b2 == 0 ? 1 : b2;
    
    return -(exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)-1.0)*((RW1*bb0)/(aWAT+a0*c0+a1*c1+a2*c2)-RW2*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,3.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)*2.0-RW1*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)+RW2*bb0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)*2.0)+(h*exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*1.0/sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*(a0*(aWAT+a0*c0+a1*c1+a2*c2)*2.0+a0*(KD0+KD1*mu0)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2)+(bb0*(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2))/b0)*(RW0+(RW1*(bbWAT+bb0*c0+bb1*c1+bb2*c2))/(aWAT+a0*c0+a1*c1+a2*c2)+RW2*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)))/mu0-(al*h*exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*1.0/sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*(a0*(aWAT+a0*c0+a1*c1+a2*c2)*2.0+a0*(KD0+KD1*mu0)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2)+(bb0*(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2))/b0))/(mu0*qf);
};

//Calculate Jacobian of the cost function with additional variable: albedo
//Returns single value of derivative on C at one wavelength
double Hydrooptics :: j_shallow_c0(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const{
    b0 = b0 == 0 ? 1 : b0;
    b1 = b1 == 0 ? 1 : b1;
    b2 = b2 == 0 ? 1 : b2;
           
    return -(exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)-1.0)*((RW1*bb0)/(aWAT+a0*c0+a1*c1+a2*c2)-RW2*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,3.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)*2.0-RW1*a0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)+RW2*bb0*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*(bbWAT+bb0*c0+bb1*c1+bb2*c2)*2.0)+(h*exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*1.0/sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*(a0*(aWAT+a0*c0+a1*c1+a2*c2)*2.0+a0*(KD0+KD1*mu0)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2)+(bb0*(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2))/b0)*(RW0+(RW1*(bbWAT+bb0*c0+bb1*c1+bb2*c2))/(aWAT+a0*c0+a1*c1+a2*c2)+RW2*1.0/pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)*pow(bbWAT+bb0*c0+bb1*c1+bb2*c2,2.0)))/mu0-(h*exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*1.0/sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*(L0+al2*ll*(1.0/2.0E3)+al1*exp(1.0/(G2*G2)*pow(G1-ll,2.0)*(-1.0/2.0))*(1.0/4.0))*(a0*(aWAT+a0*c0+a1*c1+a2*c2)*2.0+a0*(KD0+KD1*mu0)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2)+(bb0*(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2))/b0))/(mu0*qf);
};

//Calculate Jacobian of the cost function with additional variable: albedo
//Returns single value of derivative on ALR (albedo ratio) at one wavelength
double Hydrooptics :: j_shallow_al1(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const{
    b0 = b0 == 0 ? 1 : b0;
    b1 = b1 == 0 ? 1 : b1;
    b2 = b2 == 0 ? 1 : b2;

    return (exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*exp(1.0/(G2*G2)*pow(G1-ll,2.0)*(-1.0/2.0))*(1.0/4.0))/qf;
};


//Calculate Jacobian of the cost function with additional variable: albedo
//Returns single value of derivative on ALR (albedo ratio) at one wavelength
double Hydrooptics :: j_shallow_al2(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const{
    b0 = b0 == 0 ? 1 : b0;
    b1 = b1 == 0 ? 1 : b1;
    b2 = b2 == 0 ? 1 : b2;

    return (ll*exp((h*sqrt(pow(aWAT+a0*c0+a1*c1+a2*c2,2.0)+(KD0+KD1*mu0)*(aWAT+a0*c0+a1*c1+a2*c2)*(bbWAT/bWAT+(bb0*c0)/b0+(bb1*c1)/b1+(bb2*c2)/b2))*-2.0)/mu0)*(1.0/2.0E3))/qf;
};

mat Hydrooptics :: jacobian(mat c) const{
    
    mat j = zeros<mat>(3, s.n_cols);
    
    int v0n[3] = { 0, 1, 2 };
    int v1n[3] = { 1, 0, 0 };
    int v2n[3] = { 2, 2, 1 };

    //cos(sun zenith)
	double mu0 = cos(PI * theta / 180.);

    double c0, c1, c2, a0, a1, a2, bb0, bb1, bb2, b0, b1, b2, aWAT, bbWAT, bWAT;
    double tmp;
    
    mat aw = m.row(0);
    mat bbw = m.row(4);
    mat bw = bbw / B_WAT;
    
    mat a = m.rows(1, 3);
    mat bb = m.rows(5, 7);

    // specific scattering (for CHL, TSM, DOC)
    mat b = zeros<mat>(bb.n_rows, bb.n_cols);
    b.row(0) = bb.row(0) / B_CHL;
    b.row(1) = bb.row(1) / B_TSM;
    b.row(2) = bb.row(2) / B_DOC;

    for (int bn = 0; bn < s.n_cols; bn ++){

        aWAT = aw(0, bn);
        bbWAT = bbw(0, bn);
        bWAT = bw(0, bn);

        for (int vn = 0; vn < 3; vn ++){

            c0 = c(v0n[vn]);
            c1 = c(v1n[vn]);
            c2 = c(v2n[vn]);
        
            a0 = a(v0n[vn], bn);
            a1 = a(v1n[vn], bn);
            a2 = a(v2n[vn], bn);
        
            bb0 = bb(v0n[vn], bn);
            bb1 = bb(v1n[vn], bn);
            bb2 = bb(v2n[vn], bn);
    
            b0 = b(v0n[vn], bn);
            b1 = b(v1n[vn], bn);
            b2 = b(v2n[vn], bn);
        
            if (h < 0) {

                j(vn, bn) = j_deep(s(0, bn),
                    c0, c1, c2,
                    a0, a1, a2,
                    bb0, bb1, bb2,
                    b0, b1, b2,
                    aWAT, bbWAT, bWAT);

            } else {

                j(vn, bn) = j_shallow(s(0, bn),
                    c0, c1, c2,
                    a0, a1, a2,
                    bb0, bb1, bb2,
                    b0, b1, b2,
                    aWAT, bbWAT, bWAT,
                    al(0, bn));
            }
        }
    }
    
    return j;
};


//Calculate Jacobian of the cost function:
//partial derivative of the cost function by each concentration and albedo
//Returns vector for the entrie spectrum
mat Hydrooptics :: jacobian(mat c, double al1, double al2) const{

    mat j = zeros<mat>(5, s.n_cols);
    
    int v0n[3] = { 0, 1, 2 };
    int v1n[3] = { 1, 0, 0 };
    int v2n[3] = { 2, 2, 1 };

    //cos(sun zenith)
	double mu0 = cos(PI * theta / 180.);

    double c0, c1, c2, a0, a1, a2, bb0, bb1, bb2, b0, b1, b2, aWAT, bbWAT, bWAT;
    double tmp;
    
    mat aw = m.row(0);
    mat bbw = m.row(4);
    mat bw = bbw / B_WAT;
    
    mat a = m.rows(1, 3);
    mat bb = m.rows(5, 7);

    // specific scattering (for CHL, TSM, DOC)
    mat b = zeros<mat>(bb.n_rows, bb.n_cols);
    b.row(0) = bb.row(0) / B_CHL;
    b.row(1) = bb.row(1) / B_TSM;
    b.row(2) = bb.row(2) / B_DOC;

    for (int bn = 0; bn < s.n_cols; bn ++){

        aWAT = aw(0, bn);
        bbWAT = bbw(0, bn);
        bWAT = bw(0, bn);

        for (int vn = 0; vn < 3; vn ++){

            c0 = c(v0n[vn]);
            c1 = c(v1n[vn]);
            c2 = c(v2n[vn]);
        
            a0 = a(v0n[vn], bn);
            a1 = a(v1n[vn], bn);
            a2 = a(v2n[vn], bn);
        
            bb0 = bb(v0n[vn], bn);
            bb1 = bb(v1n[vn], bn);
            bb2 = bb(v2n[vn], bn);
    
            b0 = b(v0n[vn], bn);
            b1 = b(v1n[vn], bn);
            b2 = b(v2n[vn], bn);
        
            j(vn, bn) = j_shallow_c0(s(0, bn),
                c0, c1, c2,
                a0, a1, a2,
                bb0, bb1, bb2,
                b0, b1, b2,
                aWAT, bbWAT, bWAT,
                al1, al2, lambda(0, bn));

        }

        j(3, bn) = j_shallow_al1(s(0, bn),
                c0, c1, c2,
                a0, a1, a2,
                bb0, bb1, bb2,
                b0, b1, b2,
                aWAT, bbWAT, bWAT,
                al1, al2, lambda(0, bn));

        j(4, bn) = j_shallow_al2(s(0, bn),
                c0, c1, c2,
                a0, a1, a2,
                bb0, bb1, bb2,
                b0, b1, b2,
                aWAT, bbWAT, bWAT,
                al1, al2, lambda(0, bn));
    }

    //cout << j << endl;    
    return j;
}


mat startingCPA(double parameters[9], mat m){
    int ci0, ci1, ci2;
    double c0, c1, c2, dc0, dc1, dc2;
    double starts = parameters[1];
    int k, fullSize = starts * starts * starts;
    mat startC = zeros<mat>(fullSize, 3);
    
    double min0 = parameters[3], max0 = parameters[4];
    double min1 = parameters[5], max1 = parameters[6];
    double min2 = parameters[7], max2 = parameters[8];
    
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
                startC(k, 0) = c0;
                startC(k, 1) = c1;
                startC(k, 2) = c2;
                k ++;
                c2 += dc2;
            }
            c1 += dc1;
        }
        c0 += dc0;
    }
    return startC;
}

mat startingCPA_al(double parameters[9], mat m){
    int ci0, ci1, ci2, ci3, ci4;
    double c0, c1, c2, c3, c4, dc0, dc1, dc2, dc3, dc4;
    double starts = parameters[1];
    int k, fullSize = starts * starts * starts * starts * starts;
    mat startC = zeros<mat>(fullSize, 5);
    
    double min0 = parameters[3], max0 = parameters[4];
    double min1 = parameters[5], max1 = parameters[6];
    double min2 = parameters[7], max2 = parameters[8];
    double min3 = 0, max3 = 1.;
    double min4 = 0, max4 = 1.;
    
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

                        startC(k, 0) = c0;
                        startC(k, 1) = c1;
                        startC(k, 2) = c2;
                        startC(k, 3) = c3;
                        startC(k, 4) = c4;
                        
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

    return startC;
}
//iterface to python
//calculate Rrsw from given concentrations, albedo, depth, sola zenith
extern int get_rrsw(double parameters[9],
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *albedo, int albedo_n0,
         double h,
         double theta,
         double *outR, int outR_n0){

    int i, j, k, bands = outR_n0;
    //matrix for the model
    mat m = zeros<mat>(8, bands);
    //matrix for the concentrations
    mat c = zeros<mat>(1, 3);
    //matrix for the bottom albedo vecotr
    mat al = zeros<mat>(1, bands);
    //reconstructed reflectance
    mat r;

    //get model from Python
    //m.rows(0, 3) - specific absorption
    //m.rows(4, 7) - specific backscattering
    k = 0;
    for (i = 0; i < 8; i ++){
        for (j = 0; j < bands; j ++){
            m(i, j) = model[k];
            k ++;
        };
    };

    //get concentrations from Python
    c(0,0) = inC[0];
    c(0,1) = inC[1];
    c(0,2) = inC[2];

    //get albedo from Python
    for (i = 0; i < bands; i ++){
        al(0, i) = albedo[i];
    };

    //init HO-object
    Hydrooptics ho(m, al, h, theta);
    
    //use HO-object for estimation of Rrsw
    r = ho.rrsw(c);
    
    //sent Rrsw values back to Python
    for (i = 0; i < bands; i ++){
        outR[i] = r(0, i);
    };
    return 0;
};


//iterface to python
//calculate Rrsw from given concentrations, albedo, depth, sola zenith
extern int get_rrsw_al(double parameters[9],
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *lambda, int lambda_n0,
         double h,
         double theta,
         double al1,
         double al2,
         double *outR, int outR_n0){

    int i, j, k, bands = outR_n0;
    //matrix for the model
    mat m = zeros<mat>(8, bands);
    //matrix for the concentrations
    mat c = zeros<mat>(1, 3);
    //matrix for the wavelength vector
    mat lam = zeros<mat>(1, bands);
    //reconstructed reflectance
    mat r;

    //get model from Python
    //m.rows(0, 3) - specific absorption
    //m.rows(4, 7) - specific backscattering
    k = 0;
    for (i = 0; i < 8; i ++){
        for (j = 0; j < bands; j ++){
            m(i, j) = model[k];
            k ++;
        };
    };

    //get concentrations from Python
    c(0,0) = inC[0];
    c(0,1) = inC[1];
    c(0,2) = inC[2];

    //get wavelength from Python
    for (i = 0; i < bands; i ++){
        lam(0, i) = lambda[i];
    };

    //init HO-object
    Hydrooptics ho(m, h, theta, lam);
    
    //use HO-object for estimation of Rrsw
    r = ho.rrsw(c, al1, al2);
    
    //sent Rrsw values back to Python
    for (i = 0; i < bands; i ++){
        outR[i] = r(0, i);
    };
    return 0;
};

extern int get_c(double parameters[9],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *albedo, int albedo_rows, int albedo_cols,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length){

    int i, j, k, pi, bands = inR_cols, pixels = inR_rows;
    printf("Retrieval from %d bands x %d pixels...\n", bands, pixels);
    
    //result of optimization
    double xBest[4];
    //model
    mat m = zeros<mat>(8, bands);
    //concentrations
    mat c = zeros<mat>(1, 3);
    //bottom albedo
    mat al = zeros<mat>(1, bands);
    //measured reflectance
    mat s = zeros<mat>(1, bands);

    //get model from Python
    //m.rows(0, 3) - specific absorption
    //m.rows(4, 7) - specific backscattering
    k = 0;
    for (i = 0; i < 8; i ++){
        for (j = 0; j < bands; j ++){
            m(i, j) = model[k];
            k ++;
        };
    };

    //create starting CPA
    mat startC = startingCPA(parameters, m);

    //init HO-object
    Hydrooptics ho(m, al, h[i], theta[i]);

    //prepare for optimization with CMINPACK
    real x[3], fvec[10], fjac[30], tol, wa[300], fnorm;
    int info, ipvt[3], lwa = 100;
    //set tolerance to square of the machine recision
    tol = sqrt(__cminpack_func__(dpmpar)(1));

    //for all pixels
    for (i = 0; i < pixels; i ++){
        //get measured Rrsw from Python
        for (j = 0; j < bands; j ++)
            s(0, j) = inR[j + i*bands];
        
        //cout << s << endl;

        //get albedo from Python
        for (j = 0; j < bands; j ++)
            al(0, j) = albedo[j + i*bands];

        //set measured Rrsw, albedo, depth and solar zenith in the HO-object
        ho.s = s;
        ho.al = al;
        ho.h = h[i];
        ho.theta = theta[i];

        //erase xBest
        for (j = 0; j < 4; j ++)
            xBest[j] = 100;

        //start optimization several times and select the best result
        for (k = 0; k < startC.n_rows && xBest[3] > parameters[2]; k ++){
            // set initial concentrations
            //printf("start [%d]: ", k);
            for (j = 0; j < 3; j ++){
                x[j] = startC(k, j);
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
        //printf("iterations: %d\n", k);
        
        //sent values of concentrtions and residuals back to Python
        //printf("final result: %d ", i);
        for (j = 0; j < 4; j ++){
            //printf("%g ", xBest[j]);
            outC[j + i*4] = xBest[j];
        }
        if (fmod(i, 100.) == 0.)
            printf("%d\n", i);
    }

    printf("OK!\n");
    return 0;
};


extern int get_c_al(double parameters[9],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *lambda, int lambda_cols,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length){

    int i, j, k, pi, bands = inR_cols, pixels = inR_rows;
    printf("Retrieval from %d bands x %d pixels...\n", bands, pixels);
    
    double al10 = .5, al20 = .5;
    //result of optimization
    double xBest[6];
    //model
    mat m = zeros<mat>(8, bands);
    //concentrations
    mat c = zeros<mat>(1, 3);
    //measured reflectance
    mat s = zeros<mat>(1, bands);
    //wavelengths
    mat lam = zeros<mat>(1, bands);

    //get model from Python
    //m.rows(0, 3) - specific absorption
    //m.rows(4, 7) - specific backscattering
    k = 0;
    for (i = 0; i < 8; i ++){
        for (j = 0; j < bands; j ++){
            m(i, j) = model[k];
            k ++;
        };
    };

    //get wavelenth from Python
    k = 0;
    for (j = 0; j < bands; j ++){
        lam(0, j) = lambda[k];
        k ++;
    };
    
    //cout << lam << endl;

    //create starting CPA
    mat startC = startingCPA_al(parameters, m);

    //init HO-object
    Hydrooptics ho(m, h[i], theta[i], lam);

    //prepare for optimization with CMINPACK
    real x[5], fvec[10], fjac[50], tol, wa[500], fnorm;
    int info, ipvt[5], lwa = 100;
    //set tolerance to square of the machine recision
    tol = sqrt(__cminpack_func__(dpmpar)(1));
    tol = 1e-10;
    printf("tol: %g\n", tol);
    
    //for all pixels
    for (i = 0; i < pixels; i ++){
        //get measured Rrsw from Python
        for (j = 0; j < bands; j ++)
            s(0, j) = inR[j + i*bands];
        
        //cout << s << endl;

        //set measured Rrsw, depth and solar zenith in the HO-object
        ho.s = s;
        ho.h = h[i];
        ho.theta = theta[i];

        //erase xBest
        for (j = 0; j < 6; j ++)
            xBest[j] = 100;

        //start optimization several times and select the best result
        for (k = 0; k < startC.n_rows && xBest[5] > parameters[2]; k ++){
            // set initial concentrations
            //printf("start [%d]: ", k);
            for (j = 0; j < 5; j ++){
                x[j] = startC(k, j);
            };
            //x[3] = al10;
            //x[4] = al20;

            //for (j = 0; j < 5; j ++){
            //    printf("%5.2g ", (double)x[j]);
            //};

            //perform optimization
            info = __cminpack_func__(lmder1)(fcn_al, &ho, bands, 5, x, fvec, fjac, bands, tol, ipvt, wa, lwa);
            //estimate norm of residuals
            fnorm = __cminpack_func__(enorm)(bands, fvec);
    
            //printf(" ==> ");
            //for (j = 0; j < 5; j ++)
            //    printf("%5.2g ", (double)x[j]);
            //printf("%7.4g\n", (double)fnorm);
            
            if (fnorm < xBest[5]){
                for (j = 0; j < 5; j ++)
                    xBest[j] = (double)x[j];
                xBest[5] = (double)fnorm;
            }
        }
        printf("iterations: %d\n", k);
        
        //sent values of concentrtions and residuals back to Python
        printf("final result: %d ", i);
        for (j = 0; j < 6; j ++){
            printf("%7.4g ", xBest[j]);
            outC[j + i*6] = xBest[j];
        }
        if (fmod(i, 100.) == 0.)
            printf("%d\n", i);
    }

    printf("OK!\n");
    return 0;
};

int fcn(void *p, int m, int n, const real *x, real *fvec, real *fjac, 
	 int ldfjac, int iflag){

    int i0, i1;
    const Hydrooptics *ho = (Hydrooptics *)p;
    
    mat c;
    mat rs;
    mat j;
    
    c << x[0] << x[1] << x[2];
    //printf("fcn: c:\n");
    //cout << c << endl;
    
    
    if (iflag == 1){
        rs = ho -> rs(c);
        
        //printf("fcn: rs:\n");
        //cout << rs << endl;

        for (i0 = 0; i0 < m; i0 ++){
            fvec[i0] = rs(0, i0);
            //printf("fcn: fvec:%d %f\n", i0, fvec[i0]);
        };
    } else {

        j = ho->jacobian(c);

        //printf("fcn: j:\n");
        //cout << j << endl;
        
        for (i0 = 0; i0 < m; i0 ++){
            fjac[i0 + ldfjac*0] = j(0, i0);
            fjac[i0 + ldfjac*1] = j(1, i0);
            fjac[i0 + ldfjac*2] = j(2, i0);
            
            //printf("fcn: fjac:%d %f %f %f\n", i0, fjac[i0 + ldfjac*0], fjac[i0 + ldfjac*1], fjac[i0 + ldfjac*2]);
        };
    };

    return 0;
};

int fcn_al(void *p, int m, int n, const real *x, real *fvec, real *fjac, 
	 int ldfjac, int iflag){

    int i0, i1;
    const Hydrooptics *ho = (Hydrooptics *)p;
    double al1, al2;
    
    mat c;
    mat rs;
    mat j;
    
    c << x[0] << x[1] << x[2];
    al1 = x[3];
    al2 = x[4];
    //printf("fcn: c:\n");
    //cout << c << endl;
    
    
    if (iflag == 1){
        rs = ho -> rs(c, al1, al2);
        
        //printf("fcn: rs:\n");
        //cout << rs << endl;

        for (i0 = 0; i0 < m; i0 ++){
            fvec[i0] = rs(0, i0);
            //printf("fcn: fvec:%d %f\n", i0, fvec[i0]);
        };
    } else {

        j = ho->jacobian(c, al1, al2);

        //printf("fcn: j:\n");
        //cout << j << endl;
        
        for (i0 = 0; i0 < m; i0 ++){
            for (i1 = 0; i1 < 5; i1 ++)
                fjac[i0 + ldfjac*i1] = j(i1, i0);
                                    
            //printf("fcn: fjac:%d %f %f %f\n", i0, fjac[i0 + ldfjac*0], fjac[i0 + ldfjac*1], fjac[i0 + ldfjac*2]);
        };
    };

    return 0;
};
