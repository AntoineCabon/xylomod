#include <Rcpp.h>
using namespace Rcpp;

//////
////// C++ implementation of cell enlargement and division model functions
//////

/// Constants
double Rn = 8.314; // The perfect gas constant
double T0 = -273.15; // Absolute 0 temperature in degC
double Tref = 15; // Reference temperature in degC

////// Effect of temperature (on metabolic rate and microtubule stability)
double temp_microT(double x, double inflection, double scale=5){
  double out = 1/(1+exp((-x+inflection)*scale));
  return out;
}

double temp_metR(double Tc, double DHa, double DSd, double DHd){
  double Tk = Tc-T0;
  double out = Tk*exp(-DHa/(Rn*Tk)) / (1+exp(DSd/Rn*(1-(DHd/(DSd*Tk)))));
  return out;
}

// [[Rcpp::export("temp_fun")]]
double temp_fun(double Tc, double Y_T=5, double DHa=87.5e3, double DSd=1.09e3, double DHd=333e3){
  double out = temp_metR(Tc, DHa, DSd, DHd);
  out = out/temp_metR(Tref, DHa, DSd, DHd); // the output is equal to 1 at Tref degC
  out = out*temp_microT(Tc, Y_T);
  // out = 1;
  return out;
}


////// Cell expansion model

// [[Rcpp::export("r_fun")]]
double r_fun(double psi, double Tc, double pi, double phi, double Y_P, double Y_T){
  double out = phi*(psi-pi-Y_P);
  if(out<0) out=0;
  out = out*temp_fun(Tc,Y_T);
  return out;
}

double pi2n(double pi, double CRD, double Tc){
  double n = -pi*CRD/(Rn*(Tc-T0));
  return n;
}

double n2pi(double n, double CRD, double Tc){
  double pi = -n*Rn*(Tc-T0)/CRD;
  return pi;
}

// [[Rcpp::export("expansion")]]
DataFrame expansion(NumericVector psi, NumericVector Tc, int start = 1, double phi0=0.13, double pi0=-0.8,
                    double CRD0=8.3, double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){
  double _Y_T = Y_T;
  int l_psi = psi.size();
  int l_Tc = Tc.size();
  int l;
  if(l_psi>l_Tc) {
    l=l_psi;
    double Tc0=Tc[0];
    NumericVector Tc_(l,Tc0);
    Tc = Tc_;
  } else if(l_psi<l_Tc) {
    l=l_Tc;
    double psi0=psi[0];
    NumericVector psi_(l,psi0);
    psi=psi_;
  } else {
    l=l_psi;
  }

  NumericVector CRD(l); CRD[0] = CRD0;
  NumericVector phi(l); phi[0] = phi0;
  NumericVector n(l); n[0] = pi2n(pi0, CRD0, Tref);
  NumericVector pi(l); pi[0] = n2pi(n[0], CRD0, Tc[0]);
  NumericVector r(l, NA_REAL);

  for(int i=0; i<l-1; i++){

    if(i<start-1){
      r[i] = 0;

      CRD[i+1] = CRD[i];
      n[i+1] = n[i];
      pi[i+1] = pi[i];
      phi[i+1] = phi[i];
    } else {
      r[i] = r_fun(psi[i], Tc[i], pi[i], phi[i], Y_P, _Y_T);

      // Variable update
      CRD[i+1] = CRD[i]*(1+r[i]); // cell volume increment
      n[i+1] = n[i]; //*(1+f); // influx of solutes
      pi[i+1] = n2pi(n[i+1], CRD[i+1], Tc[i+1]);
      phi[i+1] = phi[i] + phi[i]*(s*r[i]-h*temp_fun(Tc[i], Y_T=-999)); //changes in cell wall properties. Hardening (thickening and lignification) is temperature sensitive but not threshold prone because lignification does not need microtubules
      if(phi[i+1]<0) {
        phi[i+1]=0;
      } else {}
    }
  }

  // return outputs
  return(DataFrame::create(Named("start")=start,
                           Named("CRD")=CRD,
                           Named("r")=r,
                           Named("phi")=phi,
                           Named("pi")=pi,
                           Named("n")=n));
}

// [[Rcpp::export]]
List expansion_seq(NumericVector psi, NumericVector Tc, IntegerVector start_vec,
                   double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                   double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){
  int l_psi = psi.size();
  int l_Tc = Tc.size();
  if(l_psi>l_Tc) {
    Tc = NumericVector (l_psi,Tc[0]);
  } else if(l_psi<l_Tc) {
    psi = NumericVector (l_Tc,psi[0]);
  } else {}

  int n = start_vec.size();
  List g(n);
  for(int i=0; i<n; i++){
    g[i] = expansion(psi=psi, Tc=Tc, start_vec[i], phi0=phi0, pi0=pi0, CRD0=CRD0, Y_P=Y_P, Y_T=Y_T, h=h, s=s);
  }
  return g;
}


////// Cell division model

// [[Rcpp::export]]
DataFrame division(NumericVector psi, NumericVector Tc, double Nc = 8.85,
                       double phi=0.13, double pi=-0.8, double Y_P=0.05, double Y_T=5){
  int l_psi = psi.size();
  int l_Tc = Tc.size();
  int l;
  if(l_psi>l_Tc) {
    Tc = NumericVector (l_psi,Tc[0]);
    l = l_psi;
  } else if(l_psi<l_Tc) {
    psi = NumericVector (l_Tc,psi[0]);
    l = l_Tc;
  } else {
    l = l_Tc;
  }

  NumericVector r(l, NA_REAL);
  NumericVector P(l, NA_REAL);
  for(int i = 0; i<l; i++){
    double pi_Tcorr = n2pi(pi2n(pi,1,Tref),1,Tc[i]);
    r[i] = r_fun(psi[i], Tc[i], pi_Tcorr, phi, Y_P, Y_T);
    P[i] = r[i]/log(2)*Nc;
  }
  return(DataFrame::create(Named("r") = r,
                           Named("P") = P));
}

