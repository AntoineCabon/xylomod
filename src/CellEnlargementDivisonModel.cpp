#include <Rcpp.h>
using namespace Rcpp;

//////
////// C++ implementation of cell enlargement and division model functions
//////

/// Constants
double Rn = 8.314; // The perfect gas constant
double T0 = -273.15; // Absolute 0 temperature in degC
double Tref = 15; // Reference temperature in degC

/// Utility functions
List _length_equalizer(NumericVector V1, NumericVector V2){
  // Increase vectors length by repeating values in order to have final length equal to the longest
  int l1 = V1.size();
  int l2 = V2.size();
  if(l1>l2) {
    double V20=V2[0];
    NumericVector _V2(l1,V20);
    V2 = _V2;
  } else {
    double V10=V1[0];
    NumericVector _V1(l2,V10);
    V1=_V1;
  }
  return List::create(V1, V2);
}

NumericMatrix _matrix_enhancer(NumericMatrix mat, IntegerVector add_dim){
  /// Take a matrix and increases the number of rows and/or columns. New rows and columns are filled with zeros
  int ncols = mat.ncol();
  int nrows = mat.nrow();
  NumericMatrix temp(nrows+add_dim[0], ncols+add_dim[1]);

  // fill new matrix with values from previous matrix
  for(int i=0; i<nrows; i++){
    for(int j=0; j<ncols; j++){
      temp(i,j) = mat(i,j);
    }
  }

  mat = temp;
  return mat;
}

////// Effect of temperature (on metabolic rate and microtubule stability)
double _microT(double Tc, double inflection, double scale=5){
  double out = 1/(1+exp((-Tc+inflection)*scale));
  return out;
}

double _metR(double Tc, double DHa, double DSd, double DHd){
  double Tk = Tc-T0;
  double out = Tk*exp(-DHa/(Rn*Tk)) / (1+exp(DSd/Rn*(1-(DHd/(DSd*Tk)))));
  return out;
}

// [[Rcpp::export("T_fun")]]
double T_fun(double Tc, double Y_T=8, double DHa=87.5e3, double DSd=1.09e3, double DHd=333e3){
  double out = _metR(Tc, DHa, DSd, DHd);
  out = out/_metR(Tref, DHa, DSd, DHd); // the output is equal to 1 at Tref degC
  out = out*_microT(Tc, Y_T);
  // out = 1;
  return out;
}

//// Convert osmotic potential to osmolyte quantity and back
double _pi2n(double pi, double V, double Tc){
  double n = -pi*V/(Rn*(Tc-T0));
  return n;
}

double _n2pi(double n, double V, double Tc){
  double pi = -n*Rn*(Tc-T0)/V;
  return pi;
}


////// Cell expansion model

double _r(double psi, double Tc, double pi, double phi, double Y_P, double Y_T){
  double out = phi*(psi-pi-Y_P);
  if(out<0) out=0;
  out = out*T_fun(Tc,Y_T);
  return out;
}

NumericVector _expand(double psi, double Tc,
                      double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                      double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){
  // default parameters from Cabon et al. New Phytologist 2020. h is different because of potential error in the ref value.

  double n = _pi2n(pi0, CRD0, Tref);
  pi0 = _n2pi(n, CRD0, Tc); // updates the value of pi0 which is given at Tref for the current temperature Tc

  // Calculate relative volume expansion rate
  double r = _r(psi, Tc, pi0, phi0, Y_P, Y_T);

  // Variable update
  double CRD1 = CRD0*(1+r); // cell diameter (volume) increment
  double pi1 = _n2pi(n, CRD1, Tref); // pi is returned at Tref in order to be consistent with input
  double phi1 = phi0 + phi0*(s*r - h*T_fun(Tc, -999)); //changes in cell wall properties. Hardening (thickening and lignification) is temperature sensitive but not threshold prone because lignification does not need microtubules
  if(phi1<0) {
    phi1=0;
  } else {}

  // return outputs
  return(NumericVector::create(Named("phi")=phi1,
                               Named("pi")=pi1,
                               Named("CRD")=CRD1));
}

// [[Rcpp::export]]
DataFrame expand(NumericVector psi, NumericVector Tc,
                 double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                 double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  List temp_list = _length_equalizer(psi, Tc);
  psi = temp_list[0];
  Tc = temp_list[1];
  int l = psi.size();

  NumericVector CRD(l);
  NumericVector phi(l);
  NumericVector pi(l);
  // NumericVector r(l);
  NumericVector temp_vec;

  for(int i=0; i<l; i++){
    if(i==0){
      temp_vec = _expand(psi[i], Tc[i], phi0, pi0, CRD0, Y_P, Y_T, h, s);
    } else {
      temp_vec = _expand(psi[i], Tc[i], phi[i-1], pi[i-1], CRD[i-1], Y_P, Y_T, h, s);
    }
    // Variable update
    // r[i] = temp_vec["r"];
    pi[i] = temp_vec["pi"];
    phi[i] = temp_vec["phi"];
    CRD[i] = temp_vec["CRD"];
  }

  // return outputs
  return(DataFrame::create(Named("phi")=phi,
                           Named("pi")=pi,
                           Named("CRD")=CRD));
}

List _expand_ring(List ring, double psi, double Tc,
                 double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                 double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  if(ring[0] == R_NilValue){
    // Create new matrices
    NumericMatrix phi(1,1);
    NumericMatrix pi(1,1);
    NumericMatrix CRD(1,1);

    NumericVector temp = _expand(psi, Tc, phi0, pi0, CRD0, Y_P, Y_T, h, s);

    phi(0,0) = temp["phi"];
    pi(0,0) = temp["pi"];
    CRD(0,0) = temp["CRD"];

    ring = List::create(_["phi"]=phi,
                        _["pi"]=pi,
                        _["CRD"]=CRD);
  } else {
    // Increase number of rows by one
    NumericMatrix phi = _matrix_enhancer(ring["phi"], IntegerVector {1,0});
    NumericMatrix pi = _matrix_enhancer(ring["pi"], IntegerVector {1,0});
    NumericMatrix CRD = _matrix_enhancer(ring["CRD"], IntegerVector {1,0});
    int nrows = CRD.nrow();
    int ncols = CRD.ncol();

    // fill new row with new simulated values
    for(int i=0; i<ncols; i++){
      NumericVector temp = _expand(psi, Tc,
                                   phi(nrows-2,i),
                                   pi(nrows-2,i),
                                   CRD(nrows-2,i),
                                   Y_P, Y_T, h, s);

      phi(nrows-1,i) = temp["phi"];
      pi(nrows-1,i) = temp["pi"];
      CRD(nrows-1,i) = temp["CRD"];
    }
    ring["phi"] = phi;
    ring["pi"] = pi;
    ring["CRD"] = CRD;
  }
  return ring;
}

// [[Rcpp::export]]
List expand_ring(List ring, NumericVector psi, NumericVector Tc,
                  double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                  double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  List temp = _length_equalizer(psi, Tc);
  psi = temp[0];
  Tc = temp[1];
  int l = psi.size();

  for(int i=0; i<l; i++){
    ring = _expand_ring(ring, psi[i], Tc[i], phi0, pi0, CRD0, Y_P, Y_T, h, s);
  }

  return ring;
}


////// Cell division model

NumericVector _divide(double psi, double Tc,
                      double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
                      double Y_P=0.05, double Y_T=8){
  // Default parameters from Cabon et al New Phytologist (2020)

  double r; //  Cell relative growth rate
  double P; // Cell production rate
  double pi_Tcorr = _n2pi(_pi2n(pi0,1,Tref),1,Tc);
  r = _r(psi, Tc, pi_Tcorr, phi0, Y_P, Y_T);
  P = r/log(2)*Nc;

  return(NumericVector::create(Named("r") = r,
                               Named("P") = P));
}

// [[Rcpp::export]]
DataFrame divide(NumericVector psi, NumericVector Tc,
                 double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
                 double Y_P=0.05, double Y_T=5){
  List temp_list = _length_equalizer(psi, Tc);
  psi = temp_list[0];
  Tc = temp_list[1];
  int l = psi.size();

  NumericVector r(l);
  NumericVector P(l);
  for(int i = 0; i<l; i++){
    NumericVector temp_vec = _divide(psi[i], Tc[i], Nc, phi0, pi0, Y_P, Y_T);
    r[i] = temp_vec["r"];
    P[i] = temp_vec["P"];
  }
  return(DataFrame::create(Named("r") = r,
                           Named("P") = P));
}


//// Combine expansion and division functions to simulate ring growth

List _grow_ring(List ring, double psi, double Tc,
                double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0 = 8.3,
                double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){
  List exp = _expand_ring(ring, psi, Tc, phi0, pi0, CRD0, Y_P, Y_T, h, s);
  NumericVector div = _divide(psi, Tc, Nc, phi0, pi0, Y_P, Y_T);

  if(ring[0]==R_NilValue){
    ring["P"] = div["P"];
  }else{
    NumericVector P = ring["P"];
    P.push_back(div["P"]);
    ring["P"] = P;
  }

  ring["phi"] = exp["phi"];
  ring["pi"] = exp["pi"];
  ring["CRD"] = exp["CRD"];

  return(ring);
}

// [[Rcpp::export]]
List grow_ring(List ring, NumericVector psi, NumericVector Tc,
                double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0 = 8.3,
                double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){

  List temp_list = _length_equalizer(psi, Tc);
  psi = temp_list[0];
  Tc = temp_list[1];
  int l = psi.size();

  for(int i=0; i<l; i++){
    ring = _grow_ring(ring, psi[i], Tc[i],
                      Nc, phi0, pi0, CRD0,
                      Y_P, Y_T, h, s);
  }
  return(ring);
}
