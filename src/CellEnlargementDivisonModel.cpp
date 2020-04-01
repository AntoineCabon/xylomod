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
// List _length_equalizer(NumericVector V1, NumericVector V2){
//   // Increase vectors length by repeating values in order to have final length equal to the longest
//   int l1 = V1.size();
//   int l2 = V2.size();
//   if(l1>l2) {
//     double V20=V2[0];
//     NumericVector _V2(l1,V20);
//     V2 = _V2;
//   } else {
//     double V10=V1[0];
//     NumericVector _V1(l2,V10);
//     V1=_V1;
//   }
//   return List::create(V1, V2);
// }

// NumericMatrix _matrix_enhancer(NumericMatrix mat, IntegerVector add_dim){
//   /// Take a matrix and increases the number of rows and/or columns. New rows and columns are filled with zeros
//   int ncols = mat.ncol();
//   int nrows = mat.nrow();
//   NumericMatrix temp(nrows+add_dim[0], ncols+add_dim[1]);
//
//   // fill new matrix with values from previous matrix
//   for(int i=0; i<nrows; i++){
//     for(int j=0; j<ncols; j++){
//       temp(i,j) = mat(i,j);
//     }
//   }
//
//   mat = temp;
//   return mat;
// }

double _sum(NumericVector v){
  int l = v.size();
  double res;
  for(int i=0; i<l; i++){
    res = res+v[i];
  }
  return res;
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

DataFrame _expand(double psi, double Tc,
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
  return(DataFrame::create(_["phi"]=phi1,
                           _["pi"]=pi1,
                           _["CRD"]=CRD1));
}

// [[Rcpp::export]]
DataFrame expand(DataFrame data,
                 double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                 double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  NumericVector psi = data["psi"];
  NumericVector Tc = data["Tc"];
  StringVector date = data["date"];
  int l = data.nrow();

  NumericVector CRD(l);
  NumericVector phi(l);
  NumericVector pi(l);

  DataFrame temp_df;

  for(int i=0; i<l; i++){
    if(i==0){
      temp_df = _expand(psi[i], Tc[i], phi0, pi0, CRD0, Y_P, Y_T, h, s);
    } else {
      temp_df = _expand(psi[i], Tc[i], phi[i-1], pi[i-1], CRD[i-1], Y_P, Y_T, h, s);
    }
    // Variable update
    // r[i] = temp_vec["r"];
    pi[i] = temp_df["pi"];
    phi[i] = temp_df["phi"];
    CRD[i] = temp_df["CRD"];
  }

  // return outputs
  return(DataFrame::create(Named("date")=date,
                           Named("phi")=phi,
                           Named("pi")=pi,
                           Named("CRD")=CRD));
}

List _expand_ring(List ring, double psi, double Tc,
                       double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  int l = cells.nrow();

  for(int i=0; i<l; i++){
    DataFrame temp = _expand(psi, Tc,
                             phi[i], pi[i], CRD[i],
                             Y_P, Y_T, h, s);

    phi[i] = temp["phi"];
    pi[i] = temp["pi"];
    CRD[i] = temp["CRD"];
  }

  cells["phi"] = phi;
  cells["pi"] = pi;
  cells["CRD"] = CRD;

  ring["cells"] = cells;
  return ring;
}

// [[Rcpp::export]]
List expand_ring(List ring, DataFrame data,
                  double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){

  NumericVector psi = data["psi"];
  NumericVector Tc = data["Tc"];
  StringVector date = data["date"];
  int nrows = data.nrow();

  DataFrame cells = as<DataFrame>(ring["cells"]);
  int ncols = cells.nrow();

  if(ring.attr("historic")){

    NumericMatrix phi(nrows, ncols);
    NumericMatrix pi(nrows, ncols);
    NumericMatrix CRD(nrows, ncols);

    for(int i=0; i<nrows; i++){
      ring = _expand_ring(ring, psi[i], Tc[i], Y_P, Y_T, h, s);

      cells = as<DataFrame>(ring["cells"]);
      NumericVector phi_i = cells["phi"];
      NumericVector pi_i = cells["pi"];
      NumericVector CRD_i = cells["CRD"];

      for(int j=0; j<ncols; j++){
        phi(i,j) = phi_i[j];
        pi(i,j) = pi_i[j];
        CRD(i,j) = CRD_i[j];
      }
    }

    StringVector colnm  = as<StringVector>(cells["formation_date"]);
    for(int j=0; j<ncols; j++){
      String str = "C";
      str += j+1;
      str += "_";
      str += colnm[j];
      colnm[j] = str;
    }

    rownames(phi) = date;
    colnames(phi) = colnm;
    rownames(pi) = date;
    colnames(pi) = colnm;
    rownames(CRD) = date;
    colnames(CRD) = colnm;

    List cells_historic;
    cells_historic["phi"] = phi;
    cells_historic["pi"] = pi;
    cells_historic["CRD"] = CRD;
    ring["cells_historic"] = cells_historic;

  } else {
    for(int i=0; i<nrows; i++){

      ring = _expand_ring(ring, psi[i], Tc[i], Y_P, Y_T, h, s);

    }
  }

  return ring;
}


////// Cell division model
double _divide(double psi, double Tc,
                      double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
                      double Y_P=0.05, double Y_T=8){
  // Default parameters from Cabon et al New Phytologist (2020)

  double r; //  Cell relative growth rate
  double P; // Cell production rate
  double pi_Tcorr = _n2pi(_pi2n(pi0,1,Tref),1,Tc);
  r = _r(psi, Tc, pi_Tcorr, phi0, Y_P, Y_T);
  P = r/log(2)*Nc;

  return(P);
}

// [[Rcpp::export]]
DataFrame divide(DataFrame data,
                     double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
                     double Y_P=0.05, double Y_T=5){

  NumericVector psi = data["psi"];
  NumericVector Tc = data["Tc"];
  StringVector date = data["date"];
  int l = data.nrow();

  NumericVector P(l);
  for(int i = 0; i<l; i++){
    P[i] = _divide(psi[i], Tc[i], Nc, phi0, pi0, Y_P, Y_T);
  }

  return(DataFrame::create(_["date"] = date, _["P"] = P));
}


//// Combine expansion and division functions to simulate ring growth
List _grow_ring(List ring, double psi, double Tc, String date,
                 double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                 double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){

  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  StringVector formation_date = cells["formation_date"];

  DataFrame divisions = as<DataFrame>(ring["divisions"]);
  NumericVector P = divisions["P"];
  StringVector P_dates = divisions["date"];

  // Calculate cell production
  double P_i = _divide(psi, Tc, Nc, phi0, pi0, Y_P, Y_T);
  // Update the cell production vector
  P.push_back(P_i);
  P_dates.push_back(date);
  divisions = DataFrame::create(_["date"] = P_dates,
                                _["P"] = P);

  if (ring.attr("cell_wise")){

    // Calculate the whole number of cells formed at the current and previous timestep
    double Pold = _sum(divisions["P"]); int Pold_int = floor(Pold);
    double Pnew = Pold+P_i; int Pnew_int = floor(Pnew);
    // Add a new value in the cell expansion vectors for each whole cell number increment
    for (int i=0; i < Pnew_int-Pold_int; i++){
      // Update cell expansion vectors
      phi.push_back(phi0);
      pi.push_back(pi0);
      CRD.push_back(CRD0);
      formation_date.push_back(date);
    }

  } else {
    // Update cell expansion vectors
    phi.push_back(phi0);
    pi.push_back(pi0);
    CRD.push_back(CRD0);
    formation_date.push_back(date);
  }

  // Create new cells data frame
  cells = DataFrame::create(_["formation_date"] = formation_date,
                            _["phi"] = phi,
                            _["pi"] = pi,
                            _["CRD"] = CRD);
  // Update ring object
  ring["divisions"] = divisions;
  ring["cells"] = cells;
  // Calculate cell expansion
  ring = _expand_ring(ring, psi, Tc, Y_P, Y_T, h, s);

  return(ring);
}

// [[Rcpp::export]]
List grow_ring(List ring, DataFrame data,
               double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
               double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){

  NumericVector psi = data["psi"];
  NumericVector Tc = data["Tc"];
  StringVector date = data["date"];
  int nrows = data.nrow();
  int ncols;

  DataFrame cells = as<DataFrame>(ring["cells"]);
  int previous_cols = cells.nrow();

  DataFrame divisions = as<DataFrame>(ring["divisions"]);

  if (ring.attr("cell_wise")){

    double previous_P = _sum(divisions["P"]);
    double remaining_P = previous_P - floor(previous_P);
    DataFrame new_divisions = divide(data, Nc, phi0, pi0, Y_P, Y_T);
    double new_P = _sum(new_divisions["P"]);
    int new_cols = floor(new_P+remaining_P);
    ncols = previous_cols+new_cols;

  } else {
    ncols = previous_cols+nrows;
  }

  if(ring.attr("historic")){

    NumericMatrix phi(nrows, ncols);
    NumericMatrix pi(nrows, ncols);
    NumericMatrix CRD(nrows, ncols);

    for(int i=0; i<nrows; i++){
      ring = _grow_ring(ring, psi[i], Tc[i], date[i],
                         Nc, phi0, pi0, CRD0,
                         Y_P, Y_T, h, s);

      cells = ring["cells"];
      NumericVector phi_i = cells["phi"];
      NumericVector pi_i = cells["pi"];
      NumericVector CRD_i = cells["CRD"];
      int l = CRD_i.size();

      for(int j=0; j<l; j++){
        phi(i,j) = phi_i[j];
        pi(i,j) = pi_i[j];
        CRD(i,j) = CRD_i[j];
      }
    }

    StringVector colnm  = as<StringVector>(cells["formation_date"]);
    for(int j=0; j<ncols; j++){
      String str = "C";
      str += j+1;
      str += "_";
      str += colnm[j];
      colnm[j] = str;
    }

    rownames(phi) = date;
    colnames(phi) = colnm;
    rownames(pi) = date;
    colnames(pi) = colnm;
    rownames(CRD) = date;
    colnames(CRD) = colnm;

    List cells_historic;
    cells_historic["phi"] = phi;
    cells_historic["pi"] = pi;
    cells_historic["CRD"] = CRD;
    ring["cells_historic"] = cells_historic;

  } else {
    for(int i=0; i<nrows; i++){

      ring = _grow_ring(ring, psi[i], Tc[i], date[i],
                         Nc, phi0, pi0, CRD0,
                         Y_P, Y_T, h, s);

    }
  }
  return ring;
}


// [[Rcpp::export]]
List initialize_ring(bool cell_wise = 0, bool historic = 0){

  StringVector date;
  NumericVector P;
  NumericVector phi;
  NumericVector pi;
  NumericVector CRD;

  DataFrame cells = DataFrame::create(_["formation_date"] = date,
                                      _["phi"] = phi,
                                      _["pi"] = pi,
                                      _["CRD"] = CRD);

  DataFrame divisions = DataFrame::create(_["date"] = date,
                                          _["P"] = P);

  List ring = List::create(_["divisions"] = divisions,
                           _["cells"] = cells);

  if(historic){
    List cells_historic;
    ring["cells_historic"] = cells_historic;
  }

  ring.attr("cell_wise") = cell_wise;
  ring.attr("historic") = historic;

  return ring;
}
