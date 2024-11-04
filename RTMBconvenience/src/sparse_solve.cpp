// [[Rcpp::depends(TMB)]]
#include <RTMB.h>
#include "RTMBconvenience_types.hpp"

// [[Rcpp::export]]
Rcpp::ComplexMatrix sparse_solve_ad(Rcpp::S4 &x, Rcpp::ComplexMatrix &y){

  CHECK_INPUT(y);
  ConstMapMatrix Ey = MatrixInput(y);
  Eigen::SparseMatrix<ad> Ex = SparseInput(x);
  Eigen::SparseLU< Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> > solver;
  solver.compute(Ex);
  matrix<ad> Ey_ad(Ey);
  matrix<ad> r = solver.solve(Ey_ad);
  return MatrixOutput(r);
}

// [[Rcpp::export]]
Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> >> sparse_solve_ptr_ad(Rcpp::S4 &x){

  Eigen::SparseMatrix<ad> Ex = SparseInput(x);
  typedef Eigen::SparseLU< Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> > SpSolver;
  SpSolver* ss = new SpSolver();
  ss->analyzePattern(Ex);
  ss->factorize(Ex);
  return Rcpp::XPtr<SpSolver> (ss); 
}

// [[Rcpp::export]]
void sparse_solve_ptr_update_ad(Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> >> &ss, Rcpp::S4 &x){
  Eigen::SparseMatrix<ad> Ex = SparseInput(x);  
  ss->factorize(Ex);
  return;
}


// [[Rcpp::export]]
Rcpp::ComplexMatrix sparse_solve_ptr_eval_ad(Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> >> &ss, Rcpp::ComplexMatrix &y){
  CHECK_INPUT(y);
  ConstMapMatrix Ey = MatrixInput(y);
  matrix<ad> Ey_ad(Ey);
  matrix<ad> r = ss->solve(Ey_ad);
  return MatrixOutput(r);
}

ConstMapMatrix_NotAD MatrixInput_NotAD(const Rcpp::NumericMatrix &x) {
  return ConstMapMatrix_NotAD ((double*) x.begin(), x.nrow(), x.ncol());
}

Rcpp::NumericMatrix MatrixOutput_NotAD(const matrix<double> &X) {
  //MapMatrix_NotAD Z((double*) z.begin(), z.nrow(), z.ncol());
  Rcpp::NumericMatrix z(X.rows(), X.cols(), (double*)X.data());  
  return z;
}

Eigen::SparseMatrix<double> SparseInput_NotAD(Rcpp::S4 S) {
  Rcpp::NumericVector x(S.slot("x"));
  Rcpp::IntegerVector i = S.slot("i");
  Rcpp::IntegerVector p = S.slot("p");
  Rcpp::IntegerVector Dim = S.slot("Dim");
  return Eigen::Map<const Eigen::SparseMatrix<double> > (Dim[0], // rows()
							 Dim[1], // cols()
							 i.size(), // nonZeros()
							 p.begin(), // outerIndexPtr()
							 i.begin(), // innerIndexPtr()
							 (double*) x.begin(), // data()
							 NULL); // innerNonZeroPtr();
}

Rcpp::S4 SparseOutput_NotAD (const Eigen::SparseMatrix<double> &S) {
  size_t nnz  = S.nonZeros();
  Rcpp::IntegerVector Dim(2);
  Dim[0] = S.rows();
  Dim[1] = S.cols();
  Rcpp::IntegerVector i(S.innerIndexPtr(), S.innerIndexPtr() + nnz);
  Rcpp::IntegerVector p(S.outerIndexPtr(), S.outerIndexPtr() + Dim[1] + 1);
  Rcpp::NumericVector x( (S.valuePtr()),
                         (S.valuePtr() + nnz));
  Rcpp::S4 ans("dgCMatrix");
  ans.slot("x") = x;
  ans.slot("i") = i;
  ans.attr("p") = p;
  ans.attr("Dim") = Dim;
  return ans;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix sparse_solve_NotAD(Rcpp::S4 &x, Rcpp::NumericMatrix &y){
  ConstMapMatrix_NotAD Ey = MatrixInput_NotAD(y);
  Eigen::SparseMatrix<double> Ex = SparseInput_NotAD(x);
  Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
  solver.compute(Ex);
  matrix<double> Ey_ad(Ey);
  matrix<double> r = solver.solve(Ey_ad);
  return MatrixOutput_NotAD(r);
}


// [[Rcpp::export]]
Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >> sparse_solve_ptr_NotAD(Rcpp::S4 &x){

  Eigen::SparseMatrix<double> Ex = SparseInput_NotAD(x);
  typedef Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SpSolver;
  SpSolver* ss = new SpSolver();
  ss->analyzePattern(Ex);
  ss->factorize(Ex);
  return Rcpp::XPtr<SpSolver> (ss); 
}

// [[Rcpp::export]]
void sparse_solve_ptr_update_NotAD(Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >> &ss, Rcpp::S4 &x){
  Eigen::SparseMatrix<double> Ex = SparseInput_NotAD(x);  
  ss->factorize(Ex);
  return;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix sparse_solve_ptr_eval_NotAD(Rcpp::XPtr<Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >> &ss, Rcpp::NumericMatrix &y){
  ConstMapMatrix_NotAD Ey = MatrixInput_NotAD(y);
  matrix<double> Ey_ad(Ey);
  matrix<double> r = ss->solve(Ey_ad);
  return MatrixOutput_NotAD(r);
}
