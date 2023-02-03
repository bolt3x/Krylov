#include <iostream>

#include "Krylov.hpp"

int main(){
 std::cout << "Chebyshev" << std::endl;     
 Krylov::SparseMatrix<double> A = Krylov::readCSRMatrix<double>("../matrixes/diffreact.mtx");

 double n = A.cols();
 Krylov::Vector<double> xe(n,1);
 Krylov::Vector<double> x(n);

 Krylov::Vector<double> b = A * xe;
 Krylov::SpaiPreconditioner<double,Krylov::PATTERN::DYNAMIC> spaiD(A.transpose());
 Krylov::SpaiPreconditioner<double> spaiS(A.transpose());
 Krylov::IdentityPreconditioner<double> id(A);
 Krylov::DiagPreconditioner<double> diag(A);
 int res;
 int max_iter;
 double tol;
 double a_max;
 double a_min;

 max_iter = 1000;
 tol = 1e-15;
 Krylov::Vector<double> v_max(n,1);
 res = Power(A,a_max,v_max,max_iter,tol);

 max_iter = 1000;
 tol = 1e-15;
 Krylov::Vector<double> v_min(n,1);
 res = InversePower(A,a_min,v_min,Krylov::CG,id,max_iter,tol);
 
 tol = 1e-9;
 max_iter = 1000;
 res = Chebyshev(A,x,b,id,max_iter,tol,a_min,a_max);

 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Computed extreme eigenvalues: " << a_min << ", " << a_max << std::endl;
 std::cout << "Identity Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 res = Krylov::plotCSRMatrix<double>(A,"original");
 res = Krylov::plotCSRMatrix<double>(spaiD.getM(),"dynamic");
 res = Krylov::plotCSRMatrix<double>(spaiS.getM(),"static");

 
}
