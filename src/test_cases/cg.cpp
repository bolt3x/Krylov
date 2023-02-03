#include <iostream>

#include "Krylov.hpp"

int main(){
 std::cout << "CG" << std::endl;
 Krylov::SparseMatrix<double> A = Krylov::readCSRMatrix<double>("../matrixes/nos4.mtx");
 double n = A.cols(); 
 Krylov::Vector<double> xe(n,1);
 Krylov::Vector<double> x(n);	
 
 Krylov::Vector<double> b = A * xe;
 double tol = 0.1;
 Krylov::SpaiPreconditioner<double,Krylov::PATTERN::DYNAMIC> spaiD(A,tol);
 tol = 0.3;
 Krylov::SpaiPreconditioner<double> spaiS(A,tol);
 Krylov::IdentityPreconditioner<double> id;
 Krylov::DiagPreconditioner<double> diag(A);
 int res;
 int max_iter;
 
 tol = 1e-9;
 max_iter = 1000;
 res = CG(A,x,b,id,max_iter,tol);

 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Identity Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 res = CG(A,x,b,diag,max_iter,tol);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Diagonal Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 res = Krylov::plotCSRMatrix<double>(A,"original");
 res = Krylov::plotCSRMatrix<double>(spaiD.getM(),"dynamic");
 res = Krylov::plotCSRMatrix<double>(spaiS.getM(),"static");
 return 0;
}
