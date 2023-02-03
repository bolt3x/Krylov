#include <iostream>

#include "Krylov.hpp"

int main(){
 std::cout << "Richardson" << std::endl;
 Krylov::SparseMatrix<double> A = Krylov::readCSRMatrix<double>("../matrixes/bcsstm22.mtx");
 
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
 double alpha_max;
 
 max_iter = 1000;
 tol = 1e-9;
 Krylov::Vector<double> v_max1(n,1);
 res = Power(A,alpha_max,v_max1,max_iter,tol);

 double alpha = std::abs(2 / (alpha_max * 2));
 tol = 1e-9;
 max_iter = 2000;
 res = Richardson(A,x,b,id,alpha,max_iter,tol);

 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Identity Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 max_iter = 1000;
 tol = 1e-9;
 Krylov::Vector<double> v_max2(n,1);
 Krylov::SparseMatrix<double> P1 = A * spaiD.getM();
 res = Power(A,alpha_max,v_max2,max_iter,tol);

 x = Krylov::Vector<double>(n,0);
 alpha = std::abs(2 / (alpha_max * 2));
 tol = 1e-9;
 max_iter = 5000;
 res = Richardson(A,x,b,spaiD,alpha,max_iter,tol);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Sparse Approximate Preconditioner (DYNAMIC)" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 max_iter = 5000;
 tol = 1e-9;
 Krylov::Vector<double> v_max3(n,1);
 Krylov::SparseMatrix<double> P2 = A * spaiS.getM();
 res = Power(P2,alpha_max,v_max3,max_iter,tol);
 
 x = Krylov::Vector<double>(n,0);
 alpha = std::abs(2 / (alpha_max * 2));
 tol = 1e-9;
 max_iter = 5000;
 res = Richardson(A,x,b,spaiS,alpha,max_iter,tol);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Sparse Approximate Preconditioner (STATIC)" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "-------------------------" << std::endl;

 res = Krylov::plotCSRMatrix<double>(A,"original");
 res = Krylov::plotCSRMatrix<double>(spaiD.getM(),"dynamic");
 res = Krylov::plotCSRMatrix<double>(spaiS.getM(),"static");
 return 0;
}
