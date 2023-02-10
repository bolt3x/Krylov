#include <iostream>
#include <chrono>
#include "Krylov.hpp"
using namespace std::chrono;
int main(int argc,char**argv){
 std::cout << "CG" << std::endl;
 Krylov::SparseMatrix<double> A = Krylov::readCSRMatrix<double>(argv[1]);
 double n = A.cols(); 
 Krylov::Vector<double> xe(n,1);
 Krylov::Vector<double> x(n);	
 
 Krylov::Vector<double> b = A * xe;
 //Krylov::SpaiPreconditioner<double,Krylov::PATTERN::DYNAMIC> spaiD(A);
 Krylov::SpaiPreconditioner<double> spaiS(A);
 Krylov::IdentityPreconditioner<double> id(A);
 Krylov::DiagPreconditioner<double> diag(A);
 Krylov::SymmetricSGSPreconditioner<double,Krylov::SparseMatrix<double>> sgs(A);

 int res;
 int max_iter;
 double tol;
 
 tol = 1e-9;
 max_iter = 1000;
 auto start = high_resolution_clock::now();
 res = CG(A,x,b,id,max_iter,tol);
 auto stop = high_resolution_clock::now();
 auto duration = duration_cast<microseconds>(stop - start);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "No Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "Total Time taken: " << duration.count() << std::endl;
 std::cout << "-------------------------" << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 start = high_resolution_clock::now();
 res = CG(A,x,b,diag,max_iter,tol);
 stop = high_resolution_clock::now();
 duration = duration_cast<microseconds>(stop - start);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Diagonal Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "Total Time taken: " << duration.count() << std::endl;
 std::cout << "-------------------------" << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 start = high_resolution_clock::now();
 res = CG(A,x,b,sgs,max_iter,tol);
 stop = high_resolution_clock::now();
 duration = duration_cast<microseconds>(stop - start);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Symmetric Gauss Seidel Preconditioner" << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "Total Time taken: " << duration.count() << std::endl;
 std::cout << "-------------------------" << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 2000;
 start = high_resolution_clock::now();
 res = CG(A,x,b,spaiS,max_iter,tol);
 stop = high_resolution_clock::now();
 duration = duration_cast<microseconds>(stop - start);
 std::cout << std::endl << "-------------------------" << std::endl;
 std::cout << "Sparse Approximate Preconditioner (STATIC)" << std::endl;
 std::cout << "Solver result: " << res << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 std::cout << "Total Time taken: " << duration.count() << std::endl;
 std::cout << "-------------------------" << std::endl;

 res = Krylov::plotCSRMatrix<double>(A,"original");
 ///res = Krylov::plotCSRMatrix<double>(spaiD.getM(),"dynamic");
 res = Krylov::plotCSRMatrix<double>(spaiS.getM(),"static");
 
 return 0;
}
