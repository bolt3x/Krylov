#include <iostream>

#include  "./src/base_classes/matrix.hpp"
#include  "./src/base_classes/vector.hpp"
#include "./src/iterative_solvers/cg.hpp"
#include "./src/preconditioners/identity.hpp"
#include "./src/iterative_solvers/bcgstab.hpp"
#include "./src/base_classes/sparse_matrix.hpp"
#include "./src/preconditioners/diag.hpp"
#include "./src/preconditioners/spai.hpp"
#include "./src/direct_solvers/qr_solver.hpp"
#include "./src/eigen_solvers/power.hpp"
#include "./src/utils/mm_reader.hpp"
#include "./src/utils/plotter.hpp"

int main(){
 Krylov::SparseMatrix<double> A = Krylov::readCSRMatrix<double>("./matrixes/nos4.mtx");
 double n = A.cols(); 
 Krylov::Vector<double> xe(n,1);
 Krylov::Vector<double> x(n);
 
 Krylov::Vector<double> b = A * xe;
 double tol = 0.1;
 Krylov::SpaiPreconditioner<double,Krylov::PATTERN::DYNAMIC> spaiD(A.transpose(),tol);
 tol = 0.3;
 Krylov::SpaiPreconditioner<double> spaiS(A.transpose(),tol);
 Krylov::IdentityPreconditioner<double> id;
 Krylov::DiagPreconditioner<double> diag(A);
 int res;

 int max_iter;
 
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,id,max_iter,tol);
 std::cout << "Identity Preconditioner: " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,diag,max_iter,tol);
 std::cout << "Diagonal Preconditioner: " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,spaiD,max_iter,tol);
 std::cout << "Sparse Approximate Preconditioner (DYNAMIC): " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 x = Krylov::Vector<double>(n,0);
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,spaiS,max_iter,tol);
 std::cout << "Sparse Approximate Preconditioner (STATIC): " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 
 Krylov::Vector<double> y(n,1);
 double a;
 tol = 1e-9;
 max_iter = 1000;
 
 res = Power(A,a,y,max_iter,tol);
 std::cout << "Power Method: " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Distance: " << tol << std::endl;
 std::cout << "Eigenvalue: " << a << std::endl;
 std::cout << "Error norm: " << (A*y - a*y).norm() << std::endl;

 res = Krylov::plotCSRMatrix<double>(spaiD.getM(),"dynamic");
 res = Krylov::plotCSRMatrix<double>(spaiS.getM(),"static");
 return 0;
}