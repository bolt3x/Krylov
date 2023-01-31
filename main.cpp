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
int main(){
 double n = 100; 

 Krylov::SparseMatrix<double> A(n,n);
 Krylov::Matrix<double> B(n,n);
 Krylov::Matrix<double> I(n,n);
 Krylov::Vector<double> xe(n);
 Krylov::Vector<double> x(n);

 for (int i=0; i<n; i++) {
    A.set(i, i) = 2.0*(i+1);
    if(i>0) A.set(i, i-1) = -i;
    if(i<n-1) A.set(i, i+1) = -(i+1);
    xe[i] = 1;
}

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

 x = 0 * x;
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,diag,max_iter,tol);
 std::cout << "Diagonal Preconditioner: " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 x = 0 * x;
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,spaiD,max_iter,tol);
 std::cout << "Sparse Approximate Preconditioner (DYNAMIC): " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;

 x = 0 * x;
 tol = 1e-9;
 max_iter = 1000;
 res = BiCGStab(A,x,b,spaiS,max_iter,tol);
 std::cout << "Sparse Approximate Preconditioner (STATIC): " << std::endl;
 std::cout << "Iterations: " << max_iter << std::endl;
 std::cout << "Residual: " << tol << std::endl;
 std::cout << "Error norm: " << (x-xe).norm() << std::endl;
 return 0;
}