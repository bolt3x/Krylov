#include <iostream>

#include  "./src/base_classes/matrix.hpp"
#include  "./src/base_classes/vector.hpp"
#include "./src/iterative_solvers/cg.hpp"
#include "./src/preconditioners/identity.hpp"
#include "./src/iterative_solvers/bcgstab.hpp"
#include "./src/base_classes/sparse_matrix.hpp"
#include "./src/preconditioners/diag.hpp"
#include "./src/direct_solvers/qr_solver.hpp"
int main(){
 double n = 3; 

 Krylov::SparseMatrix<double> A(n,n);
 Krylov::Matrix<double> B(n,2);
 Krylov::Matrix<double> I(n,n);
 Krylov::Vector<double> xe(n);
 Krylov::Vector<double> x(2);

 B(0,0) = 2;
 B(0,1) = 3;
 B(1,0) = 5;
 B(1,1) = 30;
 B(2,0) = 10;
 B(2,1) = 1;

 Krylov::Vector<double> b = B * xe;


 std::cout << B*I << std::endl;
 Krylov::QRSolver<double> qr(B);
 
 std::cout << qr.getR() << std::endl;
 std::cout << qr.solve(b) << std::endl;

 return 0;
}