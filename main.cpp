#include <iostream>

#include "matrix.hpp"
#include "cg.hpp"
#include "identity.hpp"

int main(){

  double n = 1000;
  Krylov::Vector<double> x(n);
	Krylov::Vector<double> x_e(n);
  Krylov::Vector<double> b(n);
  Krylov::Matrix<double> A(n,n);

  for(int i = 0; i < n; i++){
        
    x[i] = 1;
        
    A(i, i) = 2.0*(i+1);
    if(i>0) A(i, i-1) = -i;
    if(i<n-1) A(i, i+1) = -(i+1);

  }

	b = A * x;

	double tol = 1e-10;
	int max_iter = 1000;
  Krylov::IdentityPreconditioner<double> id;
	int res = Krylov::CG(A,x_e,b,id,max_iter,tol);

	std::cout << max_iter << std::endl;
  std::cout << tol << std::endl;
  std::cout << (x - x_e).norm() << std::endl;

}