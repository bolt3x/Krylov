#include <iostream>

#include "matrix.hpp"
#include "cg.hpp"

int main(){

  double n = 10;
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

	double tol = 1e-2;
	int max_iter = 10;
	int res = Krylov::CG<Krylov::Matrix<double>,Krylov::Vector<double>>(A,x_e,b,max_iter,tol);

	std::cout << max_iter << std::endl;
}