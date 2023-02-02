//*****************************************************************
// Iterative template routine -- Power Method
//
// Power method is used to compute the biggest in modulus 
// eigenvalue of A
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        a  --  approximate eigenvalue
//        x  --  approximate eigenvector corresponding to a
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the distance between the final iteration and the previous one
//
//*****************************************************************

#ifndef HH_POWER_HH
#define HH_POWER_HH

namespace Krylov
{

template<class Matrix,class Vector> 
int
Power(const Matrix &A,typename Matrix::Scalar &a, Vector &x,int &max_iter,typename Matrix::Scalar &tol)
{
    using Real = typename Matrix::Scalar;

    Real distance;

    for(std::size_t i = 0; i < max_iter; i++)
    {
        Vector x_new = A * x;
        Real u = 1/(x_new.norm());
        x_new = u * x_new;
        distance = (x_new - x).norm();
        if( distance <= tol)
        {
            x = x_new;
            a = x.dot(A*x);
            max_iter = i;
            tol = distance;
            return 0;
        }

        x = x_new;
    }

    a = x.dot(A*x) / (x.dot(x));
    tol = distance;
    return 1;
    
}

}//namespace Krylov
#endif