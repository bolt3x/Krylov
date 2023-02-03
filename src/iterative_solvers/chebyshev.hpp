
//*****************************************************************
// Iterative template routine -- CHEBYSHEV
//
// CHEBYSHEV solves the symmetric positive definite linear
// system Ax = b using the Preconditioned Chebyshev Method
//
// CHEBYSHEV follows the algorithm described on p. 30 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************
#ifndef HH_CHEBYSHEV__HH
#define HH_CHEBYSHEV__HH

namespace Krylov
{
template <class Matrix, class Vector, class Preconditioner, class Type>
int
Chebyshev(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
      int &max_iter, typename Vector::Scalar &tol, Type const &a_min,
      Type const &a_max)
{
  using Real = typename Vector::Scalar;


  Real   resid;
  Type   alpha, beta, c, d;
  Vector p, q, z;

  Real   normb = b.norm();
  Vector r = b - A * x;

  if(normb == 0.0)
    normb = 1;

  if((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  c = (a_max - a_min) / 2.0;
  d = (a_max + a_min) / 2.0;

  for(int i = 1; i <= max_iter; i++)
    {
      z = M.solve(r); 

      if(i == 1)
        {
          p = z;
          alpha = 2.0 / d;
        }
      else
        {
          beta = c * alpha / 2.0;
          beta = beta * beta;
          alpha = 1.0 / (d - beta); 
          p = z + beta * p;        
        }

      q = A * p;
      x += alpha * p; 
      r -= alpha * q; 

      if((resid = r.norm()/ normb) <= tol)
        {
          tol = resid;
          max_iter = i;
          return 0; 
        }
    }

  tol = resid;
  return 1; // no convergence
}
} // namespace Krylov
#endif