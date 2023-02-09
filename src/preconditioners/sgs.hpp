/*
 * SGS.hpp
 *
 *  Created on: Feb 08, 2023
 *      Author: bolt3x
 */

#ifndef HH_SGS_HH
#define HH_SGS_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <cmath>
#include <optional>

// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{

/*!
 * A symmetric gauss seidel preconditioner class
 * @tparam Scalar element
 */

template<typename SCALAR,class Matrix> class SymmetricSGSPreconditioner
{
public:
	using Scalar=SCALAR;

	/*! Constructor takes a sparse matrix
	 *  @param A matrix from which we compute the inverse
	 */
	SymmetricSGSPreconditioner(Matrix const &A) : diag(A.rows()), SGS(A)
	{
		compute(A);
	}

	void compute(Matrix const &A)
	{
#pragma omp parallel for
		for(std::size_t i = 0; i < A.rows(); i++)
        {
			
	    	diag[i] = A(i,i); 
        }	
    }

    template<class Vector>
	Vector backSolve(Vector const &v) const
    {   
        Vector res(v.size());

#pragma omp parallel for
		for(std::size_t i = 0; i < res.size(); i++)
		{
			res[i] = v[i];
		}

        
		for(int i = res.size() -1; i >= 0; i--)
		{
		    res[i] = res[i] / SGS(i,i);
            
#pragma omp parallel for
			for(int j = 0; j < i; j++)
			{
				res[j] -= SGS(j,i) * res[i];
			}
		}
		return res;
    }

    template<class Vector>
	Vector forwardSolve(Vector const &v) const
    {   
        Vector res(v.size());
        res[0] = v[0] / SGS(0,0);

        for(std::size_t i = 1; i < res.size(); i++)
        {
            res[i] = v[i];
            for(std::size_t j = 0; j < i; j++)
            {
                res[i] -= res[j] * SGS(i,j);
            }

            res[i] = res[i] / SGS(i,i);
        }
        
		return res;
    }

	/*!
	 * solve method to apply the preconditioner 
	 * @tparam Vector
	 * @param v 
	 * @return the preconditioned vector
	 */
	template<class Vector>
	Vector solve(Vector const &v) const {

        Vector res1 = forwardSolve(v);

#pragma omp parallel for
        for(std::size_t i = 0; i < v.size(); i++)
				{
					res1[i] = diag[i] * res1[i];
				}

        Vector res2 = backSolve(res1);
        
		return res2;
	}


protected:
	std::vector<Scalar> diag;
    const Matrix& SGS;
    
};

}//namespace Krylov
#endif