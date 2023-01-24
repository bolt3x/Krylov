/*
 * Diag.hpp
 *
 *  Created on: Jan 18, 2023
 *      Author: bolt3x
 */

#ifndef HH_DIAG_HH
#define HH_DIAG_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <cmath>
#include "sparse_matrix.hpp"
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{

/*!
 * A simple diagonal preconditioner class
 * @tparam Scalar element
 * @tparam Vector The vector class to "solve" with
 */

template<typename SCALAR,class Matrix = Krylov::SparseMatrix<SCALAR>,class Vector = Krylov::Vector<SCALAR>> class DiagPreconditioner
{
public:
	using Scalar=SCALAR;

	/*! Constructor takes a sparse matrix
	 *  @param A matrix from which we compute the inverse
	 */
	DiagPreconditioner(Matrix const &A) : diag(A.rows())
	{
		compute(A);
	}

	void compute(Matrix const &A)
	{
#pragma omp parallel for
		for(std::size_t i = 0; i < A.rows(); i++)
        {
			if(A(i,i))
				diag[i] = 1 / A(i,i);  
        }	
    }
	/*!
	 * solve method to apply the preconditioner 
	 * @param v 
	 * @return the preconditioned vector
	 */
	Vector solve(Vector const &v) const {
        Vector result(v.size());
#pragma omp parallel for
        for(std::size_t i = 0; i < v.size(); i++)
				{
					result[i] = diag[i] * v[i];
				}
		return result;
	}

protected:
	std::vector<Scalar> diag;
};

}//namespace Krylov
#endif