/*
 * SPAI.hpp
 *
 *  Created on: Jan 18, 2023
 *      Author: bolt3x
 */

#ifndef HH_SPAI_HH
#define HH_SPAI_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <cmath>
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{

/*!
 * An identity "preconditioner" class
 * @tparam Scalar element
 * @tparam Vector The vector class to "solve" with
 */

template<typename SCALAR> class SpaiPreconditioner
{
public:
	using Scalar=SCALAR;

	/*! Constructor takes a sparse matrix
	 *  @tparam Matrix
	 *  @param mat matrix from which we compute the inverse
	 */
	template<class Matrix>
	SpaiPreconditioner(Matrix const &A) 
	{
		compute(A);
	}

	template<class Matrix>
	void compute(Matrix const &A)
	{
		std::cout << A << std::endl;
	}
	/*!
	 * solve method to apply the preconditioner 
	 * @param v 
	 * @return the preconditioned vector
	 */
	template<class Vector>
	Vector solve(Vector const &v) const {
		return v;
	}
protected:
	std
};

}//namespace Krylov
#endif