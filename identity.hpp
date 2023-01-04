/*
 * Preconditioner.hpp
 *
 *  Created on: Dec 31, 2022
 *      Author: bolt3x
 */

#ifndef HH_PRECONDITIONER_HH
#define HH_PRECONDITIONER_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <cmath>
#include "matrix.hpp"
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

template<typename SCALAR,class Vector = Krylov::Vector<SCALAR>> class IdentityPreconditioner
{
public:
	using Scalar=SCALAR;

	/*! Constructor doesn't take anything
	 */
	IdentityPreconditioner() = default;

	/*!
	 * solve method to apply the preconditioner on the residual
	 * @param v 
	 * @return the preconditioned vector
	 */
	Vector solve(Vector const &v) const {
		return v;
	}
};

}//namespace Krylov
#endif