/*
 * Identity.hpp
 *
 *  Created on: Jan 04, 2023
 *      Author: bolt3x
 */

#ifndef HH_ID_HH
#define HH_ID_HH
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
 */

template<typename SCALAR> class IdentityPreconditioner
{
public:
	using Scalar=SCALAR;

	/*! Constructor doesn't take anything
	 */
	template<class Matrix>
	IdentityPreconditioner(const Matrix &A){}

	/*!
	 * solve method to apply the preconditioner
	 * @tparam Vector class 
	 * @param v 
	 * @return the preconditioned vector
	 */
	template<class Vector>
	Vector solve(Vector const &v) const {
		return v;
	}
};

}//namespace Krylov
#endif