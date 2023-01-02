/*
 * Preconditioner.hpp
 *
 *  Created on: Dec 31, 2022
 *      Author: bolt3x
 */

#ifndef HH_VECTOR_HH
#define HH_VECTOR_HH
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
 * A preconditioner abstract class
 * @tparam Scalar The type of the element
 * @tparam Matrix The matrix class from which the preconditioner is built;
 * @tparam Vector The vector class to solved with
 */

template<typename SCALAR,class Matrix = Krylov::Matrix<SCALAR>,class Vector = Krylov::Vector<Scalar>> class Preconditioner
{
public:
	using Scalar=SCALAR;
    
};

}//namespace Krylov