/*
 * ILU.hpp
 *
 *  Created on: Jan 14, 2023
 *      Author: bolt3x
 */

#ifndef HH_ILU_HH
#define HH_ILU_HH
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
 * An ilu preconditioner class
 * @tparam Scalar element
 * @tparam Vector The vector class to "solve" with
 */

template<typename SCALAR, class Matrix = Krylov::SparseMatrix<SCALAR>, class Vector = Krylov::Vector<SCALAR>> class ILUPreconditioner
{ 
public:
	ILUPreconditioner(const Matrix &A)
}
}//namespace Krylov
#endif