/*
 * QRSolver.hpp
 *
 *  Created on: Jan 18, 2023
 *      Author: bolt3x
 */

#ifndef HH_QRSOLVER_HH
#define HH_QRSOLVER_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <cmath>
#include "../base_classes/matrix.hpp"
#include "../base_classes/vector.hpp"
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{

template<typename SCALAR> class QRSolver
{
public:
	using Scalar = SCALAR;
	using Matrix = Krylov::Matrix<Scalar>;
	using Vector = Krylov::Vector<Scalar>;

	template<class MatrixType>
	QRSolver(MatrixType const &A) : Q(A.rows(),A.rows()), R(A.rows(),A.cols())
	{
		compute(A);
	}

	template<class MatrixType>
	void compute(MatrixType const &A);

	Vector solve(Vector const &b) const;

template<class VectorType>
Matrix householder(VectorType const &v) const;

	Matrix getQ()
	{
		return Q;
	}

	Matrix getR()
	{
		return R;
	}
protected:
	Matrix Q;
	Matrix R;
};

template<typename Scalar>
template<class MatrixType>
void
QRSolver<Scalar>::compute(MatrixType const &A)
{
	std::size_t m = A.rows();
	std::size_t n = A.cols();
	
	std::vector<Matrix> qv(m);
		
	Matrix z(m,n);
	for(std::size_t i = 0; i < m; i++)
		for(std::size_t j = 0; j < n; j++)
			z(i,j) = A(i,j);
		
	Matrix z1;
	
	for(std::size_t k = 0; k < n && k < m -1; k++)
	{
		Vector e(m), x(m);
		Scalar a;
		z1 = z.compute_minor(k);

		z1.extract_column(x,k);

		a = x.norm();
		if(A(k,k) > 0) a = -a;

		for(std::size_t i = 0; i < e.size(); i++)
			e[i] = i == k ? 1 : 0;
			
		e = x + a * e;
			
		Scalar enorm = e.norm();
			for(std::size_t i = 0; i < e.size(); i++)
				e[i] = e[i] / enorm;
			qv[k] = householder(e);
			z = qv[k] * z1;
		}
		Q = qv[0];

  		for (int i = 1; i < n && i < m - 1; i++) {

    		z1 = qv[i] * Q;
    		Q = z1;
			}

  	R = Q * A;
  	Q = Q.transpose();
	
}
	
template<typename Scalar>
Krylov::Vector<Scalar>
QRSolver<Scalar>::solve(Krylov::Vector<Scalar> const &b) const
{
		
	Matrix I(R.cols(),Q.rows());
	for(std::size_t i = 0; i < R.cols(); i++)
	{
		for(std::size_t j = 0; j < Q.rows(); j++)
		{
			if(i == j) I(i,j) = 1;	
		}
	}
	Krylov::Vector<Scalar> b_ =  I * Q.transpose() * b;
	Krylov::Vector<Scalar> res(b_.size());


#pragma omp parallel for
		for(std::size_t i = 0; i < b_.size(); i++)
		{
			res[i] = b_[i];
		}
		for(int i = b_.size() -1; i >= 0; i--)
		{
			res[i] = res[i] / R(i,i);

#pragma omp parallel for
			for(int j = 0; j < i; j++)
			{
				res[j] -= R(j,i) * res[i];
			}
		}
		return res;
	}

template<typename Scalar>
template<class VectorType>
	Krylov::Matrix<Scalar>
	QRSolver<Scalar>::householder(VectorType const &v) const
	{	
		std::size_t n = v.size();
		Krylov::Matrix<Scalar> H(n,n);
		for (std::size_t i = 0; i < n; i++)
		{
			
    		for (std::size_t  j = 0; j < n; j++)
			{
      			H(i,j) = -2 *  v[i] * v[j];
			}
		}

		for (std::size_t  i = 0; i < n; i++)
    		H(i,i) += 1; 
		return H;
		
	}
}//namespace Krylov
#endif