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
#include <map>
#include <algorithm>
#include "../direct_solvers/qr_solver.hpp"
#include "../base_classes/sparse_matrix.hpp"
#include "../utils/plotter.hpp"


// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{
/*!
 * To specify how the spartsity pattern of the Inverse is chosen
 */
enum class PATTERN
{
  STATIC = 0,
  DYNAMIC = 1
};

enum class SPARSITY
{
	LOW = 1,
	MEDIUM = 2,
	HIGH = 4
};
/*!
 * An identity "preconditioner" class
 * @tparam Scalar element
 * @tparam Pattern how the pattern of the inverse will be construtcted
 */

template<typename SCALAR,PATTERN Sparsity = PATTERN::STATIC> class SpaiPreconditioner
{
public:
	using Scalar=SCALAR;
	using Matrix = Krylov::Matrix<Scalar>;
	using Vector = Krylov::Vector<Scalar>;
	using SparseMatrix = Krylov::SparseMatrix<Scalar>;

	/*! Constructor takes a sparse matrix
	 *  @tparam Matrix
	 *  @param A matrix from which we compute the inverse
	 */
	template<class MatrixType>
	SpaiPreconditioner(MatrixType const &A,SPARSITY const &sp_ = SPARSITY::MEDIUM) : M(A.rows(),A.cols()),sp(sp_)
	{	

		staticPattern(A);
	
	
		M.rowPtrs.resize(1);
		compute(A,tol);
	}

	template<class MatrixType>
	void compute(MatrixType const &A,Scalar &tol);


	/*!
	 * this method is used in the dynamic pattern
	 * case to update the pattern
	 * @param A
	 * @param J the pattern to be updated
	 * @param r the residual is used to compute the indexes to be added to J
	 */
	template<class MatrixType>
	void update(MatrixType const &A,std::vector<std::size_t> &J,Krylov::Vector<Scalar> const &r);

	/*! this method computes M from
	 * a vector of vectors m
	 * @param m a std::vector of Vectors
	 */
	void build(std::vector<Vector> const &m);

	/*!
	 * solve method to apply the preconditioner 
	 * @param v 
	 * @return the preconditioned vector
	 */
	template<class Vector>
	Vector solve(Vector const &v) const { return M*v; }

	/*!
	 * this method return the computed 
	 * preconditioner M
	 */
	SparseMatrix getM() const {return M;}

	/*!
	 * this method sets the sparsity pattern
	 * based on A and a tolerance;
	 * @param A
	 */
	template<class MatrixType>
	void staticPattern(MatrixType const &A);
	/*!
	 * this method sets the sparsity pattern to diagonal,
	 * in case of dynamic sparsity pattern, it will start from here
	 * @param diag the lenght of the diagonal
	 */
	void diagPattern(std::size_t const &d);

	void setTol(Scalar &tol_){tol = tol_;}
	void setSparsity(SPARSITY &sp_){sp = sp_;}

protected:
	std::map<std::size_t,std::vector<std::size_t>> pattern;
	SparseMatrix M;
	Scalar tol = 0.2;
	SPARSITY sp;
};

template<typename Scalar,PATTERN Sparsity>
template<class MatrixType>
void
SpaiPreconditioner<Scalar,Sparsity>::staticPattern(MatrixType const &A)
{
	for(std::size_t i = 0; i < A.rows(); i++)
	{
		for(std::size_t j = 0; j < A.cols(); j++)
		{	
			Scalar di = std::sqrt(std::abs(1 / (A(i,i) + (1 * !A(i,i) ))));
			Scalar dj = std::sqrt(std::abs(1 / (A(j,j) + (1 * !A(j,j) ))));
			
			if(i == j || std::abs(di * A(i,j) * dj) > tol)
			{
				pattern[i].push_back(j);
			}
		}
	}
}
template<typename Scalar,PATTERN Sparsity>
void
SpaiPreconditioner<Scalar,Sparsity>::diagPattern(std::size_t const &d)
{
	for(std::size_t i = 0; i < d; i++)
	{
		pattern[i].push_back(i);
	}
}
template<typename Scalar,PATTERN Sparsity>
template<class MatrixType>
void
SpaiPreconditioner<Scalar,Sparsity>::compute(MatrixType const &A,Scalar &tol)
{	
	std::vector<Vector> m(A.rows());
#pragma omp parallel for
	for(std::size_t k = 0; k < A.rows();k++)
	{

		std::vector<std::size_t> J = pattern[k];
		double res = 1;
		while(1)
		{	
			std::vector<std::size_t> I;

			for(std::size_t i = 0; i < A.rows(); i++)
			{
				for(std::size_t j = 0; j < J.size(); j++)
				{
					if(A(i,J[j]) && std::find(I.begin(),I.end(),i) == I.end()) I.push_back(i);
				}
			}

			
			Vector reducedE_k(I.size());
			Vector e_k(A.rows());
			e_k[k] = 1;
			for(std::size_t i = 0; i < I.size(); i++)
			{
				reducedE_k[i] = e_k[I[i]];
			}	
			Matrix reducedA(I.size(),J.size());
			for(std::size_t i = 0; i < reducedA.rows(); i++)
			{
				for(std::size_t j = 0; j < reducedA.cols(); j++)
				{
					reducedA(i,j) = A(I[i],J[j]);
				}
			}
			
			if(reducedA.rows() > A.rows())
				break;
			

			if(reducedA.rows() > 1)
			{
				QRSolver<Scalar> qr(reducedA);
				
				m[k] = qr.solve(reducedE_k);
				
			}
			
			if(reducedA.rows() == 1)
			{
				m[k] = Vector(1,0);
				for(std::size_t i = 0; i < m[k].size(); i++)
				{
					if(reducedA(0,i))
					{
						m[k][0] = 1/reducedA(0,i);
					}
				}
				break;
			}
			if(!reducedA.rows() || !reducedA.cols())
				continue;

			Matrix AJ(A.rows(),J.size());
			for(std::size_t i  = 0; i < A.rows(); i++)
			{
				for(std::size_t j = 0; j < J.size(); j++)
				{
					AJ(i,j) = A(i,J[j]);
				}	
			}
			Vector r = AJ * m[k] - e_k ;

			res = r.norm();
			
			Scalar sparsity = Scalar(sp); 
			if(res < tol || m[k].size() >= A.nonzero() / (A.rows()* sparsity) || Sparsity == PATTERN::STATIC)
			{
				
				pattern[k] = J;
			
				break;
			}

			update(A,J,r);
			
		}
		
	}

	build(m);
		
}

template<typename Scalar,PATTERN Sparsity>
template<class MatrixType>
void
SpaiPreconditioner<Scalar,Sparsity>::update(MatrixType const &A,std::vector<std::size_t> &J,Vector const &r)
{
	std::vector<std::size_t> L;
	std::vector<std::size_t> Jnew;
	for(std::size_t i = 0; i < r.size(); i++)
	{
		if(r[i] && std::find(J.begin(),J.end(),i) == J.end()) 		
			Jnew.push_back(i); 
	}			
	double p_j = r.norm();
	double min = -1;
	for(auto &j : Jnew)
	{
		Vector e_j(A.cols());
		e_j[j] = 1;

		double p_jnew = r.norm() - std::pow(r.dot(A*e_j),2)/((A*e_j).norm());
		
		if(p_jnew < p_j)
		{	
			p_j	= p_jnew;
			min = j;
		}
		
	}
	
	J.push_back(min);
	std::sort(J.begin(),J.end());
	
} 

template<typename Scalar,PATTERN Sparsity>
void
SpaiPreconditioner<Scalar,Sparsity>::build(std::vector<Vector> const &m)
{		
	
	for(std::size_t k = 0; k <  pattern.size(); k++)
	{
		std::vector<std::size_t> J = pattern[k];
		Vector reducedM_k = m[k];
		
		for(std::size_t i = 0; i < J.size(); i++)
		{
			M.colInd.push_back(J[i]);
			if(!std::isnan(reducedM_k[i]))
				M.buffer.push_back(reducedM_k[i]);
			else 
				M.buffer.push_back(1);
			M.nnz++;
		}
		M.rowPtrs.push_back(M.nnz);
		
	}
}

}//namespace Krylov
#endif
