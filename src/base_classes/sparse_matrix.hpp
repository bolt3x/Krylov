/*
 * SparseMatrix.hpp
 *
 *  Created on: Jan 11, 2023
 *      Author: bolt3x
 */

#ifndef HH_SPARSE_MATRIX_HH
#define HH_SPARSE_MATRIX_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{

	/*!
	 * A sparse matrix in CSR format, this is a STATIC storage system
 	 * @tparam Scalar The type of the element
   * @tparam Vector The vector class to be multiplied with
	 */
	template<typename SCALAR> class SparseMatrix
	{
	public:
		using Scalar = SCALAR;
	
	/*!
   * Constructor may take number of rows and columns
	 * @param rows number of rows
	 * @param col number of columns
	 */
  SparseMatrix(std::size_t rows = 0, std::size_t cols = 0)
    : nRows{rows}, nCols{cols}, nnz{0}
  {	
    rowPtrs.resize(nRows + 1);
  }	

	/*!
     * Constructor may take number of rows and columns and number of non zero elements
	 * @param nnzs non zero elements
	 * @param rows number of rows
	 * @param col number of columns
	 */
	SparseMatrix(std::size_t rows, std::size_t cols,std::size_t nnz_)
    : nRows{rows}, nCols{cols}, nnz{nnz_}
  	{	
    rowPtrs.resize(nRows + 1);
	buffer.resize(nnz);
	colInd.resize(nnz);
  	}	

	auto 
	operator()(std::size_t i,std::size_t j) const;

	auto&
	set(std::size_t i,std::size_t j);
	/*!
   * Number of rows
   * @return
   */
  auto
  rows() const
  {
    return nRows;
  };
  /*!
   * Number of columns
   * @return
   */
  auto
  cols() const
  {
    return nCols;
  };
  /*!
   * Number of columns
   * @return
   */
  auto
 nonzero() const
  {
    return nnz;
  };
	/*!
   * Multiplication with a std::vector
   *
   * @note to be complete I should also add the multiplication with a matrix
   * with just one culumn
   *
	 * @tparam Vector class
   * @param v a vector
   * @return The result of A*v
   */
	template<class Vector>
  Vector operator*(Vector const &v) const;
  /*!
   * Multiplication with another Sparse Matrix
   *
   * @param matrix a sparse matrix
   * @return The result of A*B
   */
  SparseMatrix operator*(SparseMatrix const &B) const;
  /*!
   * Method to extract a column
   *
   * @tparam the column that will be copied 
   * @param v vector where the column is stored
   * @param index the column index
   */
  template<class Vector>
  void extract_column(Vector &v,std::size_t const &index) const;

	/*!
   * Method to compute the transpose
   *
   * @return the transpose of the matrix
   */
   SparseMatrix transpose() const;

	public:
		std::size_t nRows;
		std::size_t nCols;
		std::size_t nnz;
		std::vector<int> rowPtrs;
		std::vector<int> colInd;
		std::vector<Scalar> buffer;
};
/*!
 * To write the matrix the output stream
 * @tparam Scalar
 * @tparam ORDER
 * @param out
 * @param mat
 * @return
 */
template <typename Scalar>
std::ostream &operator<<(std::ostream &out, SparseMatrix<Scalar> const &mat);

/*
 * ***************************************************************************
 * Definitions
 * ***************************************************************************
 */

template <typename Scalar>
auto 
SparseMatrix<Scalar>::operator()(std::size_t i,std::size_t j) const
	{
		int pos = rowPtrs[i];	
		int currCol = -1;
		for (; pos < rowPtrs[i + 1]; pos++) 
		{	
			
			currCol = colInd[pos];
			if(currCol == j)
				return buffer[pos];
		
		}
		return Scalar();
	}

template <typename Scalar>
auto&
SparseMatrix<Scalar>::set(std::size_t i,std::size_t j)
{
		
		int pos = rowPtrs[i] - 1 ;
		int currCol = -1;
		
		for (; pos < rowPtrs[i + 1] - 1 ; pos++) 
		{
			
			currCol = colInd[pos];
		
			if (currCol > j) 
			{
				break;
			}
		}


		if (currCol != j) 
		{
			for (std::size_t row = i+1; row <= nRows; row++) 
			{
				rowPtrs[row] += 1;
			}

			if (nnz == 0)
			{	
				colInd.push_back(j);
				buffer.resize(1);
				nnz += 1;
				return buffer[0];
			} 
			else 
			{
				buffer.insert(buffer.begin() + pos + 1, 0);
				colInd.insert(colInd.begin() + pos + 1, j);
				nnz +=1;
				return buffer[pos + 1];
			}
		} 
		else 
		{
			return buffer[pos + 1];
		}
	}

template <typename Scalar>
template<class Vector>
Vector
SparseMatrix<Scalar>::operator*(Vector const &v) const
{
	Vector res(nRows);
	for(std::size_t i = 0; i < nRows; ++i)
	{
		Scalar r{0};
#pragma omp parallel for shared(i, res) reduction(+ : r)
		for(std::size_t j = rowPtrs[i]; j < rowPtrs[i+1] ;j++)
		{
			std::size_t index = colInd[j];
			r += buffer[j] * v[index];
			
		}
		res[i] = r;
	}
	return res;
}

template<typename Scalar>
SparseMatrix<Scalar>
SparseMatrix<Scalar>::operator*(SparseMatrix<Scalar> const &B) const
{
  SparseMatrix<Scalar> C(nRows,B.cols());
  
for(std::size_t i = 0; i < nRows; i++)
    for(std::size_t j = 0; j < B.cols(); j++)
    for(std::size_t k = 0; k < nCols; k++)
		if(this->operator()(i,k) && B(k,j))
          C.set(i,j) += this->operator()(i,k) * B(k,j);
  
  return C;
}


template<typename Scalar>
template <class Vector>
void
SparseMatrix<Scalar>::extract_column(Vector &v,std::size_t const &index) const
{

#pragma omp parallel for
  for(std::size_t i = 0; i < nRows; i++)
  {
    v[i] = this->operator()(i,index);
  }
}

template <typename Scalar>
SparseMatrix<Scalar>
SparseMatrix<Scalar>::transpose() const
{
	SparseMatrix<Scalar> res(nCols,nRows,nnz);
	res.rowPtrs.resize(res.rowPtrs.size()+1);
	for(std::size_t i = 0; i < nnz; i++)
	{
		res.rowPtrs[colInd[i] + 2]++;
	}
    for (std::size_t i = 2; i < res.rowPtrs.size(); ++i) 
	{
        res.rowPtrs[i] += res.rowPtrs[i - 1];
    }

    for (std::size_t i = 0; i < nRows; i++) {
        for (std::size_t j = rowPtrs[i]; j < rowPtrs[i + 1]; ++j) {
            const int new_index = res.rowPtrs[colInd[j] + 1]++;
            res.buffer[new_index] = buffer[j];
            res.colInd[new_index] = i;
        }
    }
    res.rowPtrs.pop_back(); 

    return res;
}

template <typename Scalar>
std::ostream &
operator<<(std::ostream &out, SparseMatrix<Scalar> const &mat)
{
  out << "[" << std::endl;
  for(std::size_t i = 0u; i < mat.rows(); ++i)
    {
      for(std::size_t j = 0u; j < mat.cols(); ++j)
        out << mat(i, j) << ", ";
      out << std::endl;
    }
  out << "]";
  return out;

}

}//namespace Krylov
#endif