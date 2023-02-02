/*
 * Matrix.hpp
 *
 *  Created on: Oct 14, 2022
 *      Author: forma
 */

#ifndef HH_MATRIX_HH
#define HH_MATRIX_HH
#include <exception>
#include <iostream>
#include <random>
#include <type_traits>
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{
/*!
 * To specify the matrix storage ordering
 */
enum class ORDERING
{
  ROWMAJOR = 0,
  COLUMNMAJOR = 1
};

/*!
 * A full matrix
 * @tparam Scalar The type of the element
 * @tparam Vector The vector class to be multiplied with
 * @tparam ORDER The Storage order (default row-wise)
 */
template <typename SCALAR, ORDERING ORDER = ORDERING::ROWMAJOR> class Matrix
{
public:
  using Scalar=SCALAR;
	/*!
   * Constructor may take number of rows and columns
   * @param nrows Number of rows
   * @param ncols Number of columns
   */
  Matrix(std::size_t nrows = 0, std::size_t ncols = 0)
    : nRows{nrows}, nCols{ncols}
  {
    buffer.resize(nRows * nCols);
  }

  /*!
   * Sets new number of rows and columns
   * @param nrows
   * @param ncols
   */
  void
  resize(std::size_t nrows, std::size_t ncols)
  {
    nRows = nrows;
    nCols = ncols;
    buffer.resize(nRows * nCols);
  }

  /*!
   * Get element A(i,j). Const version
   * @param i
   * @param j
   * @return The value
   */
  auto
  operator()(std::size_t i, std::size_t j) const
  {
    if constexpr(ORDER == ORDERING::ROWMAJOR)
      {
        return buffer[j + i * nCols];
      }
    else
      {
        return buffer[i + j * nRows];
      }
  }
  /*!
   * Get element A(i,j). Non const version
   *
   * It allows to change the element.
   *
   * @param i
   * @param j
   * @return The value
   */

  auto &
  operator()(std::size_t i, std::size_t j)
  {
    if constexpr(ORDER == ORDERING::ROWMAJOR)
      {
        return buffer[j + i * nCols];
      }
    else
      {
        return buffer[i + j * nRows];
      }
  }

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
   * Returns the buffer containing the elements of the matris (non const
   * version)
   *
   * @note I use decltype(auto) becouse I want it to return exactly what
   * std::vector::data() returns.
   *
   * @return A pointer to the buffer
   */
  decltype(auto)
  data() noexcept
  {
    return buffer.data();
  }

  /*!
   * Returns the buffer containing the elements of the matris (const version)
   *
   * @note I use decltype(auto) becouse I want it to return exactly what
   * std::vector::data() returns.
   *
   * @return A pointer to the buffer
   */
  decltype(auto)
  data() const noexcept
  {
    return buffer.data();
  }

  /*!
   * Multiplication with a Vector
   *
   * @tparam Vector 
   * @param v a vector
   * @return The result of A*v
   */
  template<class Vector>
  Vector operator*(Vector const &v) const;

  /*!
   * Multiplication with another Matrix
   *
   * @param matrix a vector
   * @return The result of A*B
   */
  Matrix operator*(Matrix const &B) const;

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
   * Method to extract a row
   *
   * @tparam the row that will be copied 
   * @param v vector where the row is stored
   * @param index the row index
   */
  template<class Vector>
  void extract_row(Vector &v,std::size_t const &index) const;  

   /*!
   * Method to set a rpw
   *
   * @tparam the row that will be set
   * @param v vector from where the row is copied
   * @param index the row index
   */
  template<class Vector>
  void set_row(Vector &v,std::size_t const &index);   

  /*! 
   * Method to compute a minpr
   *
   *  @param k 
   */
  Matrix compute_minor(std::size_t const& index) const;

  /*!
   * Method to compute the transpose
   *
   * @return the transpose of the matrix
   */
  Matrix transpose() const;

  /*!
   * The storage ordering of the matrix
   *
   * @note since c++14 inline is implicit for static constexpr
   */
  static constexpr ORDERING ordering{ORDER};
  /*!
   * To read from a file (or any input stream)
   * @param input The input stream to read from
   */
  void readFromStream(std::istream &input);
 /*!
  * The size of the buffer
  * @return
  */
  auto bufferSize()const {return buffer.size();}
  /*!
   * Fills the matrix with random numbers in the interval [0.,1.) if Scalar floating point
   * or in [-10,10] id scalar is an integral type.
   * @note works only if Scalar is of arithmetic type.
   */
  void fillRandom();

  /*!
   * Clears the matrix completely and frees memory
   */
  void clear()
  {
    nRows=0u;
    nCols=0u;
    buffer.clear();
    buffer.shrink_to_fit();
  }
protected:
  std::size_t         nRows = 0u;
  std::size_t         nCols = 0u;
  int                 nnz = 0u;
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
template <typename Scalar, ORDERING ORDER>
std::ostream &operator<<(std::ostream &out, Matrix<Scalar, ORDER> const &mat);

/*
 * ***************************************************************************
 * Definitions
 * ***************************************************************************
 */

template<typename Scalar, ORDERING ORDER>
template <class Vector>
Vector
Matrix<Scalar,ORDER>::operator*(Vector const &v) const
{
  Vector res(nRows);
  if constexpr(ORDER == ORDERING::ROWMAJOR)
    {
      // loop over rows
      for(std::size_t i = 0; i < nRows; ++i)
        {
          Scalar r{0};
          auto   offset = i * nCols;
          // loop over columns (can be multithreaded)
#pragma omp parallel for shared(i, res, offset) reduction(+ : r)
          for(std::size_t j = 0; j < nCols; ++j)
            r += buffer[offset + j] * v[j];
          
          res[i] = r;
        }
    }
  else
    {
      // loop over the columns
      for(std::size_t j = 0; j < nCols; ++j)
        {
          auto c = v[j];
          auto offset = j * nRows;
          // loop over rows (can be multithreaded)
#pragma omp parallel for shared(j, res, c, offset)
          for(std::size_t i = 0; i < nRows; ++i)
            {
              res[i] += c * buffer[offset + i];
            }
        }
    }
  return res;

}

template<typename Scalar, ORDERING ORDER>
Matrix<Scalar,ORDER>
Matrix<Scalar,ORDER>::operator*(Matrix<Scalar,ORDER> const &B) const
{
  Matrix<Scalar,ORDER> C(nRows,B.cols());
  
  if constexpr(ORDER == ORDERING::ROWMAJOR)
    for(std::size_t i = 0; i < nRows; i++)
      for(std::size_t j = 0; j < B.cols(); j++)
        for(std::size_t k = 0; k < nCols; k++)
          C(i,j) += buffer[k + i  * nCols] * B(k,j);
  
  return C;
}

template<typename Scalar, ORDERING ORDER>
template <class Vector>
void
Matrix<Scalar,ORDER>::extract_column(Vector &v,std::size_t const &index) const
{

#pragma omp parallel for
  for(std::size_t i = 0; i < nRows; i++)
  {
    v[i] = this->operator()(i,index);
  }
}

template<typename Scalar, ORDERING ORDER>
template <class Vector>
void
Matrix<Scalar,ORDER>::extract_row(Vector &v,std::size_t const &index) const
{

#pragma omp parallel for
  for(std::size_t j = 0; j < nCols; j++)
  {
    v[j] = this->operator()(index,j);
  }
}

template<typename Scalar, ORDERING ORDER>
template <class Vector>
void
Matrix<Scalar,ORDER>::set_row(Vector &v,std::size_t const &index) 
{

#pragma omp parallel for
  for(std::size_t j = 0; j < nCols; j++)
  {
    this->operator()(index,j) = v[j];
  }
}

template<typename Scalar, ORDERING ORDER>
Matrix<Scalar,ORDER>
Matrix<Scalar,ORDER>::transpose() const 
{
  Matrix<Scalar,ORDER> T(nRows,nCols);
  for(std::size_t i = 0; i < nRows; i++)
  {
    for(std::size_t j = 0; j < nCols; j++)
    {
      T(i,j) = this->operator()(j,i);
    }
  }
  return T;
}

template<typename Scalar, ORDERING ORDER>
Matrix<Scalar,ORDER>
Matrix<Scalar,ORDER>::compute_minor(std::size_t const& index) const
{
  Matrix<Scalar,ORDER> minor(nRows,nCols);
  for(std::size_t i = 0; i < index; i++)
    minor(i,i) = 1;

  for(std::size_t i = index; i < nRows; i++)
    for(std::size_t j = index;j < nCols; j++)
      minor(i,j) = this->operator()(i,j);
  
  return minor;
}

template <typename Scalar, ORDERING ORDER>
void
Matrix<Scalar, ORDER>::readFromStream(std::istream &input)
{
  input >> nRows >> nCols;
  resize(nRows, nCols);
  for(std::size_t i = 0u; i < nRows; ++i)
    for(std::size_t j = 0u; j < nCols; ++j)
      input >> this->operator()(i, j);

  if(!input.good())
    {
      throw std::runtime_error("ERROR: problems while reading from stream\n");
    }
}

template <typename Scalar, ORDERING ORDER>
std::ostream &
operator<<(std::ostream &out, Matrix<Scalar, ORDER> const &mat)
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

template <typename Scalar, ORDERING ORDER>
void
Matrix<Scalar, ORDER>::fillRandom()
{
  static_assert(std::is_arithmetic_v<Scalar>,
                "fillRandom requires elements of  arithmetic type");
  std::random_device         r;
  std::default_random_engine e1(r());
  if constexpr(std::is_integral_v<Scalar>)
    {
      std::uniform_int_distribution<Scalar> dist(-10, 10);
      for(auto &v : buffer)
        {
          v = dist(e1);
        }
    }
  else
    {
      std::uniform_real_distribution<Scalar> dist(0., 1.);
      for(auto &v : buffer)
        {
          v = dist(e1);
        }
    }
}

} // namespace Krylov
#pragma GCC diagnostic pop

#endif 