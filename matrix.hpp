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
#include "vector.hpp"
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
template <typename SCALAR,class Vector = Krylov::Vector<SCALAR>, ORDERING ORDER = ORDERING::ROWMAJOR> class Matrix
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
   *It allows to change the element.
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
   * Multiplication with a std::vector
   *
   * @note to be complete I should also add the multiplication with a matrix
   * with just one culumn
   *
   * @param v a vector
   * @return The result of A*v
   */
  Vector operator*(Vector const &v) const;

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
template <typename Scalar,class Vector, ORDERING ORDER>
std::ostream &operator<<(std::ostream &out, Matrix<Scalar,Vector, ORDER> const &mat);

/*
 * ***************************************************************************
 * Definitions
 * ***************************************************************************
 */
template <typename Scalar,class Vector, ORDERING ORDER>
Vector
Matrix<Scalar,Vector, ORDER>::operator*(const Vector &v) const
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

template <typename Scalar,class Vector, ORDERING ORDER>
void
Matrix<Scalar, Vector, ORDER>::readFromStream(std::istream &input)
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

template <typename Scalar,class Vector, ORDERING ORDER>
std::ostream &
operator<<(std::ostream &out, Matrix<Scalar,Vector, ORDER> const &mat)
{
  out << "[" << std::endl;
  for(std::size_t i = 0u; i < mat.rows(); ++i)
    {
      for(std::size_t j = 0u; j < mat.cols(); ++j)
        out << mat(i, j) << ", ";
      out << std::endl;
    }
  out << "]" << std::endl;
  return out;
}

template <typename Scalar,class Vector, ORDERING ORDER>
void
Matrix<Scalar, Vector, ORDER>::fillRandom()
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