/*
 * Vector.hpp
 *
 *  Created on: Dec 30, 2022
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
// To avoid stupid warnings if I do not use openmp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace Krylov
{
/*!
 * A custom vector
 * @tparam Scalar The type of the element
 */

template<typename SCALAR> class Vector
{
public:
	using Scalar=SCALAR;
	
	/*!
   * Constructor may take size of the vector
   * @param size
   */
	Vector(std::size_t size = 0)
		: nValues{size}
	{
		buffer.resize(nValues);
	}

  /*!
   * Sets a new size
   * @param size
   */
  void
  resize(std::size_t size)
  {
    nValues = size;
    buffer.resize(nValues);
  }

	/*!
	 * Get element v[i,j]. Const Version
	 * @param i 
	 * @return the value at index i
	 */
	auto
	operator[](std::size_t i) const
	{
		return buffer[i];
	}

	/*!
	 * Get element v[i,j]. Non const version
	 * It allows to change the element
	 * @param i 
	 * @return the value at index i
	 */
	auto &
	operator[](std::size_t i) 
	{
		return buffer[i];
	}

	/*!
   * Size of the vector
   * @return
   */
  auto
  size() const
  {
    return nValues;
  }

	/*!
   * Returns the buffer containing the elements of the vector (non const
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
   * Returns the buffer containing the elements of the vector(const version)
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
   * Dot product
	 *
   * @param v a vector
   * @return The result of <v,w>
   */
	Scalar dot(Vector<Scalar> const &w) const;

	/*!
   * Norm
	 *
   * @return The Euclidean norm of v
   */
	Scalar norm() const;

	/*!
	 * Vector sum
	 * @param w vector to be summed with
	 * @return The result of v + w
	 */
	Vector operator+(Vector<Scalar> const &w) const;
	Vector& operator+=(Vector<Scalar> const &w);

	/*!
	 * Vector difference
	 * @param w 
	 * @return The result of v - w
	 */
	Vector operator-(Vector<Scalar> const &w) const;
	Vector& operator-=(Vector<Scalar> const &w);

	/*!
   * Clears the vector completely and frees memory
   */
  void clear()
  {
    nValues = 0u;
    buffer.clear();
    buffer.shrink_to_fit();
  }

protected:
	std::size_t					nValues;
	std::vector<Scalar> buffer;
};

/*!
 * To write the vector to the output stream
 * @tparam Scalar
 * @param out
 * @param v
 * @return
 */
template <typename Scalar>
std::ostream &operator<<(std::ostream &out, Vector<Scalar> const &v);

/*!
 * Scalar Vector product
 *
 * @param k scalar to be multiplied with
 * @return The result of k * v
 */
template<typename Scalar,typename Arithmetic_Type>
Vector<Scalar> operator*(Arithmetic_Type const &k,Vector<Scalar> const &v);			
/*
 * ***************************************************************************
 * Definitions
 * ***************************************************************************
 */

template<typename Scalar>
Scalar
Vector<Scalar>::dot(const Vector<Scalar> &w) const 
{
	Scalar res{0};

#pragma omp parallel for reduction(+ : res)
	for(std::size_t i = 0; i < nValues; i++){
		res += buffer[i] * w[i];
	}

	return res;
}

template<typename Scalar>
Scalar
Vector<Scalar>::norm() const
{
	Scalar res{0};
#pragma omp parallel for reduction(+ : res)
	for(std::size_t i = 0; i < nValues; i++)
	{
		res += buffer[i] * buffer[i];
	}

	return std::sqrt(res);
}

template<typename Scalar>
Vector<Scalar>
Vector<Scalar>::operator+(Vector<Scalar> const &w) const
{
	Vector<Scalar> z(nValues);
#pragma omp parallel for
	for(std::size_t i = 0; i < nValues; i++)
	{
		z[i] = buffer[i] + w[i];
	}
	return z;
}

template<typename Scalar>
Vector<Scalar>&
Vector<Scalar>::operator+=(Vector<Scalar> const &w) 
{
#pragma omp parallel for
	for(std::size_t i = 0; i < nValues; i++)
	{
		buffer[i] += w[i];
	}
	return *this;
}

template<typename Scalar>
Vector<Scalar>
Vector<Scalar>::operator-(Vector<Scalar> const &w) const
{
	Vector<Scalar> z(nValues);
#pragma omp parallel for
	for(std::size_t i = 0; i < nValues; i++)
	{
		z[i] = buffer[i] - w[i];
	}
	return z;
}

template<typename Scalar>
Vector<Scalar>&
Vector<Scalar>::operator-=(Vector<Scalar> const &w) 
{
#pragma omp parallel for
	for(std::size_t i = 0; i < nValues; i++)
	{
		buffer[i] -= w[i];
	}
	return *this;
}

template<typename Scalar>
std::ostream &
operator<<(std::ostream &out, Vector<Scalar> const &v)
{
	out << "[";
	for(std::size_t i = 0; i < v.size() - 1; i++)
	{
		out << v[i] << ", ";
	}

	out << v[v.size() - 1] << "]";

	return out;
}

template<typename Scalar,typename Arithmetic_Type>
Vector<Scalar>
operator*(Arithmetic_Type const &k,Vector<Scalar> const &v)
{
	Vector<Scalar> res(v.size());

#pragma omp parallel for
	for(std::size_t i = 0; i < v.size(); i++)
	{
		res[i] = k * v[i];
	}
	
	return res;
}

} //namespace Krylov
#pragma GCC diagnostic pop
#endif