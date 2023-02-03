#ifndef HH_MMREADER_HH
#define HH_MMREADER_HH
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../base_classes/sparse_matrix.hpp"
namespace Krylov
{


template<typename Scalar>
SparseMatrix<Scalar> readCSRMatrix(const std::string &filename)
{
    
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  // Read header line
  std::string line;
  std::getline(file, line);

	auto symmetric = line.find("symmetric") != std::string::npos ? true : false;

  // Skip comments
  do {
    std::getline(file, line);
  } while (line[0] == '%');

  // Read dimensions
  std::stringstream ss(line);
  std::size_t rows, cols, nonzeros;
  ss >> rows >> cols >> nonzeros;

  SparseMatrix<Scalar> matrix(rows,cols,nonzeros + nonzeros * symmetric);

	while (std::getline(file, line)) {
    std::size_t row, col;
    double value;
    std::stringstream ss(line);
    ss >> row >> col >> value;

    // Matrix Market uses 1-based indexing, convert to 0-based
    row--;
    col--;
		matrix.set(row,col) = value;
		if(symmetric) matrix.set(col,row) = value;
	}

	return matrix;
 
}

}//namspace Krylov
#endif