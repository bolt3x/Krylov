#ifndef HH_PLOTTER_HH
#define HH_PLOTTER_HH
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../base_classes/sparse_matrix.hpp"
namespace Krylov
{


template<typename Scalar>
int plotCSRMatrix(SparseMatrix<Scalar> const &matrix,std::string const &filename)
{
  
  std::ofstream datafile("csr_matrix.dat");
  for (std::size_t i = 0; i < matrix.rows(); i++) {
    for (std::size_t j = matrix.rowPtrs[i]; j < matrix.rowPtrs[i+1]; j++) {
      std::size_t col = matrix.colInd[j];
      datafile << i << " " << col << " " << matrix.buffer[j] << std::endl;
    }
  }
  
  datafile.close();

  std::ofstream gnufile("csr_matrix.gnu");
  gnufile << "set term png" << std::endl;
  gnufile << "set output '" + filename + ".png'" << std::endl;
  gnufile << "plot 'csr_matrix.dat' with points pt 7 lc rgb 'red'" << std::endl;
  gnufile.close();

  std::string cmd = "gnuplot csr_matrix.gnu";
  int result = system(cmd.c_str());
  return result;
}

}//namspace Krylov
#endif