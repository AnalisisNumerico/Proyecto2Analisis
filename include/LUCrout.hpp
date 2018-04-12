/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi {

  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of U
   * is composed by 1's.
   */
  template<typename T>
  void unpackCrout(const Matrix<T>& LU,
                   Matrix<T>& L,
                   Matrix<T>& U) {

    if (LU.cols() == LU.rows()) {

      int luCols = LU.cols();
      int luRows = LU.rows();

      anpi::Matrix<T> u(LU.rows(),LU.cols());

      for(int i = 0; i < luRows; i++) { //diagonal
        u[i][i] = T(1);
      }

      for(int j = 1; j < luCols; j++) {
        for(int i = 0; i < j; i++) {
          u[i][j] = LU[i][j];
        }
      }

      U = u;

      anpi::Matrix<T> l(LU.rows(),LU.cols());

      for(int i = 0; i < luRows; i++) {
        for(int j = 0; j <= i; j++) {
          l[i][j] = LU[i][j];
        }
      }

      L = l;

    }
    else {
      throw anpi::Exception("Unpack Crout: Invalid compressed matrix size");
    }

  }

  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.  
   *
   * Crout's way of packing assumes a diagonal of
   * 1's in the U matrix.
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
  template<typename T>
  void luCrout(const Matrix<T>& A,
               Matrix<T>& LU,
               std::vector<size_t>& permut) {

    if(A.rows() == A.cols()) {

      LU = A;
      int n = A.rows();
      std::vector<size_t > index (A.rows());

      for(int i = 0; i < n; i++) { //relleno vector indice
        index[i] = i;
      }

      for(int j = 1; j < n; j++) {
        LU[0][j] /= LU[0][0];
      }

      T sum;

      for (int j = 0; j < n-1; j++) {

        int bigI = j;
        for(int i = j; i < n; i++) { //busco el mayor numero en la columna
          if(std::abs(LU[bigI][j]) < std::abs(LU[i][j])) { //cambio el indice de la fila
            bigI = i;
          }
        }
        if(bigI != j) { //si se encontro un pivote mayor se intercambian filas
          T matrixTmp;
          for(int k = 0; k < n; k++) { //intercambio fila en matriz LU
            matrixTmp = LU[j][k];
            LU[j][k] = LU [bigI][k];
            LU[bigI][k] = matrixTmp;
          }
          int indexTmp;
          indexTmp = index[j];    //intercambio en vector indice
          index[j] = index[bigI];
          index[bigI] = indexTmp;
        }

        if(j > 0) {
          for (int i = j; i < n; i++) {
            sum = T(0);
            for (int k = 0; k < j; k++) {
              sum += (LU[i][k] * LU[k][j]);
            }
            LU[i][j] -= sum;
          }
          for (int i = j + 1; i < n; i++) {
            sum = T(0);
            for (int k = 0; k < j; k++) {
              sum += (LU[j][k] * LU[k][i]);
            }
            LU[j][i] = (LU[j][i] - sum) / LU[j][j];
          }
        }
      }

      permut = index;

      sum = T(0);

      for(int k = 0; k < n - 1; k++) {
        sum += LU[n-1][k] * LU[k][n-1];
      }

      LU[n-1][n-1] -= sum;


/*
      LU = A;
      int n = A.rows();
      std::vector<size_t > index (A.rows());

      for(int i = 0; i < n; i++) { //relleno vector indice
        index[i] = i;
      }

      T sum;

      for(int j = 0; j < n; j++) {

        for(int i = 0; i < j; i++) {
          sum = LU[i][j];
          for(int k = 0; k < i; k++) {
            sum -= LU[i][k] * LU[k][j];
          }
          LU[i][j] = sum;
        }

        int bigI = j;

        for(int i = j; i < n; i++) { //busco el mayor numero en la columna
         
          sum = LU[i][j];
          for(int k = 0; k < j; k++) {
            sum -= LU[i][k] * LU[k][j];
          }
          LU[i][j] = sum;

          if(std::abs(LU[bigI][j]) < std::abs(LU[i][j])) { //cambio el indice de la fila
            bigI = i;
          }
        }
        if(bigI != j) { //si se encontro un pivote mayor se intercambian filas
          T matrixTmp;
          for(int k = 0; k < n; k++) { //intercambio fila en matriz LU
            matrixTmp = LU[j][k];
            LU[j][k] = LU [bigI][k];
            LU[bigI][k] = matrixTmp;
          }
          size_t indexTmp;
          indexTmp = index[j];    //intercambio en vector indice
          index[j] = index[bigI];
          index[bigI] = indexTmp;
        }

        if (LU[j][j] == 0) {
          throw anpi::Exception("Crout: zero division");
        }

        if(j < (n-1)) {
          sum = T(1)/(LU[j][j]);
          for(int i = j+1; i < n; i++) {
            LU[i][j] *= sum;
          }
        }


      }

      permut = index;

*/
/*
      LU = A;
      int n = A.rows();
      std::vector<size_t > index (A.rows());

      for(int i = 0; i < index.size(); i++) { //relleno vector indice
        index[i] = i;
      }

      for(int j = 0; j < n-1; j++) {
        int bigI = j;
        for(int i = j; i < n; i++) { //busco el mayor numero en la columna
          if(std::abs(LU[bigI][j]) < std::abs(LU[i][j])) { //cambio el indice de la fila
            bigI = i;
          }
        }
        if(bigI != j) { //si se encontro un pivote mayor se intercambian filas
          T matrixTmp;
          for(int k = 0; k < n; k++) { //intercambio fila en matriz LU
            matrixTmp = LU[j][k];
            LU[j][k] = LU [bigI][k];
            LU[bigI][k] = matrixTmp;
          }
          int indexTmp;
          indexTmp = index[j];    //intercambio en vector indice
          index[j] = index[bigI];
          index[bigI] = indexTmp;
        }
      }

      permut = index;

      for(int j = 1; j < n; j++) {
        LU[0][j] /= LU[0][0];
      }

      T sum;

      for (int j = 1; j < n-1; j++) {

        for (int i = j; i < n; i++) {
          sum = T(0);
          for (int k = 0; k < j; k++) {
            sum += (LU[i][k] * LU[k][j]);
          }
          LU[i][j] -= sum;
        }
        for (int i = j + 1; i < n; i++) {
          sum = T(0);
          for(int k = 0; k < j; k++) {
            sum += (LU[j][k] * LU[k][i]);
          }
          if (LU[j][j] == 0) {
            throw anpi::Exception("Crout: zero division");
          }
          LU[j][i] = (LU[j][i] - sum) / LU[j][j];
        }
      }

      sum = T(0);

      for(int k = 0; k < n - 1; k++) {
        sum += LU[n-1][k] * LU[k][n-1];
      }

      LU[n-1][n-1] -= sum;
*/
    }
    else {
      throw anpi::Exception("Crout: invalid decomposition matrix size");
    }

  }

}

#endif

