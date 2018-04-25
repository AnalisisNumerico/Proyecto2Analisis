#ifndef ANPI_APLICACION_HPP
#define ANPI_APLICACION_HPP

#include "Matrix.hpp"
#include "Solver.hpp"
#include <opencv2/opencv.hpp>

namespace anpi {

  void mapper(const size_t    rows,
              const size_t    cols,
              const int          n,
              const int          m,
              const int          i,
              const int          j,
              int&               x) {

    if (n < rows && i < rows && m < cols && j < cols) {
      if ((n - i) == 0 && std::abs(m - j) == 1) { //horizontal
        int min = m;
        if (j < m) {
          min = j;
        }

        x = (cols - 1) * n + cols * n + min;

      } else if ((m - j) == 0 && std::abs(n - i) == 1) { //vertical
        int min = n;
        if (i < n) {
          min = i;
        }

        x = (cols - 1) * (min + 1) + min * cols + m;

      } else {
        throw anpi::Exception("anpi::mapper error: Invalid nodes coordinates");
      }

    } else {
      throw anpi::Exception("anpi::mapper error: Nodes coordinates out of matrix bounds");
    }
  }

  void inverseMapper(const size_t rows,
                     const size_t cols,
                     const int&      x,
                     int&            n,
                     int&            m,
                     int&            i,
                     int&            j) {
    
    if(x >= (2 * rows * cols - (rows + cols))){
      throw anpi::Exception("anpi::inverse mapper error: Index out of bounds");
    }
    else {
    
      int resistIndex = 0;
      int k = 0;
      bool loop  = true;
      while(loop){
        resistIndex += (2 * cols) - 1;
        if(x < resistIndex) {
          n = k;
          resistIndex -= cols;
          if(x >= resistIndex) { //vertical
            m = x - resistIndex;
            i = n + 1;
            j = m;
          }
          else{ //horizontal
            resistIndex -= (cols - 1);
            m = x - resistIndex;
            i = n;
            j = m + 1;
          }
          loop  = false;
        }
        k++;
      }
    }
    
  }
/*
  template<typename T>
  void matrixFiller(const size_t     rows,
                    const size_t     cols,
                    anpi::Matrix<T>&    A) {
    anpi::Matrix<T> a(rows * cols,(2 * rows * cols) - (rows + cols));
    int count = 0;
    int j;
    for(int n = 0; n < rows; n++) {
      for(int m = 0; m < cols; m++) {
        if(m - 1 >= 0) {
          anpi::mapper<float>(rows, cols, n, m, n, m - 1, j);
          a[count][j] = T(1);
        }
        if(n - 1 >= 0) {
          anpi::mapper<float>(rows, cols, n, m, n - 1, m, j);
          a[count][j] = T(1);
        }
        if(m + 1 < cols) {
          anpi::mapper<float>(rows, cols, n, m, n, m + 1, j);
          a[count][j] = T(-1);
        }
        if(n + 1 < rows) {
          anpi::mapper<float>(rows, cols, n, m, n + 1, m, j);
          a[count][j] = T(-1);
        }
        count++;
      }
    }
    A = a;
  }
*/
  template<typename T>
  void vectorFiller(const size_t rows,
                    const size_t cols,
                    const int       n,
                    const int       m,
                    const int       i,
                    const int       j,
                    std::vector<T>& b) {
    if(n >= rows || m >= cols ||
       i >= rows || j >= cols ||
       n < 0 || m < 0 ||
       i < 0 || j < 0) {
      throw anpi::Exception("anpi::vector filler error: Index out of bounds");
    }
    else {
      b.clear();
      b.resize(2 * rows * cols - (rows + cols));
      b[cols * n + m] = T( 1);
      b[cols * i + j] = T(-1);
    }
  }

  template<typename T>
  void mapCreator(std::string      path,
                  anpi::Matrix<T>&  map) {
    cv::Mat image = cv::imread(path,0);
    map = anpi::Matrix<T>(image.rows,image.cols);
    int pixelColor;
    for (int x = 0;x < image.rows; x++) {
      for (int y = 0; y < image.cols; y++) {
        pixelColor = image.at<uchar>(x,y);
        if(pixelColor > 127) {
          map[x][y] = T(0);
        }
        else {
          map[x][y] = T(1);
        }
        std::cout << map[x][y];
      }
      std::cout << std::endl;
    }
  }

  template<typename T>
  void resistVector(const anpi::Matrix<T>&          map,
                    std::vector<T>&        resistVector) {
    resistVector.clear();//////////////////////////////////////////////////////////////////////////////////////////////////
    resistVector.resize(2 * map.rows() * map.cols() - (map.rows() + map.cols()),T(1));
    int n,m,i,j;
    for(int k = 0; k < resistVector.size(); k++) {
      anpi::inverseMapper(map.rows(),map.cols(),k,n,m,i,j);
      if(map[n][m] == 1 || map[i][j] == 1) {
        resistVector[k] = T(1000000);
      }
    }
  }

  /**
   *
   * @tparam T
   * @param rows Resistive grid rows
   * @param cols Resistive grid cols
   * @param A    Matrix to be filled
   */
  template<typename T>
  void nodos(const size_t     rows,
             const size_t     cols,
             anpi::Matrix<T>&    A) {

    if (A.rows() != (2 * rows * cols - (rows + cols)) &&
        A.cols() != (2 * rows * cols - (rows + cols))) {
      throw anpi::Exception("anpi::nodos not matching sizes");
    }

    else {
      int count = 0;
      int j;
      for (int n = 0; n < rows; n++) {
        for (int m = 0; m < cols; m++) {
          if (m - 1 >= 0) {
            anpi::mapper(rows, cols, n, m, n, m - 1, j);
            A[count][j] = T(1);
          }
          if (n - 1 >= 0) {
            anpi::mapper(rows, cols, n, m, n - 1, m, j);
            A[count][j] = T(1);
          }
          if (m + 1 < cols) {
            anpi::mapper(rows, cols, n, m, n, m + 1, j);
            A[count][j] = T(-1);
          }
          if (n + 1 < rows) {
            anpi::mapper(rows, cols, n, m, n + 1, m, j);
            A[count][j] = T(-1);
          }
          count++;
        }
      }
    }
  }

  template<typename T>
  void mallas(const size_t rows,
              const size_t cols,
              const std::vector<T> &resistVector,
              const std::vector<T> &b,
              anpi::Matrix<T> &A) {

    int size = (2 * rows * cols - (rows + cols));

    if (A.rows()            != size &&
        A.cols()            != size &&
        resistVector.size() != size &&
        b.size()            != size ) {
      throw anpi::Exception("anpi::mallas not matching sizes");
    } else {

      int x, y;
      bool first = true;

      for (int n = 0; n < rows - 1; n++) {
        for (int m = 0; m < cols - 1; m++) {
          if (first) {
            // reemplaza un ecuacion de nodo
            if (b[0] == 0) { //nodo superior izquierdo
              x = 0;
            } else if (b[cols - 1] == 0) { // nodo superior derecho
              x = cols - 1;
              mapper(rows, cols, 0, cols - 2, 0, cols - 1, y);
              A[x][y] = T(0);
              mapper(rows, cols, 0, cols - 1, 1, cols - 1, y);
              A[x][y] = T(0);
            } else { // nodo inferior derecho
              x = (rows * cols) - 1;
              mapper(rows, cols, rows - 1, cols - 1, rows - 2, cols - 1, y);
              A[x][y] = T(0);
              mapper(rows, cols, rows - 1, cols - 1, rows - 1, cols - 2, y);
              A[x][y] = T(0);
            }

          }
          mapper(rows, cols, n, m, n, m + 1, y);
          A[x][y] = resistVector[y] * T(-1);

          mapper(rows, cols, n, m + 1, n + 1, m + 1, y);
          A[x][y] = resistVector[y] * T(-1);

          mapper(rows, cols, n + 1, m + 1, n + 1, m, y);
          A[x][y] = resistVector[y] * T(1);

          mapper(rows, cols, n + 1, m, n, m, y);
          A[x][y] = resistVector[y] * T(1);

          x++;

          if (first) {
            x = rows * cols;
            first = false;
          }

        }
      }

    }

  }

  template<typename T>
  void firstMethod(const std::string mapPath,
                   const int           in,
                   const int           im,
                   const int           on,
                   const int           om,
                   std::vector<T>&   path) {

    anpi::Matrix<T> map;

    mapCreator(mapPath, map);

    std::vector<T> b;

    anpi::vectorFiller<T>(map.rows(), map.cols(), in, im, on, om, b);

    int size = 2 * map.rows() * map.cols() - (map.rows() + map.cols());

    anpi::Matrix<T> A(size,size);

    anpi::nodos<T>(map.rows(), map.cols(), A);

    std::vector<T> resistVector;

    anpi::resistVector(map,resistVector);

    anpi::mallas(map.rows(), map.cols(), resistVector, b, A);

    std::vector<T> x(size,T(1));

    anpi::solveLU<T>(A,x,b);

    path = x;
  }

}//anpi

#endif