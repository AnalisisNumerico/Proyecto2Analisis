#ifndef ANPI_APLICACION_HPP
#define ANPI_APLICACION_HPP

#include "Matrix.hpp"
#include <opencv2/opencv.hpp>

namespace anpi {

  template<typename T>
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

  template<typename T>
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
      std::vector<T> v(2*rows*cols);
      v[cols * n + m] = T( 1);
      v[cols * i + j] = T(-1);
      b = v;
    }
  }

  template<typename T>
  void mapCreator(anpi::Matrix<T>& A) {
    cv::Mat image = cv::imread("/home/jeanpaul/Downloads/45.png",0);
    anpi::Matrix<T> a(image.rows,image.cols);
    int pixelColor;
    for (int x = 0;x < image.rows; x++) {
      for (int y = 0; y < image.cols; y++) {
        pixelColor = image.at<uchar>(x,y);
        if(pixelColor < 127) {
          A[x][y] = T(0);
        }
        else {
          A[x][y] = T(1);
        }
      }
      std::cout << std::endl;
    }
    A = a;
  }

  template<typename T>
  void resistVector(const anpi::Matrix<T>&          map,
                    std::vector<T>&        resistVector) {
    std::vector<T> rVec(2 * map.rows() * map.cols() - (map.rows() + map.cols()),T(1));
    int n,m,i,j;
    for(int k = 0; k < rVec.size(); k++) {
      anpi::inverseMapper<T>(map.rows(),map.cols(),k,n,m,i,j);
      if(map[n][m] == 1 || map[i][j] == 1) {
        rVec[k] = T(1000000);
      }
    }
    resistVector = rVec;
  }

}//anpi

#endif