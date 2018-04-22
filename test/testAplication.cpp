#include <boost/test/unit_test.hpp>
#include "Aplication.hpp"

BOOST_AUTO_TEST_SUITE( Aplication )

BOOST_AUTO_TEST_CASE( Mapper ) {

  size_t rows = 3;
  size_t cols = 3;
  int x;

  anpi::mapper<float>(rows, cols, 1, 0, 1, 1, x);
  BOOST_CHECK(x == 5);

  anpi::mapper<float>(rows, cols, 1, 1, 1, 0, x);
  BOOST_CHECK(x == 5);

  anpi::mapper<float>(rows, cols, 2, 1, 2, 2, x);
  BOOST_CHECK(x == 11);

  anpi::mapper<float>(rows, cols, 2, 2, 2, 1, x);
  BOOST_CHECK(x == 11);

  cols = 4;

  anpi::mapper<float>(rows, cols, 1, 0, 1, 1, x);
  BOOST_CHECK(x == 7);

  anpi::mapper<float>(rows, cols, 1, 1, 1, 0, x);
  BOOST_CHECK(x == 7);

  rows = 4;
  cols = 3;

  anpi::mapper<float>(rows, cols, 2, 2, 2, 1, x);
  BOOST_CHECK(x == 11);

  anpi::mapper<float>(rows, cols, 3, 1, 3, 0, x);
  BOOST_CHECK(x == 15);

  try {
    anpi::mapper<float>(rows, cols, 3, 2, 3, 3, x);
    BOOST_CHECK_MESSAGE(false,"Mapper index out of bounds properly no properly catched");
  }
  catch(anpi::Exception& exc) {
    BOOST_CHECK_MESSAGE(true,"Mapper index out of bounds properly catched");
  }

  try {
    anpi::mapper<float>(rows, cols, 0, 0, 0, 2, x);
    BOOST_CHECK_MESSAGE(false,"Mapper invalid nodes bounds not properly catched");
  }
  catch(anpi::Exception& exc) {
    BOOST_CHECK_MESSAGE(true,"Mapper invalid nodes bounds properly detected");
  }

}

BOOST_AUTO_TEST_CASE( InverseMapper ) {

  int n;
  int m;
  int i;
  int j;

  anpi::inverseMapper<float>(3,3,8,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 1 && i == 2 && j == 1);

  anpi::inverseMapper<float>(3,3,6,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 1 && i == 1 && j == 2);

  anpi::inverseMapper<float>(3,4,9,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 2 && i == 1 && j == 3);

  anpi::inverseMapper<float>(3,4,13,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 3 && i == 2 && j == 3);

  try {
    anpi::inverseMapper<float>(3,3,12,n,m,i,j);
    BOOST_CHECK_MESSAGE(false,"Inverse Mapper index out of bounds properly no properly catched");
  }
   catch(anpi::Exception& exc) {
    BOOST_CHECK_MESSAGE(true,"Inverese Mapper index out of bounds properly catched");
  }

}

BOOST_AUTO_TEST_CASE( MatrixFiller ) {

  anpi::Matrix<float> a;
  size_t rows = 2;
  size_t cols = 2;
  anpi::matrixFiller<float>(rows, cols, a);
  anpi::Matrix<float> A = { {-1,-1, 0, 0},
                            { 1, 0,-1, 0},
                            { 0, 1, 0,-1},
                            { 0, 0, 1, 1} };
  BOOST_CHECK(a == A);

  rows = 3;
  cols = 3;
  anpi::matrixFiller<float>(rows, cols, a);

  A = { {-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 1,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 1, 0, 0,-1, 0,-1, 0, 0, 0, 0},
        { 0, 0, 0, 1, 0, 1,-1, 0,-1, 0, 0, 0},
        { 0, 0, 0, 0, 1, 0, 1, 0, 0,-1, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,-1},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}};

  BOOST_CHECK(a == A);

}

BOOST_AUTO_TEST_CASE( VectorFiller ) {

    std::vector<float> b;
    size_t rows = 3;
    size_t cols = 3;
    anpi::vectorFiller(rows,cols, 0, 0, 2, 2, b);
    std::vector<float> v {1,0,0,0,0,0,0,0,-1};
    BOOST_CHECK(b == v);
}

BOOST_AUTO_TEST_SUITE_END()