#include <boost/test/unit_test.hpp>
#include "Aplication.hpp"

BOOST_AUTO_TEST_SUITE( Aplication )

BOOST_AUTO_TEST_CASE( Mapper ) {

  size_t rows = 3;
  size_t cols = 3;
  int x;

  anpi::mapper(rows, cols, 1, 0, 1, 1, x);
  BOOST_CHECK(x == 5);

  anpi::mapper(rows, cols, 1, 1, 1, 0, x);
  BOOST_CHECK(x == 5);

  anpi::mapper(rows, cols, 2, 1, 2, 2, x);
  BOOST_CHECK(x == 11);

  anpi::mapper(rows, cols, 2, 2, 2, 1, x);
  BOOST_CHECK(x == 11);

  cols = 4;

  anpi::mapper(rows, cols, 1, 0, 1, 1, x);
  BOOST_CHECK(x == 7);

  anpi::mapper(rows, cols, 1, 1, 1, 0, x);
  BOOST_CHECK(x == 7);

  rows = 4;
  cols = 3;

  anpi::mapper(rows, cols, 2, 2, 2, 1, x);
  BOOST_CHECK(x == 11);

  anpi::mapper(rows, cols, 3, 1, 3, 0, x);
  BOOST_CHECK(x == 15);

  try {
    anpi::mapper(rows, cols, 3, 2, 3, 3, x);
    BOOST_CHECK_MESSAGE(false,"Mapper index out of bounds properly no properly catched");
  }
  catch(anpi::Exception& exc) {
    BOOST_CHECK_MESSAGE(true,"Mapper index out of bounds properly catched");
  }

  try {
    anpi::mapper(rows, cols, 0, 0, 0, 2, x);
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

  anpi::inverseMapper(3,3,8,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 1 && i == 2 && j == 1);

  anpi::inverseMapper(3,3,6,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 1 && i == 1 && j == 2);

  anpi::inverseMapper(3,4,9,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 2 && i == 1 && j == 3);

  anpi::inverseMapper(3,4,13,n,m,i,j);
  BOOST_CHECK(n == 1 && m == 3 && i == 2 && j == 3);

  try {
    anpi::inverseMapper(3,3,12,n,m,i,j);
    BOOST_CHECK_MESSAGE(false,"Inverse Mapper index out of bounds properly no properly catched");
  }
   catch(anpi::Exception& exc) {
    BOOST_CHECK_MESSAGE(true,"Inverese Mapper index out of bounds properly catched");
  }

}
/*
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
*/
BOOST_AUTO_TEST_CASE( VectorFiller ) {

    std::vector<float> b;
    size_t rows = 3;
    size_t cols = 3;
    anpi::vectorFiller(rows,cols, 0, 0, 2, 2, b);
    std::vector<float> v { 1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0};
    BOOST_CHECK(b == v);

    rows = 3;
    cols = 4;
    anpi::vectorFiller(rows,cols, 2, 0, 1, 2, b);
    v = { 0, 0, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    BOOST_CHECK(b == v);

    try {
      anpi::vectorFiller(rows,cols, 3, 0, 1, 2, b);
      BOOST_CHECK_MESSAGE(false,"Vector Filler index out of bounds properly no properly catched");
    }
    catch(anpi::Exception& exc) {
      BOOST_CHECK_MESSAGE(true,"Vector Filler index out of bounds properly catched");
    }

    try {
      anpi::vectorFiller(rows,cols, 2, -1, 1, 2, b);
      BOOST_CHECK_MESSAGE(false,"Vector Filler index out of bounds properly no properly catched");
    }
    catch(anpi::Exception& exc) {
      BOOST_CHECK_MESSAGE(true,"Vector Filler index out of bounds properly catched");
    }

    try {
      anpi::vectorFiller(rows,cols, 0, 0, 4, 2, b);
      BOOST_CHECK_MESSAGE(false,"Vector Filler index out of bounds properly no properly catched");
    }
    catch(anpi::Exception& exc) {
      BOOST_CHECK_MESSAGE(true,"Vector Filler index out of bounds properly catched");
    }

    try {
      anpi::vectorFiller(rows,cols, 1, 0, 1, 4, b);
      BOOST_CHECK_MESSAGE(false,"Vector Filler index out of bounds properly no properly catched");
    }
    catch(anpi::Exception& exc) {
      BOOST_CHECK_MESSAGE(true,"Vector Filler index out of bounds properly catched");
    }

}

BOOST_AUTO_TEST_CASE( MapCreator ) {
    anpi::Matrix<float> map;
    anpi::mapCreator<float>("../../images/2x2.png",map);

    anpi::Matrix<float> MAP = { {0, 1},
                                {0, 0}  };

    BOOST_CHECK(MAP == map);

    anpi::mapCreator<float>("../../images/3x3.png",map);

    MAP = { {0, 1, 0},
            {0, 1, 0},
            {0, 0, 0} };

    BOOST_CHECK(MAP == map);

}

BOOST_AUTO_TEST_CASE( ResistVector ) {
    anpi::Matrix<float> map = { { 1, 0, 1, 0},
                              { 0, 0, 0, 0},
                              { 0, 1, 0, 1},
                              { 0, 0, 1, 1} };
    std::vector<float> rvect;
    anpi::resistVector(map,rvect);
    std::vector<float> RVECT = {1000000, 1000000, 1000000,
                                1000000,       1, 1000000,       1,
                                      1,       1,       1,
                                      1, 1000000,       1, 1000000,
                                1000000, 1000000, 1000000,
                                      1, 1000000, 1000000, 1000000,
                                      1, 1000000, 1000000};
    BOOST_CHECK(rvect == RVECT);
}

BOOST_AUTO_TEST_CASE( FirstMethod ) {
  std::vector<float> x;
  anpi::firstMethod<float>("../../images/2x2.png", 0, 0, 1, 1, x);
  for(int i = 0; i < x.size(); i++) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl;

  anpi::firstMethod<float>("../../images/3x3.png", 0, 0, 0, 2, x);
  for(int i = 0; i < x.size(); i++) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl;

  anpi::firstMethod<float>("../../images/50x50.png", 1, 1, 49, 49, x);
  for(int i = 0; i < x.size(); i++) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl;

}


BOOST_AUTO_TEST_SUITE_END()