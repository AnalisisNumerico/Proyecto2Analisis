//
// Created by ger534 on 11/04/18.
//

#include <iostream>
#include <Matrix.hpp>
#include <Intrinsics.hpp>
#include <bits/MatrixArithmetic.hpp>
#include "LUDoolittle.hpp"
#include "PlotPyTNSHA.hpp"

int main() {

    typedef anpi::aligned_allocator<float> aalloc;

    anpi::Matrix<float, aalloc> B = { {  1, 2, 3, 4},
                                      {  5, 6, 7, 8},
                                      {  9,10,11,12},
                                      { 13,14,15,16 },
                                      { 17,18,19,20 } };

    anpi::Matrix<float, aalloc> b(B.rows(),B.cols());

    typename sse2_traits<float>::reg_type* ptr1  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(B.data());

    typename sse2_traits<float>::reg_type* ptr2  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(b.data());


    typename sse2_traits<float>::reg_type escalar = anpi::aimpl::mm_set1<float, typename sse2_traits<float>::reg_type>(float(0));

    const float fin = B.cols()/(sizeof(typename sse2_traits<float>::reg_type)/ (sizeof(float)));
    float * vector_ptr = (float*)&escalar;




    for(int j = 0; j < B.rows(); j++) {

        for(int i = 0; i <= fin; i++) {
           // *ptr2 = anpi::aimpl::mm_add<float, typename sse2_traits<float>::reg_type>(escalar, *ptr1);
           // *vector_ptr =
        }
        ptr2++;
        ptr1++;
    }

    ///METODO imprime matriz
    std::cout << "b = [";
    for(int i =0; i < b.rows(); i++){
        for(int j =0; j < b.cols(); j++){
            std::cout << b[i][j] << ", " ;
        }
        std::cout << "" << std::endl;
    }
    std::cout << "]"<< std::endl;


    int size = 4;
    anpi::Matrix<float> L1(size,size);

    for(int j = 0; j < size; j++) {
        for(int i = 0; i <= j; i++) {
            L1[i][j] = std::rand() % (( 101 ) );
        }
    }

    ///METODO imprime matriz
    std::cout << "L1 = [";
    for(int i =0; i < L1.rows(); i++){
        for(int j =0; j < L1.cols(); j++){
            std::cout <<L1[i][j] << ", " ;
        }
        std::cout << "" << std::endl;
    }
    std::cout << "]"<< std::endl;


    anpi::Matrix<float> a1 = { {1,2,3},{ 6, 5, 4} };
    std::vector<float> b1 = { 1,2,-3 };
    std::vector<float> c1 (a1.rows());

    for(size_t i=0;i<a1.rows();++i){
        float sum =0;
        for(size_t k=0;k<b1.size();++k){
            //sum = sum + a1(i,k)*b1[k];
            sum = a1(i,k)*b1[k];
        }
        c1[i]=sum;
    }

    for(int i = 0; i<c1.size();i++){
        std::cout<<c1[i]<< std::endl;
    }
    std::cout<<"++++++++"<< std::endl;

    typename sse2_traits<float>::reg_type* aptr  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(a1.data());
    typename sse2_traits<float>::reg_type* bptr  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(b1.data());
    typename sse2_traits<float>::reg_type* cptr  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(c1.data());

//    for(size_t i=0;i<a1.rows();++i){
    for(size_t i=0;i<1;++i){

        //T sum=T();
        typename sse2_traits<float>::reg_type escalar = anpi::aimpl::mm_set1<float, typename sse2_traits<float>::reg_type>(b1[i]);
        for(size_t k=0;k<fin;k++){
            //sum = sum + a(i,k)*b[k];
            //*cptr = anpi::aimpl::mm_add<float, typename sse2_traits<float>::reg_type>(*cptr, anpi::aimpl::mm_mul<float, typename sse2_traits<float>::reg_type>(*aptr, *bptr));
            //*cptr = anpi::aimpl::mm_add<float, typename sse2_traits<float>::reg_type>(*cptr, anpi::aimpl::mm_mul<float, typename sse2_traits<float>::reg_type>(*aptr++, *bptr++));
            *cptr = anpi::aimpl::mm_mul<float, typename sse2_traits<float>::reg_type>(*aptr, *bptr);
            //*cptr = anpi::aimpl::mm_add<float, typename sse2_traits<float>::reg_type>(anpi::aimpl::mm_mul<float, typename sse2_traits<float>::reg_type>(*aptr++, *bptr++),*cptr);

        }
        //*cptr++;
        //c[i]=sum;
    }

    for(int i = 0; i<c1.size();i++){
        std::cout<<c1[i]<< std::endl;
    }

    anpi::PlotTNSHA<double> plotter;
    plotter.initialize();

    std::vector<double> bb = {0,2,2,1,1,8};
    std::vector<double> cc = {0,0,4,4,6,6};

    std::vector<double> uu = {0,2,2,1,1,8,0,0,4,4,6,6};
    std::vector<double> vv = {0,0,4,4,6,6,0,2,2,1,1,8};

    plotter.plot(bb,cc,uu,vv);
    plotter.show();

/*
 b = [1, 2, 3, 4,
    0, 6, 7, 8,
    0, 0, 11, 12,
    0, 0, 0, 16]
 */

    return 0;
}
