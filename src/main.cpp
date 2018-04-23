//
// Created by ger534 on 11/04/18.
//

#include <iostream>
#include <Matrix.hpp>
#include "LUDoolittle.hpp"

int main() {

    std::cout << "corriendo el main" << std::endl;

    typedef anpi::aligned_allocator<float> aalloc;


    anpi::Matrix<float, aalloc> A = { {  1, 2, 3, 4},
                          {  5, 6, 7, 8},
                          {  9,10,11,12},
                          { 13,14,15,16 } };

    anpi::Matrix<float> ParaLU = { {  1, 2, 3, 4 },
                                   {  5, 6, 7, 8 },
                                   {  9,10,11,12 },
                                   { 13,14,15,16 } };


    anpi::Matrix<float> LU;

    std::vector<size_t> permut;

    //anpi::SIMD1::luDoolittle(ParaLU, LU, permut);

    anpi::fallback1::luDoolittle(ParaLU, LU, permut);


    anpi::Matrix<float> L;
    anpi::Matrix<float> U;


    anpi::fallback1::unpackDoolittle(LU,L,U);

    //anpi::SIMD1::unpackDoolittle(LU,L,U);

    ///METODO imprime matriz
    std::cout << "LU = [";
    for(int i =0; i < LU.rows(); i++){
        for(int j =0; j < LU.cols(); j++){
            std::cout << LU[i][j] << ", " ;
        }
        std::cout << "" << std::endl;
    }
    std::cout << "]"<< std::endl;




    ///PARA BENCHMARK
    //anpi::fallback::fill(float(2), A);
    //anpi::simd::fill(float(2), A);

    int j=0;
    int n = A.rows();
    typename sse2_traits<float>::reg_type* LUptr1  = reinterpret_cast<typename sse2_traits<float>::reg_type*>(A.data());

    ///TEMPORAL SOLO PARA HACER PRUEBAS
    for(int i = 0; i < j; i++) {
        LUptr1++;
    }

    typename sse2_traits<float>::reg_type* LUptr2;

    LUptr2 = LUptr1;

    for(int i = j+1; i<n; i++) {
        LUptr2++;
        typename sse2_traits<float>::reg_type escalar = anpi::aimpl::mm_set1<float, typename sse2_traits<float>::reg_type>(A[i][j]);
        A[i][j] = A[i][j] / A[j][j];
        const float fin = LU.cols()/(sizeof(typename sse2_traits<float>::reg_type)/ (sizeof(float)));
        for(int k = 0; k < fin; k++) {
            float temp = A[i][k];
            *LUptr2 = anpi::aimpl::mm_sub<float, typename sse2_traits<float>::reg_type>(*LUptr2, anpi::aimpl::mm_mul<float, typename sse2_traits<float>::reg_type>(escalar, *LUptr1));
            A[i][k]=temp;
        }
    }





    return 0;
}
