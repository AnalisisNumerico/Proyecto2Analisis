//
// Created by ger534 on 22/04/18.
//

/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
//#include "bits/MatrixArithmetic.hpp"

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"
#include "LUDoolittle.hpp"


BOOST_AUTO_TEST_SUITE( Matrix )

/// Benchmark for addition operations
    template<typename T>
    class benchAdd {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _size;

        /// A large matrix holding
        anpi::Matrix<T> A;

        anpi::Matrix<T> L;
        anpi::Matrix<T> U;

        anpi::Matrix<T> LU;

    public:
        /// Construct
        benchAdd(const size_t size)
                : _size(size),A(size,size,anpi::DoNotInitialize) {
        }


        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            U.allocate(size,size);
            L.allocate(size,size);

            for(int j = 0; j < size; j++) {
                for(int i = 0; i <= j; i++) {
                    U[i][j] = std::rand() % (( 101 ) );
                }
            }

            for(int i = 1; i < size; i++) {
                for(int j = 0; j < i; j++) {
                    L[i][j] = std::rand() % (( 101 ) );
                }
            }

            for(int i = 0; i < size; i++) { //diagonal
                L[i][i] = T(1);
            }


            A = L*U;

            LU.allocate(size,size);

            std::cout << size << std::endl;
        }
    };

    template<typename T,class Alloc>
    class benchDooLittleSIMD {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _size;

        /// A large matrix holding
        anpi::Matrix<T, Alloc> A;

        anpi::Matrix<T, Alloc> L;
        anpi::Matrix<T, Alloc> U;

        anpi::Matrix<T, Alloc> LU;

    public:
        /// Construct
        benchDooLittleSIMD(const size_t size) : _size(size), A(size,size,anpi::DoNotInitialize) {
        }



        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
         //assert (size<=this->_size);

            U.allocate(size,size);
            L.allocate(size,size);

            for(int j = 0; j < size; j++) {
                for(int i = 0; i <= j; i++) {
                    U[i][j] = std::rand() % (( 101 ) );
                }
            }

            for(int i = 1; i < size; i++) {
                for(int j = 0; j < i; j++) {
                    L[i][j] = std::rand() % (( 101 ) );
                }
            }

            for(int i = 0; i < size; i++) { //diagonal
                L[i][i] = T(1);
            }


            A = L*U;

            LU.allocate(size,size);

            std::cout << size << std::endl;

        }
    };


/// Provide the evaluation method for in-place addition
    template<typename T>
    class luDoolittleFallback : public benchAdd<T> {
    public:
        /// Constructor
        luDoolittleFallback(const size_t n) : benchAdd<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            std::vector<size_t> p;
            anpi::fallback1::luDoolittle(this->A,this->LU,p);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T,class Alloc>
    class luDoolittleSIMD : public benchDooLittleSIMD<T,Alloc> {
    public:
        /// Constructor
        luDoolittleSIMD(const size_t n) : benchDooLittleSIMD<T, Alloc>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            std::vector<size_t> p;
            anpi::SIMD1::luDoolittle(this->A,this->LU,p);

        }
    };




    /// Benchmark for addition operations
    template<typename T>
    class benchFill {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _maxSize;

        /// A large matrix holding
        anpi::Matrix<T> _data;

        /// State of the benchmarked evaluation
        anpi::Matrix<T> _a;

        typedef anpi::aligned_allocator<T> alloc;
        anpi::Matrix<T,alloc> _aAlloc;

    public:
        /// Construct
        benchFill(const size_t maxSize)
                : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size<=this->_maxSize);
            this->_a.allocate(size,size);
            this->_aAlloc.allocate(size,size);
        }
    };

/// Provide the evaluation method for in-place addition
    template<typename T>
    class benchFillFallback : public benchFill<T> {
    public:
        /// Constructor
        benchFillFallback(const size_t n) : benchFill<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            anpi::fallback::fill(T(2),this->_aAlloc);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class benchFillSIMD : public benchFill<T> {
    public:
        /// Constructor
        benchFillSIMD(const size_t n) : benchFill<T>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            anpi::simd::fill(T(2),this->_aAlloc);
        }
    };


    BOOST_AUTO_TEST_CASE( FILL ) {

        std::vector<size_t> sizes = {  24,  32,  48,  64,
                                       96, 128, 192, 256,
                                       384, 512, 768,1024,
                                       1536,2048,3072,4096};

        const size_t n=sizes.back();
        const size_t repetitions=100;
        std::vector<anpi::benchmark::measurement> times;

        {
            benchFillFallback<float> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill sin SIMD presicion simple",times);
            ::anpi::benchmark::plotRange(times,"Fill sin SIMD presicion simple","b");
        }

        {
            benchFillFallback<double> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill sin SIMD presicion doble",times);
            ::anpi::benchmark::plotRange(times,"Fill sin SIMD presicion doble","r");
        }

        {
            benchFillSIMD<float> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill con SIMD presicion simple",times);
            ::anpi::benchmark::plotRange(times,"Fill con SIMD presicion simple","y");
        }

        {
            benchFillSIMD<double> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill con SIMD presicion doble",times);
            ::anpi::benchmark::plotRange(times,"Fill con SIMD presicion doble","g");
        }

        ::anpi::benchmark::show();
    }



/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE( LU ) {

        std::vector<size_t> sizes = {  24,  32,  48,  64,
                                       96, 128, 192, 256
                                       };

        size_t n=sizes.front();
        const size_t repetitions=100;
        std::vector<anpi::benchmark::measurement> times;

        {
            luDoolittleFallback<float> baip(n);


            ANPI_BENCHMARK(sizes,repetitions,times,baip);


            ::anpi::benchmark::write("LU Doolittle normal precision simple",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle normal precision simple","b");
        }

        {
            luDoolittleFallback<double> baip(n);


            ANPI_BENCHMARK(sizes,repetitions,times,baip);


            ::anpi::benchmark::write("LU Doolittle normal precision doble",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle normal  precision doble","r");
        }

        {
            ///CREO EL ALLOCATOR
            typedef anpi::aligned_allocator<float> aalloc;
            luDoolittleSIMD<float, aalloc> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("LU Doolittle SIMD precision simple",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle SIMD precision simple","y");
        }

        {
            ///CREO EL ALLOCATOR
            typedef anpi::aligned_allocator<double> aalloc;
            luDoolittleSIMD<double, aalloc> baip(n);

            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("LU Doolittle SIMD precision doble",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle SIMD precision doble","g");
        }

        ::anpi::benchmark::show();
    }
BOOST_AUTO_TEST_SUITE_END()