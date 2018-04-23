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

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"
#include "LUCrout.hpp"
#include "LUDoolittle.hpp"


BOOST_AUTO_TEST_SUITE( Matrix )

/// Benchmark for addition operations
    template<typename T>
    class benchAdd {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _maxSize;

        /// A large matrix holding
        anpi::Matrix<T> _data;

        /// State of the benchmarked evaluation
        anpi::Matrix<T> _a;
        anpi::Matrix<T> _b;
        anpi::Matrix<T> _c;

        ///CREO EL ALLOCATOR
        typedef anpi::aligned_allocator<float> aalloc;

        anpi::Matrix<T, aalloc> _aAlloc;
        anpi::Matrix<T, aalloc> _bAlloc;

    public:
        /// Construct
        benchAdd(const size_t maxSize)
                : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

            size_t idx=0;
            for (size_t r=0;r<_maxSize;++r) {
                for (size_t c=0;c<_maxSize;++c) {
                    _data(r,c)=idx++;
                }
            }
        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size<=this->_maxSize);
            this->_a=std::move(anpi::Matrix<T>(size,size,_data.data()));
            this->_b=this->_a;
        }
    };

/// Provide the evaluation method for in-place addition
    template<typename T>
    class benchAddInPlaceFallback : public benchAdd<T> {
    public:
        /// Constructor
        benchAddInPlaceFallback(const size_t n) : benchAdd<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            std::vector<size_t> p;
            anpi::fallback1::luDoolittle(this->_a,this->_b,p);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class benchAddOnCopyFallback : public benchAdd<T> {
    public:
        /// Constructor
        benchAddOnCopyFallback(const size_t n) : benchAdd<T>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            std::vector<size_t> p;
            //anpi::luCrout(this->_a,this->_b,p);
            anpi::SIMD1::luDoolittle(this->_aAlloc,this->_bAlloc,p);

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
        anpi::Matrix<T> _b;
        anpi::Matrix<T> _c;


    public:
        /// Construct
        benchFill(const size_t maxSize)
                : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

            size_t idx=0;
            for (size_t r=0;r<_maxSize;++r) {
                for (size_t c=0;c<_maxSize;++c) {
                    _data(r,c)=idx++;
                }
            }
        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size<=this->_maxSize);
            this->_a=std::move(anpi::Matrix<T>(size,size,_data.data()));
            this->_b=this->_a;
        }
    };

/// Provide the evaluation method for in-place addition
    template<typename T>
    class benchFillInPlaceFallback : public benchFill<T> {
    public:
        /// Constructor
        benchFillInPlaceFallback(const size_t n) : benchFill<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            anpi::fallback::fill(float(2),this->_a);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class benchFillOnCopyFallback : public benchFill<T> {
    public:
        /// Constructor
        benchFillOnCopyFallback(const size_t n) : benchFill<T>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            anpi::simd::fill(float(2),this->_a);
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
            benchFillInPlaceFallback<float> baip(n);

            // Measure in place add
            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill sin SIMD",times);
            ::anpi::benchmark::plotRange(times,"Fill sin SIMD","b");
        }

        {
            benchFillOnCopyFallback<float> baip(n);

            // Measure in place add
            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("Fill con SIMD",times);
            ::anpi::benchmark::plotRange(times,"Fill con SIMD","r");
        }

        ::anpi::benchmark::show();
    }



/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE( LU ) {

        std::vector<size_t> sizes = {  24,  32,  48,  64,
                                       96, 128, 192, 256,
                                       384, 512, 768,1024,
                                       1536,2048,3072,4096};


        const size_t n=sizes.back();
        const size_t repetitions=100;
        std::vector<anpi::benchmark::measurement> times;

        {
            benchAddInPlaceFallback<float> baip(n);

            // Measure in place add
            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("LU Doolittle normal",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle normal","b");
        }

        {
            benchAddOnCopyFallback<float> baip(n);

            // Measure in place add
            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("LU Doolittle SIMD",times);
            ::anpi::benchmark::plotRange(times,"LU Doolittle SIMD","r");
        }

        ::anpi::benchmark::show();
    }
BOOST_AUTO_TEST_SUITE_END()