/**
 * MIT License
 *
 * Copyright (c) 2021 Panagiotis Repouskos
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * ---
 *
 * The code to compute Hilbert values was ported in c++ from python code
 * found in https://github.com/galtay/hilbertcurve:
 *
 * MIT License
 *
 * Copyright (c) 2017 Gabriel Altay
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */







#ifndef HILBERT_CURVE_HILBERTCURVE_H
#define HILBERT_CURVE_HILBERTCURVE_H

#include <bitset>
#include <iostream>
#include <algorithm>
#include <vector>

namespace HilbertCurve {

    /**
     * A class that implements the Hilbert curve.
     */
    class HilbertCurve {
    public:

        /**
         * Initialize an instance of this class
         * @param[in] dimension The dimension of the integer grid lattice
         * @param[in] iterations The iterations to use in constructing the curve (recursion depth)
         */
        HilbertCurve(unsigned int dimension, unsigned int iterations);

        /**
         * Compute the point in the hypercube lattice that has the 'hilbertNumber' value.
         * @param[in] hilbertNumber The Hilbert value of the wanted point
         * @return A vector with the coordinates of the point with 'hilbertNumber' value.
         */
        std::vector<unsigned long> pointFromHilbertNumber(int hilbertNumber);

        /**
         * Compute the Hilbert value of a point in the lattice hypercube.
         * @param[in] point A point in the hypercube lattice
         * @return The Hilbert value of 'point'
         */
        unsigned long hilbertNumberFromPoint(std::vector<unsigned long> const &point);

        /**
         * Maps the data to points in the hypercube lattice and sorts w.r.t. their Hilbert value.
         * @tparam DATA_TYPE The type of the data
         * @tparam COORDS_FUNCTOR A functor that takes a DATA_TYPE, maps it to a point in the hypercube and returns a std::vector<unsigned long> with its coordinates
         * @param[in, out] data A vector with all the data
         * @param[out] hilbertValues The Hilbert values of the sorted data
         */
        template<typename DATA_TYPE, typename COORDS_FUNCTOR>
        void sortData(std::vector<DATA_TYPE *> &data, std::vector<unsigned long> &hilbertValues);

    private:

        /// The dimension of the integer grid lattice
        int dimension;

        /// The iterations to use in constructing the curve (recursion depth)
        int iterations;

        void hilbertIntegerToTranspose(int hilbertInteger, std::vector<unsigned long> &transpose);

        unsigned long transposeToHilbertInteger(std::vector<unsigned long> const &x);
    };


    HilbertCurve::HilbertCurve(unsigned int dimension, unsigned int iterations) :
            dimension(dimension), iterations(iterations) {
    }


    std::vector<unsigned long> HilbertCurve::pointFromHilbertNumber(int hilbertNumber) {
        std::vector<unsigned long> x(dimension);

        hilbertIntegerToTranspose(hilbertNumber, x);
        long z = 2 << (iterations - 1);

        long t = x[dimension - 1] >> 1;
        for (int i = dimension - 1; i > 0; i--)
            x[i] = x[i] ^ x[i - 1];
        x[0] = x[0] ^ t;

        long q = 2;

        while (q != z) {
            long p = q - 1;

            for (int i = dimension - 1; i > -1; --i) {
                if (x[i] & q)
                    x[0] = x[0] ^ p;
                else {
                    t = (x[0] ^ x[i]) & p;
                    x[0] = x[0] ^ t;
                    x[i] = x[i] ^ t;
                }
            }

            q = q << 1;
        }

        return x;
    }

    void HilbertCurve::hilbertIntegerToTranspose(int hilbertInteger, std::vector<unsigned long> &transpose) {
        std::bitset<sizeof(unsigned long) * 8> bits(hilbertInteger);
        int size = iterations * dimension;

        for (int i = 0; i < dimension; i++) {

            int coord = 0;
            int rep = iterations - 1;

            for (int j = i; j < size; j += dimension) {
                if (bits[size - 1 - j])
                    coord += static_cast<int>(1 << rep);
                rep--;
            }

            transpose[i] = coord;
        }
    }


    unsigned long HilbertCurve::transposeToHilbertInteger(const std::vector<unsigned long> &x) {
        typedef std::bitset<sizeof(unsigned long) * 8> BitSet;
        std::vector<BitSet> xBitString(dimension);

        for (int i = 0; i < dimension; ++i)
            xBitString[i] = BitSet(x[i]);

        int size = dimension * iterations;

        int at = 0;
        unsigned long hilbertInteger = 0;

        for (int i = 0; i < iterations; ++i)
            for (int j = 0; j < dimension; ++j) {
                if (xBitString[j][iterations - 1 - i])
                    hilbertInteger += static_cast<unsigned long >(1 << (size - 1 - at));
                at++;
            }

        return hilbertInteger;
    }

    unsigned long HilbertCurve::hilbertNumberFromPoint(const std::vector<unsigned long> &point) {
        std::vector<unsigned long> _point(point);

        unsigned long m = 1 << (iterations - 1);
        unsigned long q = m;

        while (q > 1) {
            unsigned long p = q - 1;

            for (int i = 0; i < dimension; ++i) {
                if (_point[i] & q)
                    _point[0] = _point[0] ^ p;
                else {
                    unsigned long t = (_point[0] ^ _point[i]) & p;
                    _point[0] = _point[0] ^ t;
                    _point[i] = _point[i] ^ t;
                }
            }

            q = q >> 1;
        }

        for (int i = 1; i < dimension; ++i)
            _point[i] = _point[i] ^ _point[i - 1];
        unsigned long t = 0;
        q = m;

        while (q > 1) {
            if (_point[dimension - 1] & q)
                t = t ^ (q - 1);

            q = q >> 1;
        }

        for (int i = 0; i < dimension; ++i)
            _point[i] = _point[i] ^ t;

        return transposeToHilbertInteger(_point);
    }

    template<typename DATA_TYPE, typename COORDS_FUNCTOR>
    void HilbertCurve::sortData(std::vector<DATA_TYPE *> &data, std::vector<unsigned long> &hilbertValues) {
        // we need this class for std::sort()
        class DATA {
        public:
            DATA_TYPE *d;
            unsigned long hilbertValue;

            DATA(DATA_TYPE *d, unsigned long hilbertValue) : d(d), hilbertValue(hilbertValue) {};

            DATA() = default;
        };

        std::vector<DATA> _data(data.size());
        COORDS_FUNCTOR functor;

        // compute the Hilbert value of all data
        int at = 0;
        for (auto &d : data) {
            _data[at] = DATA(d, hilbertNumberFromPoint(functor(*d)));
            ++at;
        }

        // sort
        struct comp {
            bool operator()(DATA &d1, DATA &d2) { return d1.hilbertValue < d2.hilbertValue; };
        } comp;

        std::sort(_data.begin(), _data.end(), comp);
        hilbertValues = std::vector<unsigned long>(data.size());
        std::vector<DATA_TYPE *> retData(data.size());

        at = 0;
        for (DATA &d : _data) {
            retData[at] = d.d;
            hilbertValues[at] = d.hilbertValue;
            ++at;
        }

        data = retData;
    }


}

#endif //HILBERT_CURVE_HILBERTCURVE_H
