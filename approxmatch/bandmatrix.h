/***************************************************************************
 *   Copyright (C) 2018 Jan Fostier (jan.fostier@ugent.be)                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include <vector>

#include <cstdio>

using namespace std;

// ============================================================================
// CLASS BANDED MATRIX
// ============================================================================

// The band matrix class can best be understood as m horizontal bands each
// with a width of 2W+1 elements. It is allocated as a single m*(2W+1) array.
// For example, for W = 2 and m = 6.
// XX|XXX...|
//  X|XXXX..|
//   |XXXXX.|
//   |.XXXXX|
//   |..XXXX|X
//   |...XXX|XX
// The actual matrix is between |.| The Xs left and right of the |.| are
// allocated but should in principle not be addressed.
// The storage order is row-major.

class BandMatrix
{
private:
        std::vector<int> matrix;
        int W;

public:
        /**
         * Constructor
         * @param m Number of rows
         * @param W Number of off-diagonal elements (one sided)
         */
        BandMatrix(int m, int W) : W(W) {
                matrix.resize(m * (2*W+1));
        }

        /**
         * Operator () overloading
         * @param i Row index
         * @param j Column index
         * @return Element at position (i, j)
         */
        int operator() (int i, int j) const {
                return matrix[i * (2*W+1) + j - i + W];
        }

        /**
         * Operator () overloading
         * @param i Row index
         * @param j Column index
         * @return Reference to element at position (i, j)
         */
        int& operator() (int i, int j) {
                return matrix[i * (2*W+1) + j - i + W];
        }

        /**
         * Get the band width
         * @return The band width
         */
        int getWidth() const { return W; }
};

#endif
