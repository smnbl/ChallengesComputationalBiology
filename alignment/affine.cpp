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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "matrix.h"

using namespace std;

/**
 * Write program usage information to the standard output
 */
void printUsage()
{
        cout << "Usage: globalNW input.fasta\n\n";
        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>" << endl;
}

/**
 * Read sequences from a FASTA file
 * @param filename FASTA file filename
 * @param sequences Vector of sequences (output)
 */
void readSequences(const string& filename, vector<string>& sequences)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        string line;
        while (ifs) {
                getline(ifs, line);
                if (line.empty())
                        continue;
                if (line.front() == '>') {
                        sequences.push_back(string());
                        continue;
                }

                sequences.back().append(line);
        }
}

/**
 * Compute the maximum of three values
 * @param a Value one
 * @param b Value two
 * @param c Value three
 * @return max(a, b, c)
 */
int max3(int a, int b, int c)
{
        return max(a, max(b, c));
}

/**
 * Compute the index containing the maximum of three values
 * @param a Value one
 * @param b Value two
 * @param c Value three
 * @return 0 when a is max, 1 when b is max, 2 when c is max
 */
int max3idx(int a, int b, int c)
{
        if (a >= b && a >= c)
                return 0;
        if (b >= a && b >= c)
                return 1;
        return 2;
}

/**
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 */
void alignAffineGap(const string& X, const string& Y)
{
        const int Go = -3;      // cost for opening a gap
        const int Ge = -1;      // cost for extending a gap
        const int M = 1;        // match
        const int I = -1;       // mismatch cost

        size_t m = X.length();
        size_t n = Y.length();

        const int minint = max(m, n) * Go;

        // initialize (m+1) x (n+1) matrices S, X and Y
        Matrix SD(m+1, n+1), SX(m+1, n+1), SY(m+1, n+1);
        SD(0, 0) = 0;
        SX(0, 0) = SY(0, 0) = minint;

        // initialize first column of S, X and Y
        for (size_t i = 1; i <= m; i++) {
                SD(i, 0) = minint;
                SX(i, 0) = minint;
                SY(i, 0) = Go + i * Ge;
        }

        // initialize first row of S, X and Y
        for (size_t j = 1; j <= n; j++) {
                SD(0, j) = minint;
                SX(0, j) = Go + j * Ge;
                SY(0, j) = minint;
        }

        // fill in the rest of the matrix
        for (size_t i = 1; i <= m; i++) {
                for (size_t j = 1; j <= n; j++) {
                        // alignment ends with two aligned characters
                        int cScore = X[i-1] == Y[j-1] ? M : I;
                        SD(i, j) = max3(SD(i-1, j-1) + cScore,
                                        SX(i-1, j-1) + cScore,
                                        SY(i-1, j-1) + cScore);

                        // alignment ends with a gap in X
                        SX(i, j) = max3(SD(i, j-1) + Go + Ge,
                                        SX(i, j-1) + Ge,
                                        SY(i, j-1) + Go + Ge);

                        // alignment ends with a gap in Y
                        SY(i, j) = max3(SD(i-1, j) + Go + Ge,
                                        SX(i-1, j) + Go + Ge,
                                        SY(i-1, j) + Ge);
                }
        }

        // create an alignment
        string alX, alY, mid;

        int i = (int)X.size();
        int j = (int)Y.size();
        int state = max3idx(SD(i, j), SX(i, j), SY(i, j));

        while (i > 0 || j > 0) {
                if (state == 0) { // diag
                        alX.push_back(X[i-1]);
                        alY.push_back(Y[j-1]);
                        mid.push_back((X[i-1] == Y[j-1]) ? '|' : '*');
                        int cScore = X[i-1] == Y[j-1] ? M : I;
                        state = max3idx(SD(i-1, j-1) + cScore,
                                        SX(i-1, j-1) + cScore,
                                        SY(i-1, j-1) + cScore);
                        i--;
                        j--;
                } else if (state == 1) { // gap in X
                        alX.push_back('-');
                        alY.push_back(Y[j-1]);
                        mid.push_back(' ');
                        state = max3idx(SD(i, j-1) + Go + Ge,
                                        SX(i, j-1) + Ge,
                                        SY(i, j-1) + Go + Ge);
                        j--;
                } else { // gap in Y
                        alX.push_back(X[i-1]);
                        alY.push_back('-');
                        mid.push_back(' ');
                        state = max3idx(SD(i-1, j) + Go + Ge,
                                        SX(i-1, j) + Go + Ge,
                                        SY(i-1, j) + Ge);
                        i--;
                }
        }

        reverse(alX.begin(), alX.end());
        reverse(alY.begin(), alY.end());
        reverse(mid.begin(), mid.end());

        for (size_t i = 0; i < alX.size(); i += 80) {
                cout << alX.substr(i, 80) << "\n"
                     << mid.substr(i, 80) << "\n"
                     << alY.substr(i, 80) << "\n\n";
        }
        cout << "Alignment score: " << max3(SD(m, n), SX(m, n), SY(m, n)) << endl;
}

int main(int argc, char** argv)
{
        if (argc != 2) {
                printUsage();
                return EXIT_FAILURE;
        }

        vector<string> sequences;
        readSequences(argv[1], sequences);

        if (sequences.size() != 2) {
                cerr << "Input FASTA file should contain only two sequences\n";
                return EXIT_FAILURE;
        }

        alignAffineGap(sequences[0], sequences[1]);

        return EXIT_SUCCESS;
}
