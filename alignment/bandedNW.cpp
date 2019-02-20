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
        cout << "Usage: bandedNW input.fasta\n\n";
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
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 * @param W number of diagonals on each side of the main diagonal (W >= 0)
 */
void alignBandedNW(const string& X, const string& Y, int W)
{
        const int G = -3;
        const int M = 1;
        const int I = -1;

        int m = (int)X.length();
        int n = (int)Y.length();

        // check the dimensions of s1 and s2
        if ((max(m, n) - min(m, n)) > W) {
                cerr << "Cannot align sequences with length " << m << " and "
                    << n << " as the maximum number of gaps is " << W << "\n";
                exit(EXIT_FAILURE);
        }

        // no need for more than n extra diagonals
        if (W > n)
                W = n;

        BandMatrix S(m+1, W);

        // initialize first column
        for (int i = 0; i <= min(W, m); i++)
                S(i, 0) = i * G;

        // initialize first row
        for (int j = 1; j <= W; j++)
                S(0, j) = j * G;

        // fill in the rest of the matrix
        for (int i = 1; i <= m; i++) {
                int startj = max(1, i - W);
                int endj = min(n, i + W);

                for (int j = startj; j <= endj; j++) {
                        int diag = S(i-1, j-1) + (X[i-1] == Y[j-1] ? M : I);
                        int gapX = (j > i - W) ? S(i, j-1) + G : diag - 1;
                        int gapY = (j < i + W) ? S(i-1, j) + G : diag - 1;
                        S(i, j) = max(max(diag, gapX), gapY);
                }
        }

        // create an alignment
        string alX, alY, mid;

        int i = (int)X.size();
        int j = (int)Y.size();

        while (i > 0 || j > 0) {
                if ((i > 0) && (j < i + W) && (S(i, j) == S(i-1, j) + G)) {
                        alX.push_back(X[i-1]);
                        alY.push_back('-');
                        mid.push_back(' ');
                        i--;
                } else if ((j > 0) && (j > i - W) && (S(i, j) == S(i, j-1) + G)) {
                        alX.push_back('-');
                        alY.push_back(Y[j-1]);
                        mid.push_back(' ');
                        j--;
                } else {
                        alX.push_back(X[i-1]);
                        alY.push_back(Y[j-1]);
                        char c = (X[i-1] == Y[j-1]) ? '|' : '*';
                        mid.push_back(c);
                        i--;
                        j--;
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
        cout << "Alignment score: " << S(m, n) << endl;
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

        alignBandedNW(sequences[0], sequences[1], 3);

        return EXIT_SUCCESS;
}
