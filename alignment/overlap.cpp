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
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 */
void alignOverlap(const string& X, const string& Y)
{
        const int G = -3;
        const int M = 1;
        const int I = -1;

        int m = (int)X.length();
        int n = (int)Y.length();

        // initialize an (m+1) x (n+1) matrix S
        Matrix S(m+1, n+1);

        // initialize first column
        for (int i = 0; i <= m; i++)
                S(i, 0) = 0;

        // initialize first row
        for (int j = 1; j <= n; j++)
                S(0, j) = 0;

        // fill in the rest of the matrix
        for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                        int diag = S(i-1, j-1) + (X[i-1] == Y[j-1] ? M : I);
                        int gapX = S(i, j-1) + G;
                        int gapY = S(i-1, j) + G;
                        S(i, j) = max(max(diag, gapX), gapY);
                }
        }

        // find the endpoint of the alignment
        int maxVal = 0, maxi = 0, maxj = 0;
        for (int i = 0; i <= m; i++) {
                if (S(i, n) > maxVal) {
                        maxVal = S(i, n);
                        maxi = i;
                        maxj = n;
                }
        }

        for (int j = 0; j <= n; j++) {
                if (S(m, j) > maxVal) {
                        maxVal = S(m, j);
                        maxi = m;
                        maxj = j;
                }
        }

        // create an alignment
        string alX, alY, mid;

        int i = maxi;
        int j = maxj;

        while (i > 0 && j > 0) {
                if ((i > 0) && (S(i, j) == S(i-1, j) + G)) {
                        alX.push_back(X[i-1]);
                        alY.push_back('-');
                        mid.push_back(' ');
                        i--;
                } else if ((j > 0) && (S(i, j) == S(i, j-1) + G)) {
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

        for (int c = 0; c < i; c++) {
                alY.push_back(' ');
                mid.push_back(' ');
        }

        for (int c = 0; c < j; c++) {
                alX.push_back(' ');
                mid.push_back(' ');
        }

        reverse(alX.begin(), alX.end());
        reverse(alY.begin(), alY.end());
        reverse(mid.begin(), mid.end());

        alX = X.substr(0, i) + alX + X.substr(maxi);
        alY = Y.substr(0, j) + alY + Y.substr(maxj);

        cout << "Overlap alignment X[" << i << "-" << maxi-1 << "], Y["
             << j << "-" << maxj-1 << "]\n";
        for (size_t i = 0; i < alX.size(); i += 80) {
                cout << alX.substr(i, 80) << "\n"
                     << mid.substr(i, 80) << "\n"
                     << alY.substr(i, 80) << "\n\n";
        }
        cout << "Alignment score: " << S(maxi, maxj) << endl;
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

        alignOverlap(sequences[0], sequences[1]);

        return EXIT_SUCCESS;
}
