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
#include <vector>
#include <algorithm>
#include <tuple>

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
 * Reverse a string
 * @param s input string
 * @return reversed string
 */
string reverse(const string& s)
{
        string rev = s;
        reverse(rev.begin(), rev.end());
        return rev;
}

/**
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 * @return Last line of the alignment matrix S
 */
vector<int> NWScore(const string& X, const string& Y)
{
        const int G = -3;
        const int M = 1;
        const int I = -1;

        size_t m = X.length();
        size_t n = Y.length();

        vector<int> cr(n+1), pr(n+1);   // cr = current row, pr = previous row

        // initialize first row
        for (size_t j = 0; j <= n; j++)
                pr[j] = j * G;

        // fill in the rest of the matrix
        for (size_t i = 1; i <= m; i++) {
                cr[0] = i * G;
                for (size_t j = 1; j <= n; j++) {
                        int diag = pr[j-1] + (X[i-1] == Y[j-1] ? M : I);
                        int gapX = cr[j-1] + G;
                        int gapY = pr[j] + G;
                        cr[j] = max(max(diag, gapX), gapY);
                }

                pr = cr;
        }

        return cr;
}

/**
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 */
tuple<string, string, string, int> hirschberg(const string& X, const string& Y)
{
        const int G = -3;
        const int M = 1;
        const int I = -1;

        // create an alignment
        string alX, alY, mid;
        int score = 0;

        if (X.empty()) { // trivial case A: align Y to an empty string
                for (size_t i = 0; i < Y.length(); i++) {
                        alX.push_back('-');
                        mid.push_back(' ');
                        alY.push_back(Y[i]);
                        score += G;
                }
        } else if (Y.empty()) { // trivial case B: align X to an empty string
                for (size_t i = 0; i < X.length(); i++) {
                        alX.push_back(X[i]);
                        mid.push_back(' ');
                        alY.push_back('-');
                        score += G;
                }
        } else if (X.length() == 1 && Y.length() == 1) { // trivial case C: align two strings of length 1
                if (I < 2*G) {
                        alX.push_back(X[0]);
                        mid.push_back(' ');
                        alY.push_back('-');
                        alX.push_back('-');
                        mid.push_back(' ');
                        alY.push_back(Y[0]);
                        score = 2*G;
                } else {
                        alX.push_back(X[0]);
                        mid.push_back(X[0] == Y[0] ? '|' : '*');
                        alY.push_back(Y[0]);
                        score = (X[0] == Y[0]) ? M : I;
                }
        } else {
                size_t xmid = X.length() / 2;

                // align the left half of X
                string Xl = X.substr(0, xmid);
                vector<int> scoreL = NWScore(Xl, Y);

                // align the right half of Y
                string Xr = X.substr(xmid);
                vector<int> scoreR = NWScore(reverse(Xr), reverse(Y));
                reverse(scoreR.begin(), scoreR.end());

                // sum the score halves, store in scoreL
                transform(scoreL.begin(), scoreL.end(), scoreR.begin(), scoreL.begin(), plus<int>());

                // find the position with the highest score
                auto maxIt = max_element(scoreL.begin(), scoreL.end());
                int ymid = distance(scoreL.begin(), maxIt);
                string Yl = Y.substr(0, ymid);
                string Yr = Y.substr(ymid);

                tuple<string, string, string, int> left = hirschberg(Xl, Yl);
                tuple<string, string, string, int> right = hirschberg(Xr, Yr);

                alX = get<0>(left).append(get<0>(right));
                mid = get<1>(left).append(get<1>(right));
                alY = get<2>(left).append(get<2>(right));
                score = get<3>(left) + get<3>(right);
        }

        return make_tuple(alX, mid, alY, score);
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

        tuple<string, string, string, int> ret = hirschberg(sequences[0], sequences[1]);

        for (size_t i = 0; i < get<0>(ret).size(); i += 80) {
                cout << get<0>(ret).substr(i, 80) << "\n"
                     << get<1>(ret).substr(i, 80) << "\n"
                     << get<2>(ret).substr(i, 80) << "\n\n";
        }
        cout << "Alignment score: " << get<3>(ret) << endl;

        return EXIT_SUCCESS;
}
