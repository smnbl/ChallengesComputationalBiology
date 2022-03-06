/*
 * Own implementation of the global NW algorithm, based around Jan Fostier's implementation.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "../matrix.h"

using namespace std;

// read FASTA file
void readSequences(const string& filename, vector<string>& sequences) {
    ifstream ifs(filename.c_str());
    if (!ifs)
        throw runtime_error("error opening fasta file");

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

// costs
const int G = -3; // gap
const int M = 1; // match
const int I = -1; // mismatch

void alignGlobalNW(const string& X, const string& Y) {
    const int m = X.length();
    const int n = Y.length();

    Matrix S(m+1, n+1);

    // initialize
    for (size_t i = 0; i <= m; i++)
        S(i, 0) = i*G;
    for (size_t j = 0; j <= n; j++)
        S(0, j) = j*G;

    // fill S with scores
    for (size_t i = 1; i <= m; i++)
        for (size_t j = 1; j <= n; j++) {
            const int diag = S(i - 1, j - 1) + (X[i-1] == Y[j-1] ? M : I); // NOTE: in the course notes X[1] is the first element of X!!
            const int gapY = S(i - 1, j) + G;
            const int gapX = S(i, j - 1) + G;
            S(i, j) = max(diag, max(gapX, gapY));
        }

    string alX = "";
    string mid = "";
    string alY = "";
    size_t i = m;
    size_t j = n;

    // extract alignment
    while (i > 0 or j > 0) { // keep iterating till border
        if (i > 0 and j > 0 and S(i, j) == S(i - 1, j - 1) + (X[i-1] == Y[j-1] ? M : I)) { // diag
            alX = X[i-1] + alX;
            alY = Y[j-1] + alY;
            mid = (X[i-1] == Y[j-1] ? '|' : '*' ) + mid;
            i--; j--;
        } else if (i > 0 and S(i, j) == S(i - 1, j) + G) { // gapY
            alX = X[i-1] + alX;
            alY = '-' + alY;
            mid = ' ' + mid;
            i--;
        } else if (j > 0 and S(i, j) == S(i, j - 1) + G) { // gapX
            alX = '-' + alX;
            alY = Y[j-1] + alY;
            mid = ' ' + mid;
            j--;
        } else {
            throw runtime_error("impossible path found!");
        }
    }

    cout << "Results:" << endl;
    cout << alX << endl;
    cout << mid << endl;
    cout << alY << endl;

    cout << "Score: " << S(m, n) << endl;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr << "Usage: ./program input.fasta" << endl;
        return EXIT_FAILURE;
    }

    vector<string> sequences;
    readSequences(argv[1], sequences);

    if (sequences.size() != 2) {
                cerr << "Input FASTA file should contain only two sequences\n";
                return EXIT_FAILURE;
        }


    alignGlobalNW(sequences[0], sequences[1]);

    return EXIT_SUCCESS;
}
