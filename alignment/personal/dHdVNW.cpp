/*
 * NW with ΔH and ΔV matrices, based around Jan Fostier's implementation.
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

void fillS(const string& X, const string& Y, Matrix& S) {
    const int m = X.length();
    const int n = Y.length();

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
}

void filldHdV(const string& X, const string& Y, Matrix& dV, Matrix& dH) {
    const int m = X.length();
    const int n = Y.length();

    // initialize
    for (size_t i = 1; i <= m; i++) {
        dV(i, 0) = G;
    }
    for (size_t j = 1; j <= n; j++) {
        dH(0, j) = G;
    }

    // fill S with scores
    for (size_t i = 1; i <= m; i++) {
        for (size_t j = 1; j <= n; j++) {
            if (X[i - 1] == Y[j - 1]) {// match
                dV(i, j) = M - dH(i - 1, j);
                dH(i, j) = M - dV(i, j - 1);
            } else if (dH(i - 1, j) <= I - G &&
                       dV(i, j - 1) <= I - G) { // mismatch
                dV(i, j) = I - dH(i - 1, j);
                dH(i, j) = I - dV(i, j - 1);
            } else if (dH(i - 1, j) >= max(I - G, dV(i, j - 1))) { // indel from above
                dV(i, j) = G;
                dH(i, j) = G + dH(i - 1, j) - dV(i, j - 1);
            } else if (dV(i, j - 1) >= max(I - G, dH(i - 1, j))) { // indel from left
                dH(i, j) = G;
                dV(i, j) = G + dV(i, j - 1) - dH(i - 1, j);
            } else {
                throw runtime_error("unrecognized direction");
            }
        }
    }
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
    
    const auto X = sequences[0];
    const auto Y = sequences[1];

    const int m = X.length();
    const int n = Y.length();

    Matrix dV(m + 1, n + 1);
    Matrix dH(m + 1, n + 1);
    Matrix S(m + 1, n + 1);

    // populate dV, dH
    filldHdV(X, Y, dV, dH);
    // populate S
    fillS(X, Y, S);

    for (size_t i = 1; i <= m; i++) {
        for (size_t j = 1; j <= n; j++) {
            if (S(i, j) - S(i, j - 1) != dH(i, j) ||
                S(i, j) - S(i - 1, j) != dV(i, j)) {
                throw runtime_error("ΔH & ΔV are not matching S!! :(");
            }
        }
    }

    cout << "ΔH & ΔV match correctly with S :)" << endl;

    return EXIT_SUCCESS;
}
