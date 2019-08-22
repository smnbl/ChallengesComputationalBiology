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

#include <iostream>
#include <fstream>
#include "suffixtree.h"

using namespace std;

void readInput(const string& filename, string& T)
{
        ifstream ifs(filename.c_str());

        while (ifs) {
                string line;
                getline(ifs, line);

                if (line.empty())
                        continue;

                if (!T.empty())
                        T.push_back('#');
                T.append(line);
        }

        T.push_back('$');
}

int main(int argc, char* argv[])
{
        if (argc != 2)
                throw runtime_error("Program requires an input text as an argument");

        cout << "Reading input..." << endl;
        string T;
        readInput(argv[1], T);

        cout << "Building suffix tree..." << endl;
        SuffixTree ST(T);

        while (true) {
                cout << "Type a word: ";
                string word;
                cin >> word;

                vector<AppMatch> appMatches;
                ST.matchPatternApprox(word, appMatches, 2);
                cout << "Found " << appMatches.size() << " occurrences." << endl;
                for (const auto appMatch : appMatches) {
                        size_t b = appMatch.begin;
                        size_t e = appMatch.end;
                        cout << "* [" << T.substr(b, e-b) << "] at edit distance " << appMatch.editDist;

                        while ( b > 0 && T[b-1] != '#')
                                b--;

                        while ( e < T.size() && T[e] != '#')
                                e++;

                        cout << " -- corresponding word: " << T.substr(b, e-b) << endl;
                }

        }

        cout << "Bye..." << endl;

        return EXIT_SUCCESS;
}
