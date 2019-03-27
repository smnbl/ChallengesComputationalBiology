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

                vector<size_t> occ;
                ST.matchPattern(word, occ);
                cout << "Found " << occ.size() << " occurrences." << endl;
                for (auto pos : occ) {
                        size_t begin = pos;
                        while (begin > 0 && T[begin-1] != '#')
                                begin--;

                        size_t end = pos;
                        while (end+1 < T.size() && T[end+1] != '#')
                                end++;

                        cout << "* " << T.substr(begin, end-begin+1) << endl;
                }

        }

        cout << "Bye..." << endl;

        return EXIT_SUCCESS;
}
