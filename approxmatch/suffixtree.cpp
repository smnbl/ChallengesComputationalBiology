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
#include <cassert>
#include <stack>
#include "suffixtree.h"

using namespace std;

// ============================================================================
// SUFFIX TREE (PRIVATE FUNCTIONS)
// ============================================================================

// ----------------------------------------------------------------------------
// ROUTINES TO MANIPULATE SUFFIX TREE POSITIONS
// ----------------------------------------------------------------------------

bool SuffixTree::advancePos(STPosition& pos, char c) const
{
        // case a) we are at a node: find a new child
        if (pos.atNode()) {
                STNode* chd = pos.node->getChild(c);
                if (chd == NULL)        // no such edge, get out
                        return false;
                pos.node = chd;
                pos.offset = 1;
                return true;
        }

        // case b) we are at an edge: try and match next character along edge
        if (T[pos.node->begin() + pos.offset] != c)
                return false;

        pos.offset++;
        return true;
}

bool SuffixTree::advancePos(STPosition& pos, const string& P,
                            size_t begin, size_t end) const
{
        for (auto itP = P.begin() + begin; itP < P.begin() + end; itP++)
                if (!advancePos(pos, *itP))
                        return false;

        return true;
}

vector<char> SuffixTree::getCharExtensions(const STPosition& pos) const
{
        vector<char> result;

        // case a) we are at a node: find a new child
        if (pos.atNode()) {
                for (int i = MAX_CHAR-1; i >= 0; i--) {
                        char c = (char)i;
                        if (pos.node->getChild(c) != NULL)
                                result.push_back(c);
                }
        } else
                result.push_back(T[pos.node->begin() + pos.offset]);

        return result;
}

void SuffixTree::advancePosSkipCount(STPosition& pos, const std::string& P,
                                     size_t begin, size_t end) const
{
        assert(begin <= end);

        // progress as much as possible on the current edge
        length_t adv = min<length_t>(end-begin, pos.node->getEdgeLength()-pos.offset);
        pos.offset += adv;
        begin += adv;

        while (begin < end) {
                // move to the correct child
                pos.node = pos.node->getChild(P[begin]);

                // progress as much as possible on the current edge
                adv = min<length_t>(end-begin, pos.node->getEdgeLength());
                pos.offset = adv;
                begin += adv;
        }
}

STPosition SuffixTree::followSuffixLink(const STPosition& pos) const
{
        // for root, do nothing
        if (pos.node == root)
                return pos;

        STNode* parent = pos.node->getParent();
        if (parent == root) {   // special case for edges originating from the root
                STPosition newPos(root);
                advancePosSkipCount(newPos, T, pos.node->begin()+1, pos.node->begin()+pos.offset);
                return newPos;
        } else {                // generic case (parent is an internal node)
                STPosition newPos(parent->getSuffixLink());
                advancePosSkipCount(newPos, T, pos.node->begin(), pos.node->begin()+pos.offset);
                return newPos;
        }
}

string SuffixTree::posToStr(const STPosition& pos) const
{
        string str;
        str.reserve(pos.getDepth());

        // create a stack with nodes from position (bottom) to root (top)
        stack<STNode*> path;
        path.push(pos.node);
        while (path.top()->getParent() != NULL)
                path.push(path.top()->getParent());

        path.pop();     // remove root (empty string);

        // traverse the stack from root (top) to position (bottom)
        while(!path.empty()) {
                STNode* node = path.top();
                path.pop();

                size_t len = path.empty() ? pos.offset : node->getEdgeLength();
                str.append(T.substr(node->begin(), len));
        }

        return str;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR PATTERN MATCHING
// ----------------------------------------------------------------------------

void SuffixTree::getOccurrences(const STPosition& pos, vector<size_t>& occ) const
{
        // if so, find all the leaves under "pos"
        stack<STNode*> stack;
        stack.push(pos.node);

        while (!stack.empty()) {
                STNode* node = stack.top();
                stack.pop();

                if (node->isLeaf())
                        occ.push_back(node->getSuffixIdx());
                else
                        for (int i = 0; i < MAX_CHAR; i++)
                                if (node->getChild(i) != NULL)
                                        stack.push(node->getChild(i));
        }
}

// ----------------------------------------------------------------------------
// ROUTINES TO MANIPULATE SUFFIX TREE
// ----------------------------------------------------------------------------

STPosition SuffixTree::splitEdge(const STPosition& pos)
{
        // check whether we are really in the middle of an edge
        if (pos.atNode())
                return pos;

        STNode *chd = pos.node;
        STNode *par = chd->getParent();
        STNode *mid = new STNode(chd->begin(), chd->begin() + pos.offset);
        chd->setBegin(mid->end());

        // set the correct pointers
        par->setChild(T[mid->begin()], mid);
        mid->setChild(T[chd->begin()], chd);

        return STPosition(mid);
}

void SuffixTree::addLeaf(const STPosition& pos, length_t suffixIndex)
{
        STNode *leaf = new STNode(suffixIndex + pos.getDepth(), T.size());
        leaf->setSuffixIdx(suffixIndex);
        pos.node->setChild(T[suffixIndex + pos.getDepth()], leaf);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR I/O
// ----------------------------------------------------------------------------

std::ostream& SuffixTree::write(std::ostream& o) const
{
        stack<pair<int, STNode*> > stack;
        stack.push(make_pair(0, root));

        while (!stack.empty()) {
                int depth = stack.top().first;
                STNode* node = stack.top().second;
                stack.pop();

                for (length_t i = 0; i < depth; i++)
                        cout << "  ";
                string nodeDescr = (node->isLeaf())? "LEAF " : "INTL ";
                if (node == root)
                        nodeDescr = "ROOT ";

                o << nodeDescr << "\""
                  << T.substr(node->begin(), node->getEdgeLength()) << "\""
                  << ", depth=" << node->getDepth();

                if (!node->isLeaf() && node != root && node->getSuffixLink() != NULL)
                        o << " (SL: \"" << posToStr(node->getSuffixLink()) << "\")";
                if (node->isLeaf())
                        o << " (" << node->getSuffixIdx() << ")";
                o << "\n";

                for (int i = MAX_CHAR-1; i >= 0; i--) {
                        char c = (char)i;
                        if (node->getChild(c) != NULL)
                                stack.push(make_pair(depth+1, node->getChild(c)));
                }
        }

        return o;
}

void SuffixTree::constructUkonen()
{
        // create root node with an empty range (it has no parent)
        root = new STNode(0, 0);

        // algorithm invariant: pos points to T[i:j-1[
        STPosition pos(root);

        // in phase j, build implicit suffix tree for prefix T[0:j[
        for (size_t j = 1, numLeaves = 0; j <= T.size(); j++) {
                STNode *prevInternal = NULL;

                // skip 'numleafs' times rule 1 (extension of leaf)
                for (size_t i = numLeaves; i < j; i++) {
                        // note that pos will always point to T[i:j-1[ at this point

                        // add a SL from the previously created internal node
                        if (prevInternal != NULL && pos.atNode()) {
                                prevInternal->setSuffixLink(pos.node);
                                prevInternal = NULL;
                        }

                        // try and add character T[j-1]
                        if (advancePos(pos, T[j-1]))
                                break; // rule 3 : do nothing + show stopper

                        // rule 2 : create internal node (optionally) and leaf
                        if (!pos.atNode()) {
                                pos = splitEdge(pos); // pos points to new node
                                if (prevInternal != NULL)
                                        prevInternal->setSuffixLink(pos.node);
                                prevInternal = pos.node;
                        }

                        addLeaf(pos, i);
                        numLeaves++;

                        // extension is complete: follow suffix link
                        pos = followSuffixLink(pos);
                }
        }
}

// ============================================================================
// SUFFIX TREE (PUBLIC FUNCTIONS)
// ============================================================================

SuffixTree::SuffixTree(const string& T) : T(T)
{
        // maximum string length = 2^32-1
        if (T.size() >= (size_t)numeric_limits<length_t>::max())
                throw runtime_error("String exceeds maximum length");

        // construct suffix tree using Ukonen's algorithm
        constructUkonen();

        // construct suffix tree using naive algorithm
        //constructNaive();
}

SuffixTree::~SuffixTree()
{
        // Depth-first traversal of the tree
        stack<STNode*> stack;
        stack.push(root);

        while (!stack.empty()) {
                STNode* node = stack.top();
                stack.pop();

                for (int i = 0; i < MAX_CHAR; i++)
                        if (node->getChild(i) != NULL)
                                stack.push(node->getChild(i));

                delete node;
        }
}

void SuffixTree::matchPattern(const string& P, vector<size_t>& occ)
{
        // clear the occurrence vector
        occ.clear();

        // can we completely match P?
        STPosition pos(root);
        if (!advancePos(pos, P, 0, P.size()))
                return;

        // get all occurrences under position
        getOccurrences(pos.node, occ);
}

void SuffixTree::recMatchApprox(STPosition pos, const string& P,
                                vector<size_t>& occ, BandMatrix& M)
{
        const int W = M.getWidth();

        // matrix range we're filling in: M[starti:endi, j]
        int j = pos.getDepth() + 1;
        int starti = max<int>(1, j - W);
        int endi = min<int>(P.size(), j + W);

        vector<char> extensions = getCharExtensions(pos);
        for (char c : extensions) {
                if (c == '#' || c == '$')       // don't pass # signs
                        continue;

                STPosition newPos = pos;        // only descend towards existing edges
                if (!advancePos(newPos, c))
                        continue;

                // compute the next edit distance matrix column
                int minVal = W + 1;
                for (int i = starti; i <= endi; i++) {
                        int diag = M(i-1, j-1) + (P[i-1] == c ? 0 : 1);
                        int gapX = (j > i - W) ? M(i, j-1) + 1 : diag + 1;
                        int gapY = (j < i + W) ? M(i-1, j) + 1 : diag + 1;

                        M(i, j) = min(min(diag, gapX), gapY);
                        minVal = min(minVal, M(i, j));
                }

                // edit distance threshold exceeded
                if (minVal > W)
                        continue;

                // have we reached the end of P?
                if (endi >= P.size())
                        if (M(P.size(), j) <= W ) {
                                getOccurrences(newPos.node, occ);
                                continue;
                        }

                // continue searching deeper
                recMatchApprox(newPos, P, occ, M);
        }
}

void SuffixTree::matchPatternApprox(const string& P, vector<size_t>& occ,
                                    int maxEditDist)
{
        // clear the occurrence vector
        occ.clear();

        STPosition pos(root);

        // intialize a band-diagonal matrix
        BandMatrix M(P.size() + 1, maxEditDist);
        for (int i = 0; i <= maxEditDist; i++) {
                M(0, i) = i;
                M(i, 0) = i;
        }

        recMatchApprox(pos, P, occ, M);
}
