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

void SuffixTree::reportMEM(const string& Q, size_t j, size_t minSize,
                           const STPosition& pos, vector<MEMOcc>& occ)
{
        vector<size_t> this_occ;
        getOccurrences(pos.node, this_occ);
        for (auto i : this_occ)
                if (j == 0 || i == 0 || Q[j-1] != T[i-1]) // left-maximal?
                        occ.push_back(MEMOcc(i, j, pos.getDepth()));

        STNode* node = pos.node->getParent();
        STNode* last = pos.node;

        while (node->getDepth() >= minSize) {
                for (size_t i = 0; i < MAX_CHAR; i++) {
                        STNode* chd = node->getChild(i);
                        if (chd == NULL || chd == last)
                                continue;

                        vector<size_t> this_occ;
                        getOccurrences(chd, this_occ);

                        for (auto i : this_occ)
                                if (i == 0 || j == 0 || Q[j-1] != T[i-1]) // left-maximal?
                                        occ.push_back(MEMOcc(i, j, node->getDepth()));
                }

                last = node;
                node = node->getParent();
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

length_t SuffixTree::recComputeSLPhase1(STNode* node, vector<STNode*>& A)
{
        // for leaves, simply return the suffix index
        if (node->isLeaf())
                return node->getSuffixIdx();

        length_t min1 = T.length() + 1; // the smallest suffix idx at this node
        length_t min2 = T.length() + 1; // the second smallest suffix idx

        for (size_t i = 0; i < MAX_CHAR; i++) {
                if (node->getChild(i) == NULL)
                        continue;

                length_t m = recComputeSLPhase1(node->getChild(i), A);
                if (m < min1) {
                        min2 = min1;
                        min1 = m;
                } else if (m < min2) {
                        min2 = m;
                }
        }

        // indicates that there is a suffix link from "node" to the internal
        // node that is caused by leaf with suffix index min2 + 1
        A[min2 + 1] = node;

        return min1;
}

void SuffixTree::recComputeSLPhase2(STNode* node, const vector<STNode*>& A,
                                    vector<STNode*>& B)
{
        if (node->isLeaf()) {
                length_t sufIdx = node->getSuffixIdx();
                if (A[sufIdx] != NULL && A[sufIdx] != root) {
                        length_t d = A[sufIdx]->getDepth();
                        A[sufIdx]->setSuffixLink(B[d-1]);
                }
        } else {
                length_t d = node->getDepth();
                B[d] = node;

                for (size_t i = 0; i < MAX_CHAR; i++)
                        if (node->getChild(i) != NULL)
                                recComputeSLPhase2(node->getChild(i), A, B);
        }
}

void SuffixTree::computeSuffixLinks()
{
        // see paper by Moritz G. Maass for definintions of cause() and branch()

        // compute A such that A[cause(node) + 1] = node
        vector<STNode*> A(T.size(), NULL);
        recComputeSLPhase1(root, A);

        // compute B such that B[depth] = branch(node, depth)
        vector<STNode*> B(T.size(), NULL);
        recComputeSLPhase2(root, A, B);
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

void SuffixTree::constructNaive()
{
        // create a node with an empty range (it has no parent)
        root = new STNode(0, 0);

        for (size_t i = 0; i < T.size(); i++) {
                // find position in ST that maximally matches the prefix of suf_i(T)
                STPosition pos(root);

                // if suf_i(T) is completely matched, T is not prefix-free
                if (advancePos(pos, T, i, T.size()))
                        break;
                        //throw runtime_error("Text T is not prefix-free");

                // middle of an edge? -> first split edge first
                if (!pos.atNode())
                        pos = splitEdge(pos);

                // add a new leaf to the suffix tree
                addLeaf(pos, i);
        }

        // Maass algorithm to compute suffix links in O(n) time
        computeSuffixLinks();
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

void SuffixTree::findMEM(const string& Q, size_t minSize, vector<MEMOcc>& occ)
{
        // clear the occurrence vector
        occ.clear();

        STPosition pos(root);
        for (size_t j = 0; j < Q.size(); j++) {

                // advance posMax as much as possible
                advancePos(pos, Q, j+pos.getDepth(), Q.size());

                if (pos.getDepth() >= minSize)
                        reportMEM(Q, j, minSize, pos, occ);

                pos = followSuffixLink(pos);
        }
}
