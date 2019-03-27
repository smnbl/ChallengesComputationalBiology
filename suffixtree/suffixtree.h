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

#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <tuple>
#include <limits>
#include <string>
#include <vector>
#include <cassert>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint32_t length_t;

// <position in T, position in Q, length of MEM>
typedef std::tuple<size_t, size_t, size_t> MEMOcc;

#define MAX_CHAR 256

// ============================================================================
// CLASS SUFFIX TREE NODE
// ============================================================================

// A suffix tree node contains links to its parent and children, a suffix link,
// the depth of the node and a suffix index in case it is a leaf node.
// For non-leaf nodes, the suffix index contains the length_t::max() value.

// A suffix tree node also contains a range [beginIdx, endIdx[ in T of its
// parent edge. The range encodes the characters implied on the edge.

class STNode {

private:
        // edge properties (between parent and current node)
        length_t beginIdx;              // begin index in T of parent edge
        length_t endIdx;                // end index in T of parent edge

        // node properties
        STNode* parent;                 // pointer to the parent (NULL for root)
        STNode* child[MAX_CHAR];        // pointers to children
        STNode* suffixLink;             // points to suffix link node
        length_t depth;                 // depth of current node
        length_t suffixIdx;             // suffix index (only for leaf nodes)

public:
        /**
         * Constructor
         * @param begin Begin index in T
         * @param end End index in T
         */
        STNode(length_t begin, length_t end) : beginIdx(begin), endIdx(end) {
                parent = NULL;
                for (int i = 0; i < MAX_CHAR; i++)
                        child[i] = NULL;
                suffixLink = NULL;
                depth = end - begin;
                suffixIdx = std::numeric_limits<length_t>::max();
        }

        /**
         * Get the begin index in T of parent edge
         * @return Begin index
         */
        length_t begin() const {
                return beginIdx;
        }

        /**
         * Set the begin index in T of parent edge
         * @param target target value
         */
        void setBegin(length_t target) {
                beginIdx = target;
        }

        /**
         * Get the end index in T of parent edge
         * @return End index
         */
        length_t end() const {
                return endIdx;
        }

        /**
         * Set the begin index in T of parent edge
         * @param target target value
         */
        void setEnd(length_t target) {
                endIdx = target;
        }

        /**
         * Get the length of the edge between parent en this node
         * @return Edge length
         */
        length_t getEdgeLength() const {
                return endIdx - beginIdx;
        }

        /**
         * Get the parent node
         * @return The parent node
         */
        STNode* getParent() const {
                return parent;
        }

        /**
         * Get a pointer to child node for which the edge starts with c
         * @param c Character c
         */
        STNode* getChild(char c) const {
                return child[static_cast<unsigned char>(c)];
        }

        /**
         * Set a (pre-allocated) child for this node
         * @param c First character on the edge to the child
         * @param chdToAdd Pointer to the child
         */
        void setChild(char c, STNode* chdToAdd) {
                chdToAdd->parent = this;
                chdToAdd->depth = depth + chdToAdd->getEdgeLength();
                child[static_cast<unsigned char>(c)] = chdToAdd;
        }

        /**
         * Set the suffix link (for internal nodes only)
         * @param target Target value
         */
        void setSuffixLink(STNode* target) {
                suffixLink = target;
        }

        /**
         * Get the suffix link from this node (only for internal nodes)
         * @return suffix link
         */
        STNode* getSuffixLink() const {
                return suffixLink;
        }

        /**
         * Set the suffix index (for leaves only)
         * @param target Target value
         */
        void setSuffixIdx(length_t target) {
                suffixIdx = target;
        }

        /**
         * Get the suffix index in T (for leaves only)
         * @return suffix index
         */
        length_t getSuffixIdx() const {
                return suffixIdx;
        }

        /**
         * Set the depth of the node
         * @param target Target value
         */
        void setDepth(length_t target) {
                depth = target;
        }

        /**
         * Get the depth of the node
         * @return Depth of node
         */
        length_t getDepth() const {
                return depth;
        }

        /**
         * Check whether the node is a leaf
         * @return true or false
         */
        bool isLeaf() const {
                return suffixIdx != std::numeric_limits<length_t>::max();
        }
};

// ============================================================================
// CLASS SUFFIX TREE POSITION
// ============================================================================

// A position in a suffix tree is a pointer to some node of the (implied)
// uncompressed suffix trie. In other words, a position could point to a leaf
// or an internal node in the suffix tree, or to a position halfway an edge.
// In that case, the position is specified by a pointer to the child node of
// that edge and an offset along the edge.

class STPosition {

private:
        STNode* node;                   // pointer to child node
        length_t offset;                // offset > 0 along the edge

public:
        /**
         * Constructor for position that points to a node itself
         * @param node Pointer to node in the ST
         */
        STPosition(STNode* node) : node(node), offset(node->getEdgeLength()) {}

        /**
         * Constructor for position that points to an edge
         * @param node Pointer to the child node
         * @param offset Offset along edge from parent to child
         */
        STPosition(STNode* node, length_t offset) : node(node), offset(offset) {
                assert(offset <= node->getEdgeLength());
        }

        /**
         * Check whether position points to a node (true) or edge (false)
         * @return True or false
         */
        bool atNode() const {
                return offset == node->getEdgeLength();
        }

        /**
         * Get the depth of the position
         * @return Depth of the position
         */
        length_t getDepth() const {
                // do not use node->getParent()->depth as root has no parent
                return node->getDepth() - node->getEdgeLength() + offset;
        }

        /**
         * Operator== overloading
         * @param rhs Right hand side
         * @return true or false
         */
        bool operator==(const STPosition& rhs) {
                if (node != rhs.node)
                        return false;
                if (offset != rhs.offset)
                        return false;
                return true;
        }

        /**
         * Operator != overloading
         * @param rhs Right hand side
         * @return true or false
         */
        bool operator!=(const STPosition& rhs) {
                return !(*this == rhs);
        }

        friend class SuffixTree;
};

// ============================================================================
// CLASS SUFFIX TREE
// ============================================================================

class SuffixTree {

private:
        // --------------------------------------------------------------------
        // ROUTINES TO MANIPULATE SUFFIX TREE POSITIONS
        // --------------------------------------------------------------------

        /**
         * Given a position, try and advance by matching a single character
         * @param pos Position to advance (input / output)
         * @param c Character c to match
         * @return true if position was advanced, false otherwise
         */
        bool advancePos(STPosition& pos, char c) const;

        /**
         * Given a position, try and advance by matching a pattern
         * This procedure advances a much as possible, even if P cannot be
         * entirely matched
         * @param pos Position to advance (input / output)
         * @param P Pattern to match str[begin, end[
         * @param begin First position in str to match
         * @param end End position in str to match
         * @return true if pattern is fully matched, false otherwise
         */
        bool advancePos(STPosition& pos, const std::string& P,
                        size_t begin, size_t end) const;

        /**
         * Faster matching using the skip/count trick (patterns *must* exist!)
         * @param pos Position to advance (input / output)
         * @param P Pattern to match str[begin, end[
         * @param begin First position in str to match
         * @param end End position in str to match
         */
        void advancePosSkipCount(STPosition& pos, const std::string& P,
                                 size_t begin, size_t end) const;

        /**
         * Given a suffix tree position, get the position by following the SL
         * @param pos Input position
         * @return Position after following the suffix link
         */
        STPosition followSuffixLink(const STPosition& pos) const;

        /**
         * Get the implied string from a suffix tree position
         * @param pos Suffix tree position
         * @return Implied string from root to pos
         */
        std::string posToStr(const STPosition& pos) const;

        // --------------------------------------------------------------------
        // ROUTINES FOR PATTERN MATCHING
        // --------------------------------------------------------------------

        /**
         * Get all occurrences under a given position
         * @param node Pointer to a node
         * @param occ Vector to store the occurrences (output)
         */
        void getOccurrences(const STPosition& pos, std::vector<size_t>& occ) const;

        /**
         * Find the Maximal Exact Matches between T and P
         * @param Q Query pattern
         * @param i Current suffix of Q under consideration
         * @param minSize Minimum size of the MEM
         * @param pos Position in ST that points to RMEM of suf_i(Q)
         * @param occ MEM occurrences (output)
         */
        void reportMEM(const std::string& Q, size_t i, size_t minSize,
                       const STPosition& pos, std::vector<MEMOcc>& occ);

        // --------------------------------------------------------------------
        // ROUTINES TO CONSTRUCT/MANIPULATE SUFFIX TREE
        // --------------------------------------------------------------------

        /**
         * Given a suffix tree position, split the edge
         * @param pos Suffix tree position
         * @return Position of the edge has been split
         */
        STPosition splitEdge(const STPosition& pos);

        /**
         * Add a leaf to the suffix tree
         * @param pos Suffix tree position
         * @param suffixIndex
         */
        void addLeaf(const STPosition& pos, length_t suffixIndex);

        /**
         * Step 1 of the Maass algorithm for the computation of suffix links
         * @param node Current node
         * @param A Vector that maps causing leaf idx to the source of a SL
         * @return The minimum suffix index for the subtree under node
         */
        length_t recComputeSLPhase1(STNode* node, std::vector<STNode*>& A);

        /**
         * Step 2 of the Maass algorithm for the computation of suffix links
         * @param node Current node
         * @param A Vector that maps causing leaf idx to the source of a SL
         * @param B Vector that maps depth to node on the current path
         */
        void recComputeSLPhase2(STNode* node, const std::vector<STNode*>& A,
                                std::vector<STNode*>& B);

        /**
         * Compute suffix links in linear time using the Maass algorithm
         */
        void computeSuffixLinks();

        /**
         * Construct the suffix tree using a naive O(n^2) time algorithm
         */
        void constructNaive();

        /**
         * Construct the suffix tree using Ukonen's algorithm in O(n) time
         */
        void constructUkonen();

        // --------------------------------------------------------------------
        // ROUTINES FOR I/O
        // --------------------------------------------------------------------

        /**
         * Write suffix tree to the output stream in a formatted way
         * @param o Output stream
         * @return Output stream
         */
        std::ostream& write(std::ostream& o) const;

        // --------------------------------------------------------------------
        const std::string T;            // text to index
        STNode* root;                   // pointer to the root node
        // --------------------------------------------------------------------

public:
        /**
         * Constructor
         * @param T Text to be indexed
         */
        SuffixTree(const std::string& T);

        /**
         * Destructor
         */
        ~SuffixTree();

        /**
         * Find occurrences of P in the suffix tree in O(m + #occ) time
         * @param P Pattern P to match (length m)
         * @param occ Start positions of the occurrences in T (output)
         */
        void matchPattern(const std::string& P, std::vector<size_t>& occ);

        /**
         * Find the Maximal Exact Matches between T and P in O(m + #RMEMs) time
         * @param Q Query pattern Q (length m)
         * @param minSize Minimum size of the MEM
         * @param occ MEM occurrences (output)
         */
        void findMEM(const std::string& Q, size_t minSize,
                     std::vector<MEMOcc>& occ);

        /**
         * Operator<< overloading
         * @param o Output stream
         * @param ST Suffix tree
         * @return Output stream with ST information
         */
        friend std::ostream& operator<< (std::ostream& o, const SuffixTree& t) {
                return t.write(o);
        }
};

#endif
