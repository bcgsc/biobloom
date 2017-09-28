#ifndef TREE_HH_INCLUDED
#define TREE_HH_INCLUDED

#include <string>
#include <vector>
#include <iosfwd>
#include <cassert>
#include <google/dense_hash_map>

namespace BiRC { namespace treelib {
    /**
     * \brief Data structure for binary trees.
     *
     * There are essentially no consistency guarantees in this
     * structure except that if the structure is a correct tree before
     * a method call, then it will always be a consistent tree after a
     * method call.  Stronger guarantees about consistency would put
     * too strong restrictions on the tree construction interface.
     */
    class Tree {
	
	/// Mapping of a node's index to its parents index.
	std::vector<int> parent_;
	/// Mapping of a node's index to its left child's index.
	std::vector<int> left_child_;
	/// Mapping of a node's index to its right child's index.
	std::vector<int> right_child_;

	/// Mapping of a node's index to its label.
	std::vector<std::string> label_;
	google::dense_hash_map<std::string, int> labelHash_;

	/// Branch lengths (measured from a node to its parent)
	std::vector<double> length_to_parent_;

	/// Index of the root of the tree.
	int root_;

	/// Helper function for printing
	void dfs_print(std::ostream &os, int node) const;

    public:

	/// Construct a new empty tree.
	Tree() : root_(-1) {labelHash_.set_empty_key("");};

	// --- tree statistics... ----------------------------------
	/**
	 * \brief Get the number of nodes in the tree
	 *
	 * Notice, this is the total number of nodes in the tree,
	 * <em>not</em> the number of leaves in the the tree.
	 *
	 * \returns Number of leaves in the tree.
	 */
	int size() { return parent_.size(); }

	// --- predicates ------------------------------------------

	/// Predicate testing if \a node is a leaf.
	bool is_leaf(int node) const
	{ 
	    assert(node < (int)left_child_.size());
	    return left_child_[node] < 0 and right_child_[node] < 0;
	}

	/// Predicate testing if \a node is an inner node.
	bool is_inner(int node) const
	{ 
	    assert(node < (int)left_child_.size());
	    return left_child_[node] >= 0 and right_child_[node] >= 0;
	}


	/// Predicate testing if \a node is the root.
	bool is_root(int node) const
	{ 
	    return node == root_;
	}

	// --- node access -----------------------------------------

	int root() const { return root_; }

	/// Get the label associated with \a node.
	const std::string &label(int node) const
	{
	    assert(node < (int)label_.size());
	    return label_[node];
	}

	void setLabel(int node, std::string val)
	{
		labelHash_[val] = node;
//		labelHash_.erase(val);
	    label_[node] = val;
	}

	int node(const std::string &label) const{
		return labelHash_.find(label)->second;
	}

	/// Get the left child of \a node.
	int left_child(int node) const
	{
	    assert(node < (int)left_child_.size());
	    return left_child_[node];
	}

	/// Get the right child of \a node.
	int right_child(int node) const
	{
	    assert(node < (int)right_child_.size());
	    return right_child_[node];
	}

	/// Get the parent of \a node.
	int parent(int node) const
	{
	    assert(node < (int)parent_.size());
	    return parent_[node];
	}

	double length_to_parent(int node) const
	{
	    assert(node < (int)length_to_parent_.size());
	    return length_to_parent_[node];
	}

	// --- construction ----------------------------------------

	/**
	 * \brief Insert a new node in the tree.
	 *
	 * Adds a node to the tree.  If children are specified, their
	 * parent index is updated (potentially deleting their
	 * previous parent).
	 *
	 * \param label		Label to give the node, if any.
	 * \param left_child	Index of left child, if any (if >= 0).
	 * \param right_child	Index of right child, if any (if >= 0).
	 * \param left_length   Branchlength from the left child to this node.
	 * \param left_length   Branchlength from the right child to this node.
	 *
	 * \returns		Index of the new node.
	 */
	int add_node(const std::string &label,
		     int left_child = -1, int right_child = -1,
		     double left_length = 0, double right_length = 0);

	/**
	 * \brief Sets the root index.
	 *
	 * The \a root index will be considered the root of the tree.
	 * There is no guarantee that the sub-tree rooted in \a root
	 * will contain the full structure; that must be guaranteed by
	 * the tree building code.
	 *
	 * \param root 		Index to be used as root.
	 *
	 * \pre \a root must be an index containing an actual node,
	 * i.e. one that was returned by add_node().
	 */
	void set_root(int root)
	{
	    assert(root < (int)parent_.size());
	    root_ = root;
	}

	
	// --- misc. -----------------------------------------------

	/**
	 * \brief Prints the tree structure to a stream.
	 *
	 * Prints the tree, both in its internal representation and in
	 * Newick format.
	 *
	 * If the tree is rooted in a leaf, that leaf will be printed
	 * as a left-most leaf in the Newick representation.
	 */
	void print(std::ostream &os) const;
    };

    inline std::ostream &operator<<(std::ostream &os, const Tree &tree)
    {
	tree.print(os);
	return os;
    }
}}



#endif
