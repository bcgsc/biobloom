#ifndef NEWICK_HH_INCLUDED
#define NEWICK_HH_INCLUDED

#include <string>
#include <memory>


namespace BiRC { namespace treelib {

    // forward decl.
    class Tree;

    /**
     * Parse a string in Newick format into a tree.
     *
     * \param str	The string to be parsed.
     * \returns	A tree corresponding to the string, or 0 in case of errors.
     */
    std::shared_ptr<Tree> parse_newick(const std::string &str);

    /**
     * Parse a string in Newick format into a tree.
     *
     * \param beg	Iterator to the beginning of the string.
     * \param end	Iterator to the end of the string.
     * \returns	A tree corresponding to the string, or 0 in case of errors.
     */
    std::shared_ptr<Tree> parse_newick(const char *beg, const char *end);
}}

#endif
