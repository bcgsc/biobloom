#include "newick.hh"
#include "tree.hh"
using namespace BiRC::treelib;

#include <boost/spirit/include/classic_core.hpp>
#include <stack>
#include <iterator>
#include <istream>
#include <cassert>

namespace {
    struct newick_grammar : public boost::spirit::classic::grammar<newick_grammar>
    {
	mutable std::shared_ptr<Tree> tree;
	mutable std::stack<int> node_stack;
	mutable std::stack<double> branch_len_stack;
	
	struct handle_leaf
	{
	    newick_grammar &self;
	    handle_leaf(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr beg, Itr end) const
	    {
		std::string label(beg,end);
		self.node_stack.push(self.tree->add_node(label));
	    }
	};

	struct handle_inner_node
	{
	    newick_grammar &self;
	    handle_inner_node(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr, Itr) const
	    {
		int right = self.node_stack.top();
		self.node_stack.pop();
		double right_time = self.branch_len_stack.top();
		self.branch_len_stack.pop();
		int left = self.node_stack.top();
		self.node_stack.pop();
		double left_time = self.branch_len_stack.top();
		self.branch_len_stack.pop();

		self.node_stack.push(self.tree->add_node("",
							 left,right,
							 left_time,right_time));
	    }
	};

	struct handle_time_1
	{
	    newick_grammar &self;
	    handle_time_1(newick_grammar &self) : self(self) {}

	    void operator()(double d) const
	    {
		self.branch_len_stack.push(d);
	    }
	};

	struct handle_time_2
	{
	    newick_grammar &self;
	    handle_time_2(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr, Itr) const
	    {
		self.branch_len_stack.push(0);
	    }
	};


	struct handle_tree_1
	{
	    newick_grammar &self;
	    handle_tree_1(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr, Itr) const
	    {
		self.tree->set_root(self.node_stack.top());
	    }
	};

	struct handle_tree_2
	{
	    newick_grammar &self;
	    handle_tree_2(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr, Itr) const
	    {
		int right = self.node_stack.top();
		self.node_stack.pop();
		double right_time = self.branch_len_stack.top();
		self.branch_len_stack.pop();
		int left = self.node_stack.top();
		self.node_stack.pop();
		double left_time = self.branch_len_stack.top();
		self.branch_len_stack.pop();

		self.tree->set_root(self.tree->add_node("", 
							left,right,
							left_time,right_time));
	    }
	};

	struct handle_tree_3
	{
	    newick_grammar &self;
	    handle_tree_3(newick_grammar &self) : self(self) {}

	    template <typename Itr>
	    void operator()(Itr, Itr) const
	    {
		int right2 = self.node_stack.top(); self.node_stack.pop();
		int right1 = self.node_stack.top(); self.node_stack.pop();
		int left = self.node_stack.top(); self.node_stack.pop();

		int subtree = self.tree->add_node("",left,right1);
		self.tree->set_root(self.tree->add_node("",subtree,right2));
	    }
	};
	
	
	handle_leaf leaf_handler;
	handle_inner_node inner_node_handler;
	handle_time_1 time_1_handler;
	handle_time_2 time_2_handler;
	handle_tree_1 tree_1_handler;
	handle_tree_2 tree_2_handler;
	handle_tree_3 tree_3_handler;
	
	newick_grammar()
	    : tree(new Tree),
	      leaf_handler(*this),
	      inner_node_handler(*this),
	      time_1_handler(*this),
	      time_2_handler(*this),
	      tree_1_handler(*this),
	      tree_2_handler(*this),
	      tree_3_handler(*this)
	{
	}
	
	template <typename scanner_t>
	struct definition
	{
	    definition(newick_grammar const &self)
	    {
		using namespace boost::spirit::classic;
		leaf_r =
		    (alnum_p >> *alnum_p)[self.leaf_handler]
		    ;
		inner_node_r =
		    '(' >> edge_r >> ',' >> edge_r >> ')' 
			>> epsilon_p[self.inner_node_handler]
		    ;

		edge_r =
		    leaf_r >> ((':' >> (real_p[self.time_1_handler]))
			       | epsilon_p[self.time_2_handler])
		    |
		    inner_node_r >> ((':' >> (real_p[self.time_1_handler]))
				     | epsilon_p[self.time_2_handler])
		    ;

		// The split in start and end is because we need to
		// follow different rules depending on the number of
		// children of the root in the newick format, and we
		// cannot distinguish between 2 and 3 children from the
		// start of the parsing.
		tree_start_r =
		    leaf_r >> ';' >> epsilon_p[self.tree_1_handler]
		    |
		    '(' >> edge_r >> ',' >> edge_r >> tree_end_r
		    ;
		tree_end_r =
		    epsilon_p >> ')' >> ';' 
			>> epsilon_p[self.tree_2_handler]
		    |
		    epsilon_p >> ',' >> edge_r >> ')' >> ';' 
			>> epsilon_p[self.tree_3_handler]
		    ;
	    }
	    
	    boost::spirit::classic::rule<scanner_t> leaf_r, inner_node_r, edge_r;
	    boost::spirit::classic::rule<scanner_t> tree_start_r, tree_end_r;
	    const boost::spirit::classic::rule<scanner_t> &
	    start() const { return tree_start_r; }
	};   
    };
}

namespace {
    template <typename Itr>
    std::shared_ptr<Tree>
    parse(Itr beg, Itr end)
    {
	newick_grammar grammar;
	boost::spirit::classic::parse_info<Itr> info =
	    boost::spirit::classic::parse(beg, end, grammar, boost::spirit::classic::space_p);
	
	if (info.full) return grammar.tree;
	else           return std::shared_ptr<Tree>(0);
    }
}
    

std::shared_ptr<Tree>
BiRC::treelib::parse_newick(const std::string &str)
{
    return ::parse(str.begin(), str.end());
}


std::shared_ptr<Tree>
BiRC::treelib::parse_newick(const char *beg, const char *end)
{
    return ::parse(beg, end);
}



