#ifndef FM_TREE_H
#define FM_TREE_H
#include <vector>

class GridPointTree {
public:
	// variables
	std::vector<unsigned> tree;						// List of indices of grid points
	const std::vector<double>* u;					// pointer to vector of u values in master class
	std::vector<int>* indexInConsidered;

	// constructors -- no default defined
	GridPointTree(const std::vector<double>* u_in, std::vector<int>* indexInConsidered_in) : u(u_in), indexInConsidered(indexInConsidered_in) {}

	// member functions
	bool isEmpty() { return !(tree.size() > 0); }
	void swapPoints(unsigned swap1, unsigned swap2);
	unsigned minChild(unsigned parent);
	unsigned takeMin();
	void reassembleTopDown();
	void reassembleBottomUp();
	void reassembleFromTarget(unsigned index);
};
#endif