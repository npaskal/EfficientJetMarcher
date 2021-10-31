#pragma once
/*
	\file:		binarytree.h
	\brief:		This file declares the binarytree class, a list containing
				a heapsort structure designed for the Considered list of the
				Fast Marching method.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include <vector>


// class containing a list of ints with heapsort structure
// TODO: refactor this
class BinaryTree {
public:
	std::vector<unsigned> tree;				// List of indices of mesh points

	const std::vector<double>* u;			// pointer to vector of u values 

	std::vector<int>* indexInConsidered;	// pointer to vector of indices

	// constructors -- no default defined
	BinaryTree(
		const std::vector<double>* u_in, 
		std::vector<int>* indexInConsidered_in
	): 
		u(u_in), 
		indexInConsidered(indexInConsidered_in) 
	{}

	// check if list is empty
	bool isEmpty() { return !(tree.size() > 0); }

	// swap points at indices of two arguments
	void swapPoints(unsigned swap1, unsigned swap2);


	unsigned minChild(unsigned parent);
	unsigned takeMin();
	void reassembleTopDown();
	void reassembleBottomUp();
	void reassembleFromTarget(unsigned index);
};
