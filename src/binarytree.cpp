/*
	\file:		binarytree.cpp
	\brief:		This file defines the binarytree class, a list containing
				a heapsort structure designed for the Considered list of the
				Fast Marching method.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <algorithm>
#include <cmath>
#include "binarytree.h"

// swap points at indices
void BinaryTree::swapPoints(unsigned swap1, unsigned swap2) {
	unsigned point2 = tree[swap2];
	(*indexInConsidered)[tree[swap1]] = swap2;
	(*indexInConsidered)[tree[swap2]] = swap1;
	tree[swap2] = tree[swap1];
	tree[swap1] = point2;
}


unsigned BinaryTree::minChild(unsigned parent) { 
	// if no children, return parent
	if (2 * parent + 1 > tree.size() - 1)
		return parent; 
	// if 1 child 
	else if (2 * parent + 2 > tree.size() - 1) {
		if ((*u)[tree[parent]] <= (*u)[tree[2 * parent + 1]])
			return parent;
		else return 2 * parent + 1;
	}
	// if 2 children
	else {
		if ((*u)[tree[parent]] <= std::min(
			(*u)[tree[2 * parent + 1]], 
			(*u)[tree[2 * parent + 2]]
		) )
			return parent;
		else if ((*u)[tree[2 * parent + 1]] <= (*u)[tree[2 * parent + 2]])
			return 2 * parent + 1;
		else return 2 * parent + 2;
	}
}

unsigned BinaryTree::takeMin() {
	unsigned minimum = tree[0];
	swapPoints(0, (unsigned)tree.size() - 1);
	tree.pop_back();
	reassembleTopDown();
	return minimum;
}

void BinaryTree::reassembleTopDown() {
	unsigned rounds = (tree.size() > 1) ? (unsigned) log2(tree.size()) : 0;
	unsigned currIndex = 0, temp = 0, round = 0;
	while (round < rounds) {
		temp = minChild(currIndex);
		if (temp == currIndex)
			break;
		else {
			swapPoints(currIndex, temp);
			currIndex = temp;
			round++;
		}
	}
}

void BinaryTree::reassembleBottomUp() {
	unsigned curr = (unsigned)tree.size() - 1; // start at bottom
	while (curr > 1 && (*u)[tree[curr]] < (*u)[tree[(curr - 1) / 2]]) {
		swapPoints(curr, (curr - 1) / 2);
		curr = (curr - 1) / 2;
	}
}

void BinaryTree::reassembleFromTarget(unsigned index) {
	unsigned curr = index; 
	while (curr >= 1 && (*u)[tree[curr]] < (*u)[tree[(curr - 1) / 2]]) {
		swapPoints(curr, (curr - 1) / 2);
		curr = (curr - 1) / 2;
	}
}
