#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include "binarytree.h"

void GridPointTree::swapPoints(unsigned swap1, unsigned swap2) {
	unsigned point2 = tree[swap2];
	(*indexInConsidered)[tree[swap1]] = swap2;
	(*indexInConsidered)[tree[swap2]] = swap1;
	tree[swap2] = tree[swap1];
	tree[swap1] = point2;
}

unsigned GridPointTree::minChild(unsigned parent) { // returns TREE index of minimum Child
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
		if ((*u)[tree[parent]] <= std::min((*u)[tree[2 * parent + 1]], (*u)[tree[2 * parent + 2]]))
			return parent;
		else if ((*u)[tree[2 * parent + 1]] <= (*u)[tree[2 * parent + 2]])
			return 2 * parent + 1;
		else return 2 * parent + 2;
	}
}
unsigned GridPointTree::takeMin() {
	unsigned minimum = tree[0];
	swapPoints(0, (unsigned)tree.size() - 1);
	tree.pop_back();
	reassembleTopDown();
	return minimum;
}
void GridPointTree::reassembleTopDown() {
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
void GridPointTree::reassembleBottomUp() {
	unsigned curr = (unsigned)tree.size() - 1; // start at bottom
	while (curr > 1 && (*u)[tree[curr]] < (*u)[tree[(curr - 1) / 2]]) {
		swapPoints(curr, (curr - 1) / 2);
		curr = (curr - 1) / 2;
	}
}
void GridPointTree::reassembleFromTarget(unsigned index) {
	unsigned curr = index; 
	while (curr >= 1 && (*u)[tree[curr]] < (*u)[tree[(curr - 1) / 2]]) {
		swapPoints(curr, (curr - 1) / 2);
		curr = (curr - 1) / 2;
	}
}
