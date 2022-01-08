#pragma once
/*
	\file:		linAlg.h
	\brief:		This file defines some vector and matrix computation classes.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#define _USE_MATH_DEFINES 
#include <math.h>
#include <iostream>

// convenient struct for 2d vector
struct PairDoub {
	double x;
	double y;

	// constructor
	PairDoub(double x_in, double y_in) :x(x_in), y(y_in) {}

	// copy constructor
	PairDoub(const PairDoub& z) : x(z.x), y(z.y) {}

	// default contructor
	PairDoub() :x(0), y(0) {}

	// calc norm
	double norm() const {
		return sqrt(x * x + y * y);	
	}

	double normsq() const {
		return x * x + y * y;
	}

	// dot product between two vectors
	friend double dot(
		const PairDoub& leftArg, 
		const PairDoub& rightArg
	) {
		return leftArg.x * rightArg.x + leftArg.y * rightArg.y;
	}

	PairDoub operator + (const PairDoub& rightArg) const {
		return PairDoub{ x + rightArg.x, y + rightArg.y };
	}

	PairDoub operator - (const PairDoub& rightArg) const {
		return PairDoub{ x - rightArg.x, y - rightArg.y };
	}

	PairDoub operator * (double rightArg) const {
		return PairDoub{ x * rightArg, y * rightArg };
	}

	PairDoub operator / (double rightArg) const {
		return PairDoub{ x / rightArg, y / rightArg };
	}

	friend PairDoub operator *(double leftArg, const PairDoub& rightArg) {
		return PairDoub{ leftArg * rightArg.x, leftArg * rightArg.y };
	}

	PairDoub perpUnit() {
		if (norm() == 0) 
			return PairDoub();
		else return 
			PairDoub(-y / norm(), x / norm());
	}

	// return vector (-y,x)
	PairDoub perp() {
		return PairDoub(-y, x);
	}

	// rotate vectory by theta
	PairDoub rotate(double theta) {
		return PairDoub(
			cos(theta) * x - sin(theta) * y, 
			sin(theta) * x + cos(theta) * y
		);
	}

	// return angle
	double angle() const {
		double theta = atan2(y, x);
		if (theta < 0) 
			theta += 2 * M_PI;
		return theta;
	}

	// compare angle with input vector
	bool operator<(const PairDoub& rightArg) {
		double tol = 1e-8;
		double my_ang = angle(), other_ang = rightArg.angle();
		if (my_ang < other_ang - tol)
			return true;
		else if (my_ang > other_ang + tol)
			return false;
		else
			return normsq() < rightArg.normsq();
	}
	void print() {
		std::cout << "(" << x << "," << y << ")";
	}
};

// struct for 2 x 2 matrix.
struct MatrixDoub {
	PairDoub row1, 
		row2;

	MatrixDoub(
		const PairDoub& row1In, 
		const PairDoub& row2In
	) : 
		row1(row1In), 
		row2(row2In) {}

	MatrixDoub() : 
		row1(PairDoub(0, 0)), 
		row2(PairDoub(0, 0)) {}

	MatrixDoub transpose() {
		return MatrixDoub(PairDoub(row1.x, row2.x), PairDoub(row1.y, row2.y));
	}

	MatrixDoub operator - (const MatrixDoub& rightArg) const {
		return MatrixDoub( 
			row1-rightArg.row1,
			row2-rightArg.row2
		);
	}

	MatrixDoub operator / (double rightArg) const {
		return MatrixDoub(
			row1 /rightArg,
			row2 / rightArg
		);
	}
	// multiply matrix by vector
	friend PairDoub rightMultiply(
		const MatrixDoub& leftArg, 
		const PairDoub& rightArg
	) {
		return PairDoub(
			dot(leftArg.row1, rightArg), 
			dot(leftArg.row2, rightArg)
		);
	}

	// multiply vector by matrix
	friend PairDoub leftMultiply(
		const PairDoub& leftArg, 
		const MatrixDoub& rightArg
	) {
		return PairDoub( 
			leftArg.x * rightArg.row1.x + leftArg.y * rightArg.row2.x, 
			leftArg.x * rightArg.row1.y + leftArg.y * rightArg.row2.y 
		);
	}
};

// Struct for 2 x 2 x 2 tensor.
struct TensorDoub {
	MatrixDoub mat1, 
		mat2;

	// return xAy where A is tensor 
	PairDoub sandwich(
		const PairDoub& left, 
		const PairDoub& right
	) {
		return PairDoub(
			dot(left, rightMultiply(mat1, right)), 
			dot(left, rightMultiply(mat2, right))
		);
	}
	TensorDoub(
		const MatrixDoub& mat1In, 
		const MatrixDoub& mat2In
	) : 
		mat1(mat1In), 
		mat2(mat2In) {}
	TensorDoub(
	) :
		mat1(MatrixDoub()), 
		mat2(MatrixDoub()) 
	{}
};

// vector of integers, for mesh indices
struct PairInt {
	int xi = 0;

	int yi = 0;

	PairInt(int x_in, int y_in) :xi(x_in), yi(y_in) {}

	PairInt(const PairInt& z) : xi(z.xi), yi(z.yi) {}

	double norm() const {
		return sqrt(xi * xi + yi * yi);
	}

	double normsq() const {
		return ((double)xi * xi + (double)yi * yi);
	}

	friend double dot(
		const PairInt& leftArg, 
		const PairDoub& rightArg
	) {
		return leftArg.xi * rightArg.x + leftArg.yi * rightArg.y;
	}

	friend double dot(
		const PairDoub& leftArg, 
		const PairInt& rightArg
	) {
		return leftArg.x * rightArg.xi + leftArg.y * rightArg.yi;
	}

	friend double dot(
		const PairInt& leftArg, 
		const PairInt& rightArg
	) {
		return static_cast<double>(leftArg.xi)* rightArg.xi 
			+ static_cast<double>(leftArg.yi)* rightArg.yi;
	}

	PairInt operator + (
		const PairInt& rightArg
		) const {
		return PairInt{ xi + rightArg.xi, yi + rightArg.yi };
	}

	PairInt operator - (
		const PairInt& rightArg
		) const {
		return PairInt{ xi - rightArg.xi, yi - rightArg.yi };
	}

	PairInt operator * (
		int rightArg
		) const {
		return PairInt{ xi * rightArg, yi * rightArg };
	}

	friend PairInt operator *(
		int leftArg, 
		const PairInt& rightArg
		) {
		return PairInt{ leftArg * rightArg.xi, leftArg * rightArg.yi };
	}

	void print() {
		std::cout << "(" << xi << "," << yi << ")";
	}

	double angle() const {
		double theta = atan2(static_cast<double>(yi), static_cast<double>(xi));
		if (theta < 0) theta += 2 * M_PI;
		return theta;
	}

	bool operator<(
		const PairInt& rightArg
		) {
		double tol = 1e-8;
		double my_ang = angle(), other_ang = rightArg.angle();
		if (my_ang < other_ang - tol)
			return true;
		else if (my_ang > other_ang + tol)
			return false;
		else
			return normsq() < rightArg.normsq();
	}

	bool operator==(
		const PairInt& rightArg
		) {
		return (xi == rightArg.xi && yi == rightArg.yi);
	}
}; 




