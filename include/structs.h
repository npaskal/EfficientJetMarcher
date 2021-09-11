#ifndef FM_STRUCTS_H
#define FM_STRUCTS_H
#define _USE_MATH_DEFINES 
#include <math.h>
#include <assert.h>
#include <tuple>
#include <iostream>
#include <algorithm>
#include <vector>


// ************ Linear Algebra Structures *************

/*
	My ordered pair struct, endowed with relative linear algebra.
*/
struct PairDoub {
	double x;
	double y;

	// Initializers:
	PairDoub(double x_in, double y_in) :x(x_in), y(y_in) {}
	PairDoub(const PairDoub& z) : x(z.x), y(z.y) {}
	PairDoub() :x(0), y(0) {}

	// Vector functions:
	double norm() const { return sqrt(x * x + y * y); }
	double normsq() const {	return x * x + y * y;}
	friend double dot(const PairDoub& leftArg, const PairDoub& rightArg) { return leftArg.x * rightArg.x + leftArg.y * rightArg.y; }
	PairDoub operator + (const PairDoub& rightArg) const { return PairDoub{ x + rightArg.x, y + rightArg.y }; }
	PairDoub operator - (const PairDoub& rightArg) const { return PairDoub{ x - rightArg.x, y - rightArg.y }; }
	PairDoub operator * (double rightArg) const { return PairDoub{ x * rightArg, y * rightArg }; }
	PairDoub operator / (double rightArg) const { return PairDoub{ x / rightArg, y / rightArg }; }
	friend PairDoub operator *(double leftArg, const PairDoub& rightArg) { return PairDoub{ leftArg * rightArg.x, leftArg * rightArg.y }; }
	PairDoub perpUnit() {
		if (norm() == 0) return PairDoub();
		else return PairDoub(-y / norm(), x / norm());
	}
	PairDoub perp() {return PairDoub(-y, x);}
	PairDoub rotate(double theta) {
		return PairDoub(cos(theta) * x - sin(theta) * y, sin(theta) * x + cos(theta) * y);
	}
	double angle() const {
		double theta = atan2(y,x);
		if (theta < 0) theta += 2 * M_PI;
		return theta;
	}
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
	void print() { std::cout << "(" << x << "," << y << ")"; }
};

// Struct for 2 x 2 matrix.
struct MatrixDoub {
	PairDoub row1, row2;
	MatrixDoub(const PairDoub& row1In, const PairDoub& row2In) : row1(row1In), row2(row2In) {}
	MatrixDoub() : row1(PairDoub(0,0)),row2(PairDoub(0,0)) {}
	friend PairDoub rightMultiply(const MatrixDoub& leftArg, const PairDoub& rightArg) { return PairDoub(dot(leftArg.row1,rightArg),dot(leftArg.row2,rightArg)); }
	friend PairDoub leftMultiply(const PairDoub& leftArg, const MatrixDoub& rightArg) {
		return PairDoub{ leftArg.x * rightArg.row1.x + leftArg.y * rightArg.row2.x, leftArg.x * rightArg.row1.y + leftArg.y * rightArg.row2.y };
	}
};

// Struct for 2 x 2 x 2 tensor.
struct TensorDoub {
	MatrixDoub mat1, mat2;
	PairDoub sandwich(const PairDoub& left, const PairDoub& right) {
		return PairDoub(dot(left, rightMultiply(mat1, right)), dot(left, rightMultiply(mat2, right)));
	}
	TensorDoub(const MatrixDoub& mat1In, const MatrixDoub& mat2In): mat1(mat1In), mat2(mat2In) {}
	TensorDoub() :mat1(MatrixDoub()), mat2(MatrixDoub()) {}
};

struct PairInt {
	int xi = 0;
	int yi = 0;
	PairInt(int x_in, int y_in) :xi(x_in), yi(y_in) {}
	PairInt(const PairInt& z) : xi(z.xi), yi(z.yi) {}

	double norm() const { return sqrt(xi * xi + yi * yi); }
	double normsq() const { return ((double)xi * xi + (double)yi * yi); }
	friend double dot(const PairInt& leftArg, const PairDoub& rightArg) { return leftArg.xi * rightArg.x + leftArg.yi * rightArg.y; }
	friend double dot(const PairDoub& leftArg, const PairInt& rightArg) { return leftArg.x * rightArg.xi + leftArg.y * rightArg.yi; }
	friend double dot(const PairInt& leftArg, const PairInt& rightArg) { return static_cast<double>(leftArg.xi) * rightArg.xi + static_cast<double>(leftArg.yi) * rightArg.yi; }
	PairInt operator + (const PairInt& rightArg) const { return PairInt{ xi + rightArg.xi, yi + rightArg.yi }; }
	PairInt operator - (const PairInt& rightArg) const { return PairInt{ xi - rightArg.xi, yi - rightArg.yi }; }
	PairInt operator * (int rightArg) const { return PairInt{ xi * rightArg, yi * rightArg }; }
	friend PairInt operator *(int leftArg, const PairInt& rightArg) { return PairInt{ leftArg * rightArg.xi, leftArg * rightArg.yi }; }
	void print() { std::cout << "(" << xi << "," << yi << ")"; }
	double angle() const {
		double theta = atan2(static_cast<double>(yi), static_cast<double>(xi));
		if (theta < 0) theta += 2 * M_PI;
		return theta;
	}
	bool operator<(const PairInt& rightArg) {
		double tol = 1e-8;
		double my_ang = angle(), other_ang = rightArg.angle();
		if (my_ang < other_ang - tol)
			return true;
		else if (my_ang > other_ang + tol)
			return false;
		else
			return normsq() < rightArg.normsq();
	}
	bool operator==(const PairInt& rightArg) {
		return (xi == rightArg.xi && yi == rightArg.yi);
	}
};

// ***************** Fast March Structures *****************
enum algorithm {
	FM = 0,
	ASR = 1,
	OLIM = 2
};

enum group {
	FAR = 0,
	CONSIDERED = 1,
	FRONT = 2,
	ACCEPTED = 3
};
enum initShape {
	ELLIPSE = 0,
	BOX = 1
};

struct StencilData {
	double alpha;
	int stencil_cutoff;
	int num_theta_cells;
	bool boxUpdate;
	bool bubble;
	std::vector<double> theta_res;
	StencilData(double alphaIn, int maxStencIn, int binsIn, bool boxUpIn, bool bubbleIn ) : 
		alpha(alphaIn), stencil_cutoff(maxStencIn), num_theta_cells(binsIn), boxUpdate(boxUpIn), bubble(bubbleIn) {
		createStencilResolution_new();
	}
	void createStencilResolution() {
		theta_res.push_back(acos((alpha + sqrt(alpha * alpha + 8)) / 4.0));
		int k = 1;
		while (1) {
			double fac = 1;
			double inc = acos(std::max(alpha*fac * cos(theta_res[k-1]), 0.0));
			if (theta_res[k - 1] + inc > M_PI) {
				theta_res.push_back(M_PI);
				break;
			}
			else {
				theta_res.push_back(theta_res[k - 1] + inc / 2);
				theta_res.push_back(theta_res[k - 1] + inc);

			}
			k += 2;
		}
		int j = k - 1;
		k++;
		while (j >= 0) {
			theta_res.push_back(2 * M_PI - theta_res[j]);
			k++; j--;
		}

	}

	void createStencilResolution_new() {
		// Note this also puts the midpoint angle in between pi/2^k and pi/2^{k+1}
		int kmax = 8;
		for (int k = kmax; k >= 1; k--) {
			theta_res.push_back(M_PI /  pow(2, k));
			if (k > 1) theta_res.push_back(M_PI * 3.0 / pow(2, k + 1));
		}
		for (int k = 1; k <= kmax; k++) {
			if (k > 1) theta_res.push_back(2 * M_PI - M_PI * 3.0 / pow(2, k + 1));
			theta_res.push_back(2*M_PI - M_PI / pow(2, k));

		}

	}
};

struct GridInfo {
public:
	int nx;
	double xl, xr;
	int ny;
	double yl, yr;
	double hx, hy;
	const std::vector<PairInt> eightPtNeighborShifts;
	const std::vector<PairInt> fourPtNeighborShifts;

	GridInfo(int nx_i, double xl_i, double xr_i, int ny_i, double yl_i, double yr_i) :
		nx(nx_i), xl(xl_i), xr(xr_i), ny(ny_i), yl(yl_i), yr(yr_i),
		hx(double(xr_i - xl_i) / (static_cast<double>(nx_i) - 1)),
		hy(double(yr_i - yl_i) / (static_cast<double>(ny_i) - 1)),
		eightPtNeighborShifts({ PairInt(1,0),PairInt(1,1),	PairInt(0,1),PairInt(-1,1),
			PairInt(-1,0),PairInt(-1,-1),PairInt(0,-1),PairInt(1,-1) }),
		fourPtNeighborShifts({ PairInt(1,0),PairInt(0,1),PairInt(-1,0),PairInt(0,-1) }) {}

	GridInfo(const GridInfo& in) : nx(in.nx), xl(in.xl), xr(in.xr), ny(in.ny), yl(in.yl), yr(in.yr), 
		hx(in.hx), hy(in.hy), 
		eightPtNeighborShifts(in.eightPtNeighborShifts),
		fourPtNeighborShifts(in.fourPtNeighborShifts) {}

	double getXValue(int index) const { return xl + hx * (index % ny); } // Recall that the const qualifier here indicates that 'this' is constant.
	double getYValue(int index) const { return yl + hy * (index / ny); }
	int getXIndex(int index) const { return index % ny; }
	int getYIndex(int index) const { return index / ny; }
	int gridLength() const { return nx * ny; }
	bool inBoundary(int index) const { return ((getXIndex(index) == 0) || (getXIndex(index) == nx - 1) || (getYIndex(index) == 0) || (getYIndex(index) == ny - 1)); }
	bool invalidIndex(int index) const { return ((index < 0) || (index >= gridLength())); }
};



// ***************** Speed and Solver Structures **************



struct SpeedParams {
	int switchKey;
	// for drift field b = [-ax - by, abx -y]
	// solution is u = ax^2 + y^2
	double a;
	double b;
	SpeedParams(int key, double aIn, double bIn) :switchKey(key), a(aIn), b(bIn) {}
	SpeedParams(const SpeedParams& spIn) :switchKey(spIn.switchKey), a(spIn.a), b(spIn.b) {}
};

// This contains all drift and slowness functions in a single structure of function pointers.
// Here are the forward declarations of the user options for the speeds.
// External functions can also be set.
PairDoub masterDrift(const PairDoub& z, const void* sp);
MatrixDoub masterGradDrift(const PairDoub& z, const void* sp);
TensorDoub masterHessDrift(const PairDoub& z, const void* sp);
double masterSol(const PairDoub& z, const void* sp);
double masterSlow(const PairDoub& z, const PairDoub& v, const void* sp);
PairDoub masterGradZ(const PairDoub& z, const PairDoub& v, const void* sp);
MatrixDoub masterHessZ(const PairDoub& z, const PairDoub& v, const void* sp);
PairDoub masterGradV(const PairDoub& z, const PairDoub& v, const void* sp);
MatrixDoub masterHessV(const PairDoub& z, const PairDoub& v, const void* sp);
PairDoub masterGradSol(const PairDoub& z, const void* sp);
double masterLaplacianSol(const PairDoub& z, const void* sp);
std::string masterStr(const SpeedParams &spdPrms);
MatrixDoub masterDriftLinearization(const void* sp);
double masterDivDrift(const PairDoub& z, const void* sp);
// The reason these functions all take const void* parameters is for ease of inputting external functions
// However, if we're already 
struct SpeedInfo {
	std::string str;
	bool quasi;
	const SpeedParams& sp;
	double (*slow) (const PairDoub& z, const PairDoub& v, const void* sp);
	PairDoub(*slowGradZ) (const PairDoub& z, const PairDoub& v, const void* sp);
	MatrixDoub(*slowHessZ) (const PairDoub& z, const PairDoub& v, const void* sp);
	PairDoub(*slowGradV) (const PairDoub& z, const PairDoub& v, const void* sp);
	MatrixDoub(*slowHessV) (const PairDoub& z, const PairDoub& v, const void* sp);
	PairDoub(*drift)(const PairDoub& z, const void* sp);
	MatrixDoub(*gradDrift)(const PairDoub& z, const void* sp);
	TensorDoub(*hessDrift)(const PairDoub& z, const void* sp);
	double (*solution) (const PairDoub& z, const void* sp);
	PairDoub(*gradSolution)(const PairDoub& z, const void* sp);
	double (*laplaceSolution) (const PairDoub& z, const void* sp);
	MatrixDoub(*driftLinearization)(const void* sp);
	double(*divDrift)(const PairDoub& z, const void* sp);

	SpeedInfo(const SpeedParams& params, bool quasiIn) : quasi(quasiIn),sp(params), drift(masterDrift),
		solution(masterSol), slow(masterSlow), laplaceSolution(masterLaplacianSol),
		slowGradZ(masterGradZ), slowHessZ(masterHessZ), slowGradV(masterGradV),
		slowHessV(masterHessV),gradSolution(masterGradSol), 
		gradDrift(masterGradDrift), str(masterStr(params)), hessDrift(masterHessDrift), 
		driftLinearization(masterDriftLinearization),
		divDrift(masterDivDrift) {}
};

// This is a struct that gets passed to any of the various solvers.
// Solvers may not use all fields.
struct SolverParams2Pt {
	PairDoub z1, z2, z3;
	double u2, u3;
	PairDoub u2del, u3del;
	const SpeedInfo& si;
};
struct SolverParams1Pt {
	PairDoub z1, z2;
	double u2;
	PairDoub u2del;
	const SpeedInfo& si;
};
struct Flags {
	int updateType{ 0 };
	double a0{ 0 };
	double a1{ 0 };
	double lambda{ 0 };

};
struct UpdateInfo {
	int ind1, ind2, ind3;
	PairDoub z1, z2, z3;
	double u1, u2, u3;
	PairDoub u1del, u2del, u3del;
	double lam, a0, a1;
	int lastUpdateType;
	int rank = 0;
	void updateFlags(const Flags& flags_in) {
		lastUpdateType = flags_in.updateType;
		lam = flags_in.lambda;
		a0 = flags_in.a0;
		a1 = flags_in.a1;
	}
	void updateInd(int ind2_in, int ind3_in) {
		ind2 = ind2_in;
		ind3 = ind3_in;
	}
};

/*
	Struct containing miscellaneous options.
*/
struct  RunOptions {
	bool quasi;
	bool verbose;
	bool debugMode;
	bool twoPtUpdateMinDist;
	bool failSafe;
	bool fakeFilter;
	int initType;
	bool z3one = false;

};

struct SummaryStats {
	double err_sup;
	double err_rms;
	double err_grad_sup;
	double err_grad_rms;
	double run_time;
	int num_acc;
	int num_acc_1pt;
	int num_acc_2pt;
	SummaryStats() : err_sup(INFINITY), err_rms(INFINITY), err_grad_sup(INFINITY), err_grad_rms(INFINITY),
		run_time(INFINITY), num_acc(0), num_acc_1pt(0), num_acc_2pt(0) {}

};


double gauss();
#endif
