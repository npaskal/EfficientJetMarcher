#ifndef _UPDATER_H
#define _UPDATER_H
#include "structs.h"

double gauss();

class Updater {
public:
	virtual ~Updater() {}
	virtual double onePointUpdateValue(const SolverParams1Pt& params, PairDoub& u1grad, Flags& flags) const =0;
	virtual double twoPointUpdateValue(const SolverParams2Pt& params, PairDoub& u1grad, bool& interiorMinimizer, Flags& flags) const =0;
};

class LinearChar : public Updater {
private:
	double hybridSolver(double a0, double a1, const SolverParams2Pt& p, bool& interiorSolution) const;
public:
	virtual double F1(const SolverParams1Pt& p) const =0;
	virtual double F2(double r, const SolverParams2Pt& p) const =0;
	virtual double F2p(double r, const SolverParams2Pt& p) const =0;
	double onePointUpdateValue(const SolverParams1Pt& p, PairDoub& u1grad, Flags& flags) const override;
	double twoPointUpdateValue(const SolverParams2Pt& p, PairDoub& u1grad, bool& interiorMinimizer, Flags& flags) const override;
};

class EndpointLinear : public LinearChar {
public:
	double F1(const SolverParams1Pt& p) const override;
	double F2(double r, const SolverParams2Pt& p) const override;
	double F2p(double r, const SolverParams2Pt& p) const override;
};

class MidpointLinear : public LinearChar {
public:
	double F1(const SolverParams1Pt& p) const override;
	double F2(double r, const SolverParams2Pt& p) const override;
	double F2p(double r, const SolverParams2Pt& p) const override;
};

class Hermite : public Updater {
private:
	// These fields set whether the action quadrature evaluation uses a larger
	// number of points than did the corresponding minimization.
	bool higherEval_1pt = true;
	bool higherEval_2pt = true;
	unsigned evalPts_1pt =13;
	unsigned evalPts_2pt = 13;

	std::vector<double> NewtonSolver2D(double s0, double t0,
		const SolverParams1Pt& p, bool &converges) const;

	std::vector<double> NewtonSolver3D(double a0, double a1, double lam,
		const SolverParams2Pt& p, bool& converges) const;

	// Here are the respective 1pt and 2pt update formulas, as well as their
	// gradients and Hessians.
	double F1(double a0, double a1, const SolverParams1Pt& p) const;
	double F1_higher(double a0, double a1, const SolverParams1Pt& p) const;
	std::vector<double> F1p(double a0, double a1,
		const SolverParams1Pt& p) const;
	std::vector<double> F1pp(double a0, double a1,
		const SolverParams1Pt& p) const;
	double F2(double a0, double a1, double r, const SolverParams2Pt& p) const;
	double F2_higher(double a0, double a1, double r, 
		const SolverParams2Pt& p) const;
	std::vector<double> F2p(double a0, double a1, double r, 
		const SolverParams2Pt& p) const;
	std::vector<double> F2pp(double a0, double a1, double r,
		const SolverParams2Pt& p) const;

public:
	double onePointUpdateValue(const SolverParams1Pt& p, PairDoub& u1grad, 
		Flags& flags) const override;
	double twoPointUpdateValue(const SolverParams2Pt& p, PairDoub& u1grad, 
		bool& interiorMinimizer, Flags& flags) const override;
	void checkFormulas(const SpeedInfo& speeds, const SpeedParams& sp) const;
};












#endif