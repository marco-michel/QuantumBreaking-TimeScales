#pragma once

#include <vector>
#include <complex>

#include "hamiltonian.h"
#include "Basis.h"





class quantumBreakingHamiltonians : public Hamiltonian
{
public:
    int C; int N; double Alpha; int nbSpecies; int numModes; int K; double Cm; double Cgap;
    std::vector<opTerm> hamiltonianString;

    virtual void createObservables(smatrix**& out, basicBasis* basis) = 0;
    virtual std::vector<opTerm> createHamiltonian() = 0;
    virtual basisVector createInitState(int state, int substate) = 0;

protected:
    double symBreak(int a, int b);

};


class periodicHamiltonian : public  quantumBreakingHamiltonians
{
public:
    periodicHamiltonian(int n, int q, int k, int c, double alpha, double cm);
    basisVector createInitState(int state, int substate);
    std::vector<opTerm> createHamiltonian();

    void createObservables(smatrix**& out, basicBasis* basis);
private:

    int momentumMode(int index);
    int indexMode(int mode);
};




class simplifiedHamiltonianQuadInt : public quantumBreakingHamiltonians
{
public:
    bool DefiniteSign, DoubleSum;
    simplifiedHamiltonianQuadInt(int n, int q, int c, double alpha, double cm, double cgap, bool definiteSign=false, bool doubleSum=true);
    basisVector createInitState(int state, int substate);
    std::vector<opTerm> createHamiltonian();

    void createObservables(smatrix**& out, basicBasis* basis);
};

