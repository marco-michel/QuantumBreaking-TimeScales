#include "quantumBreakingHamiltonians.h"
#include <cmath>



double quantumBreakingHamiltonians::symBreak(int a, int b)
{
	double coef = std::fmod(std::sqrt(2) * std::pow(a, 3) + std::sqrt(7) * std::pow(b, 5), 1.0);
	if (coef < 0.5)
		coef -= 1.0;
	return coef;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int periodicHamiltonian::momentumMode(int index)
{
	int shiftedIndex;

	if (nbSpecies > 1 && index > 2 * K + 1)
		shiftedIndex = index % (2 * K);
	else
		shiftedIndex = index;

	if (index % 2 == 0)
		return (int) -shiftedIndex / 2.0;
	else
		return (int)ceil(shiftedIndex / 2.0);
}


int periodicHamiltonian::indexMode(int mode) //only gives the first species
{
	if (mode == 0)
		return 0;
	else if (mode < 0)
		return -2 * mode;
	else
		return 2 * mode - 1;
}

periodicHamiltonian::periodicHamiltonian(int n, int q, int k, int c, double alpha, double cm)
{
	N = n; nbSpecies = q; C = c; Alpha = alpha; numModes = 1 + 2 * k; K = k; Cm = cm;
}


std::vector<opTerm> periodicHamiltonian::createHamiltonian()
{
	std::vector<opTerm> ret;
	opTerm tmp;

	std::vector<opTerm> ham;


	for (int q = 0; q < nbSpecies; q++)
	{

		int offset = q * 2 * K; //modes for species i are shifted i * numModes/species.
		//Mapping mod numbers: 
		//odd -> ceil(i/2), even -> -i/2

		for (int i = 1; i != numModes; i++)
		{
			ham.push_back(createNumberOperator(i + offset, std::pow(momentumMode(i), 2)));
		}

		for (int i = -K; i <= K; i++)
		{
			for (int j = -K; j <= K; j++)
			{
				for (int k = -K; k <= K; k++)
				{
					for (int l = -K; l <= K; l++)
					{
						if (i + j - k - l == 0)
						{
							int im, jm, km, lm;
							im = (i == 0) ? 0 : indexMode(i) + offset;
							jm = (j == 0) ? 0 : indexMode(j) + offset;
							km = (k == 0) ? 0 : indexMode(k) + offset;
							lm = (l == 0) ? 0 : indexMode(l) + offset;
							ham.push_back(createFourPointB(im, jm, km , lm, 0, 0, 0, 0, -Alpha / 4.0));
						}
					}
				}
			}
		}


		for (int qq = 0; qq < nbSpecies; qq++)
		{
			int offset2 = qq * 2 * K;
			if (offset != offset2)
			{
				double symBCoef = Cm * symBreak(offset, offset2);
				for (int ll = 1; ll <= 2*K; ll++)
				{
					ham.push_back(linInteraction(ll + offset, ll + offset2, 0, 0, true, symBCoef));
					ham.push_back(linInteraction(ll + offset, ll + offset2, 0, 0, false, symBCoef));
				}
			}
		}
	}


	hamiltonOperator = ham;
	hamiltonianString = ham;
	return ham;
}


basisVector periodicHamiltonian::createInitState(int state, int substate)
{
	basisVector ret;

	ret.length = 1 + 2 * K * nbSpecies;
	ret.e = new int[ret.length];

	if (state == 1)
	{
		ret.e[0] = N;
		for (int i = 1; i != ret.length; i++)
			ret.e[i] = 0;
	}
	else
	{
		std::cerr << "no initial state selected" << std::endl;
		exit(14);
	}

	return ret;
}


void periodicHamiltonian::createObservables(smatrix**& out, basicBasis* basis)
{
	int nbModes = basis->numberModes;

	int nbObs =  (1 + 2 * K * nbSpecies);

	out = new smatrix * [nbObs];

	for (int i = 0; i != nbModes; i++)
	{
		std::vector<opTerm> operators;
		operators.push_back(createNumberOperator(i, 1.0));
		smatrix* matrix;
		matrix = createMatrix(operators, basis);
		out[i] = matrix;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////

simplifiedHamiltonianQuadInt::simplifiedHamiltonianQuadInt(int n, int q, int c, double alpha, double cm, double cgap, bool definiteSign, bool doubleSum)
{
	N = n; nbSpecies = q; C = c; Alpha = alpha; numModes = 2; Cm = cm; K = 1; Cgap = cgap; DefiniteSign = definiteSign; DoubleSum = doubleSum;
}

basisVector simplifiedHamiltonianQuadInt::createInitState(int state, int substate)
{
	basisVector ret;

	ret.length = 1 + nbSpecies;
	ret.e = new int[ret.length];
	if (state == 1)
	{
		ret.e[0] = N;
		for (int i = 1; i != ret.length; i++)
			ret.e[i] = 0;
	}
	else
	{
		std::cerr << "no initial state selected" << std::endl;
		exit(14);
	}
	return ret;
}

std::vector<opTerm> simplifiedHamiltonianQuadInt::createHamiltonian()
{
	std::vector<opTerm> ret;
	opTerm tmp;

	std::vector<opTerm> ham;


	for (int q = 0; q < nbSpecies; q++)
	{

		int offset = q;

		ham.push_back(createNumberOperator(1 + offset, Cgap));


		ham.push_back(createFourPointB(0, 0, 1 + offset, 1 + offset, 0, 0, C, C, -Alpha / 4.0));
		ham.push_back(createFourPointB(1 + offset, 1 + offset, 0, 0, C, C, 0, 0, -Alpha / 4.0));
		ham.push_back(createFourPointA(0, 0, 1 + offset, 1 + offset, 0, 0, C, C, -Alpha / 2.0));

		for (int qq = 0; qq < nbSpecies; qq++)
		{
			int offset2 = qq;
			double symBCoef;
			if (offset != offset2 && (DoubleSum || offset < offset2))
			{
				//!! Added factor of 1/2 here
				if(! DefiniteSign)
				   symBCoef = Cm/2.0 * symBreak(offset, offset2);
				else
				   symBCoef = Cm/2.0 * std::abs(symBreak(offset, offset2));

				ham.push_back(createFourPointB(1 + offset, 1 + offset, 1 + offset2, 1 + offset2, C, C, C, C, symBCoef));
				ham.push_back(createFourPointB(1 + offset2, 1 + offset2, 1 + offset, 1 + offset, C, C, C, C, symBCoef));
			}
		}
	}


	hamiltonOperator = ham;
	hamiltonianString = ham;
	return ham;
}

void simplifiedHamiltonianQuadInt::createObservables(smatrix**& out, basicBasis* basis)
{
	int nbModes = basis->numberModes;

	int nbObs =  nbModes;

	out = new smatrix * [nbObs];

	for (int i = 0; i != nbModes; i++)
	{
		std::vector<opTerm> operators;
		operators.push_back(createNumberOperator(i, 1.0));
		smatrix* matrix;
		matrix = createMatrix(operators, basis);
		out[i] = matrix;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
