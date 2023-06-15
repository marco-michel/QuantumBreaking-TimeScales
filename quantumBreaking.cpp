#include <unordered_map>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <chrono>

#undef I
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <mkl_types.h>
#include <mkl.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif

#include "matrixDataTypes.h"
#include "krylovTimeEvolver.h"
#include "krylovHelper.h"
#include "Basis.h"
#include "hamiltonian.h"
#include "quantumBreakingHamiltonians.h"

//INPUT:
//    int N0; int Nm; int K;
//  double C0; double Cm
//  double maxT; double samplingStep; 
//   double tol; int m; int numThreads; int DeltaN; int capacity;  bool fastIntegration
//  in total, argc can take up to 13 arguments
int main(int argc, char* argv[])
{

    //Determine parameters
    int N;
    int K;
    double Alpha, Cm, Cgap;
    int capacity;
    int Ham;
    int Q;

    double maxT; double samplingStep;
    double tol; int m; int numThreads;
    bool fastIntegration;
    bool progressBar;

    po::options_description desc("Allowed options");
    po::variables_map vm;

    try {

        desc.add_options()
            ("help", "produce help message")
            ("N", po::value<int>(&N)->default_value(20), "Number of particles")
            ("K", po::value<int>(&K)->default_value(1), "Highest momentum mode")
            ("Q", po::value<int>(&Q)->default_value(3), "Number of species")
            ("Ham", po::value<int>(&Ham)->default_value(2), "Hamiltonian")
            ("Alpha", po::value<double>(&Alpha)->default_value(0.025), "Coupling")
            ("Cm", po::value<double>(&Cm)->default_value(0.025), "Interspecies coupling")
            ("Cgap", po::value<double>(&Cgap)->default_value(1.0), "Interspecies coupling")
            ("maxT", po::value<double>(&maxT)->default_value(10), "Simulation-time")
            ("samplingStep", po::value<double>(&samplingStep)->default_value(0.01), "Time interval of sampling")
            ("tol", po::value<double>(&tol)->default_value(1.0e-8), "Numerical tolerance")
            ("m", po::value<int>(&m)->default_value(40), "Dimension of Krylov-Space")
            ("threads", po::value<int>(&numThreads)->default_value(2), "Number of OpenMP Threads for Intel MKL")
            ("capacity", po::value<int>(&capacity)->default_value(20), "Capacity")
            ("fastIntegration", po::value<bool>(&fastIntegration)->default_value(true), "Use faster and less accurate integration")
            ("progressBar", po::value<bool>(&progressBar)->default_value(false), "Display a progress bar during time evolution")
            ;

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
    }
    catch (const boost::program_options::required_option& e) {
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        else {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    //Create Hamiltonian
    quantumBreakingHamiltonians* ham;

    if (Ham == 1)
        ham = new periodicHamiltonian(N, Q, K, capacity, Alpha, Cm);
    else if (Ham == 2)
        ham = new simplifiedHamiltonianQuadInt(N, Q, capacity, Alpha, Cm, Cgap, false, false);
    else
    {
        std::cerr << "No valid Hamiltonian selected" << std::endl;
        return 1;
    }


    int nbObservables;
    int nbModes;

    if (Ham == 1) //periodic
        nbModes = 2 * K * Q + 1 ;
    else if (Ham == 2) //simplified
    {
        if (K != 1)
        {
            std::cerr << "Only one mode allowed for the selected Hamiltonian" << std::endl;
            return 1;
        }
        nbModes = Q + 1;
    }



    mkl_set_num_threads(numThreads);
    //end determine parameters


    std::cout << "QuantumBreaking" << std::endl;

    //Create basis

    basis* bas;
    if (Ham==2) 
	bas = new basis(N, nbModes, nbModes-1, capacity,2);
    else
	bas = new basis(N, nbModes, nbModes-1, capacity);
    std::cout << "Creating basis with " << bas->numberElements << " elements..." << std::endl;





    ham->createHamiltonian();
    //std::cout << ham->toString() << std::endl;

    //Create Hamiltonian matrix
    std::cout << "Creating Hamiltonian matrix..." << std::endl;
    smatrix* hamMatrix;
    ham->createHamiltonMatrix(hamMatrix, bas);



    //Create initial state
    basisVector init = ham->createInitState(1, 0);
    std::complex<double>* vec = new std::complex<double>[bas->numberElements];
    //find init state in hash table
    int entry = bas->hashTable.find(init)->second;
    vec[entry].real(1.0);

    //Create matrices for observables
    std::cout << "Creating observables..." << std::endl;
    smatrix* complexObs = 0;
    smatrix** observables;
    smatrix** observablesNumOp;
    nbObservables =  nbModes;

    ham->createObservables(observablesNumOp, bas);
    std::vector<std::unique_ptr<krylovBasicObservable>> awu;

    for (int i = 0; i != bas->numberModes; i++)
    {
        awu.push_back(std::make_unique<krylovSpMatrixObservable>("mode" + std::to_string(i), observablesNumOp[i]));
    }

    delete[] observablesNumOp;


    //Start of actual time evolution   
    std::cout << "Starting time evolution..." << std::endl;
    std::complex<double> imaginaryMinus;
    imaginaryMinus.imag(-1);

    auto begin = std::chrono::high_resolution_clock::now();

    krylovTimeEvolver timeEvolver(maxT, bas->numberElements, vec, samplingStep, tol, m, std::move(awu), hamMatrix, imaginaryMinus, true, fastIntegration, progressBar);
    krylovReturn* results = timeEvolver.timeEvolve();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    std::cout << "Time: " << elapsed.count() * 1e-9 << " sec" << std::endl;

    std::cout << "------------------------------------------------------\n" << std::endl;
    std::cout << "Number of steps: " << results->n_steps << "    error: " << results->err << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    //End of actual time evolution 




    //Create string with information about input parameters
    std::cout << "Writing results to file..." << std::endl;

    parameter_list parameters;

    parameters.push_back(paraPush("Ham", Ham));
    parameters.push_back(paraPush("N", N));
    parameters.push_back(paraPush("Q", Q));
    parameters.push_back(paraPush("K", K));
    parameters.push_back(paraPush("Alpha", Alpha));
    parameters.push_back(paraPush("Cm", Cm));
    parameters.push_back(paraPush("Cgap", Cgap));
    parameters.push_back(paraPush("C", capacity));
    parameters.push_back(paraPush("maxT", maxT));
    parameters.push_back(paraPush("tol", tol));
    parameters.push_back(paraPush("samplingStep", samplingStep));
    parameters.push_back(paraPush("m", m));
    parameters.push_back(paraPush("fastIntegration", fastIntegration));



    observable_list obsus; //correct order is still important here for correct naming
    for (int i = 0; i != nbModes; i++)
        obsus.push_back(obsPush("mode" + std::to_string(i), obsType::SPARSE_MATRIX_TYPE_OBS));


    outputHelper fileHelper(results, parameters, obsus, "TimeScalesQB");
    fileHelper.saveResult();

    //clean up

    delete results; delete hamMatrix; delete[] vec; 
    delete ham;

    return 0;
}
