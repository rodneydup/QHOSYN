#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <stack>

#ifdef __APPLE__
#include <boost/math/special_functions/hermite.hpp>
#endif

#include "utility.hpp"

using namespace std::complex_literals;

// The state space of the wave function
class HilbertSpace {
 public:
  HilbertSpace(int dimensions, double hamiltonianPotential(double)) {
    v = hamiltonianPotential;
    dim = dimensions;
    calcEigenBasis();
  }

  HilbertSpace(int dimensions) {
    dim = dimensions;
    calcEigenBasis();
  }    
  // pi^1/4 / sqrt(2^n * n!)  * e^(-(x^2/2)) * Hn
#ifdef __APPLE__
  double qhoBasis(int n, double x) {
    return (pow(M_PI, 0.25) / sqrt(pow(2, n) * factorial(n))) * exp(pow(x, 2) / -2.0) *
           boost::math::hermite(n, x);
  };
#else
  double qhoBasis(int n, double x) {
    return (pow(M_PI, 0.25) / sqrt(pow(2, n) * factorial(n))) * exp(pow(x, 2) / -2.0) *
           std::hermite(n, x);
  };
#endif

  void calcEigenBasis() {
    for (int i = 0; i < dim; i++)
      qhoBasisApprox.push_back(std::make_unique<SampleFunction>(
        [=](double x) -> double { return qhoBasis(i, x); }, -10, 10, 1024));
  }

  double eigenbasis(int n, float x) { return qhoBasisApprox[n]->lookup(x); }

  template <typename T>
  T eigenvalues(T x) {
    return 0.5 + x;
  }

  std::vector<std::unique_ptr<SampleFunction>> qhoBasisApprox;
  int dim;              // dimensions of the hilbert space (effectively the number of superposed
                        // energy states)
  double (*v)(double);  // hamiltonian operator
};

// The QHO wavefunction
class WaveFunction {
 public:
  WaveFunction(
    HilbertSpace *hs, int dimensions,
    std::function<double(double)> initWaveFunc = [](double x) { return 1.0; },
    std::vector<double> coeff = {}, bool project = 1) {
    newHilbertSpace(hs, dimensions, initWaveFunc, coeff, project);
  }

  // calculate phase
  std::complex<double> phaseFactor(double n, double t) {
    return exp(-1i * hilbertSpace->eigenvalues(n) * t);
  }

  // Lookup value from wavefuntion for time t at position value x, summing over
  // eigenstates
  std::complex<double> wavefunc(double x, double t) {
    std::complex<double> sum = (0, 0);
    for (int n = 0; n < eigenStates; n++)
      sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, t);
    return sum;
  }

  // calculate normalize value
  void normalize() {
    auto ptr = [=](double x) { return pow(abs(wavefunc(x, 0)), 2); };
    gsl_function_pp<decltype(ptr)> Fp(ptr);
    gsl_function *normFunc = static_cast<gsl_function *>(&Fp);
    gsl_integration_cquad_workspace *normW = gsl_integration_cquad_workspace_alloc(1000);
    gsl_integration_cquad(normFunc, -1, 1, 1.49e-08, 1.49e-08, normW, &integrationResult, NULL,
                          NULL);
    gsl_integration_cquad_workspace_free(normW);
    normF = sqrt(integrationResult);
  }

  // Return wavefunc normalized
  std::complex<double> evaluate(double x, double t) { return wavefunc(x, t) / normF; }

  // get the coefficients from an orthogonal basis projection of the initial function
  std::vector<double> orthogonalBasisProjection(std::function<double(double x)> waveFunc) {
    std::vector<double> tempCoeff;
    tempCoeff.resize(eigenStates);

    gsl_integration_cquad_workspace *orthoW = gsl_integration_cquad_workspace_alloc(1000);
    for (int i = 0; i < tempCoeff.size(); i++) {
      auto ptr = [=](double x) {
        return (std::conj(hilbertSpace->eigenbasis(i, x)) * waveFunc(x)).real();
      };
      gsl_function_pp<decltype(ptr)> Fp(ptr);
      gsl_function *orthoFunc = static_cast<gsl_function *>(&Fp);
      gsl_integration_cquad(orthoFunc, -10, 10, 1.49e-08, 1.49e-08, orthoW, &integrationResult,
                            NULL, NULL);
      tempCoeff[i] = integrationResult;
    }
    gsl_integration_cquad_workspace_free(orthoW);
    return tempCoeff;
  }

  // instantiate a new hilbert space
  void newHilbertSpace(
    HilbertSpace *hs, int dimensions,
    std::function<double(double x)> initWaveFunc = [](double x) { return 1.0; },
    std::vector<double> coeff = {}, bool project = 1) {
    eigenStates = dimensions;
    hilbertSpace = hs;
    if (!coeff.empty()) {
      coefficients = coeff;
    } else {
      coefficients.resize(dimensions);
      if (project) {
        coefficients = orthogonalBasisProjection(initWaveFunc);
      } else {
        for (int i = 0; i < eigenStates; i++) coefficients[i] = initWaveFunc(i);
      }
    }
    normalize();
  }

  // collapse the wavefunction, would need to find a numerical solution to make this work
  // in real-time and continue evolving.
  void collapse() {}

  int eigenStates;
  double integrationResult;
  double normF;
  HilbertSpace *hilbertSpace;
  std::vector<double> coefficients;
};