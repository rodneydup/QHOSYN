#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <stack>

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
  double qhoBasis(int n, double x) {
    return (pow(M_PI, 0.25) / sqrt(pow(2, n) * factorial(n))) * exp(pow(x, 2) / -2.0) *
           std::hermite(n, x);
  };

  void calcEigenBasis() {
    for (int i = 0; i < dim; i++)
      qhoBasisApprox.push_back(std::make_unique<SampleFunction>(
        [=](double x) -> double { return qhoBasis(i, x); }, -15, 15, 2000));
  }

  double eigenbasis(int n, float x) { return qhoBasisApprox[n]->lookup(x); }

  template <typename T>
  T eigenvalues(T x) {
    return 0.5 + x;
  }

  std::vector<std::unique_ptr<SampleFunction>> qhoBasisApprox;
  int dim;  // dimensions of the hilbert space (effectively the number of superposed energy states)
  double (*v)(double);  // hamiltonian operator
};

// The QHO wavefunction
class WaveFunction {
 public:
  WaveFunction(HilbertSpace *hs, double initWaveFunc(double), std::vector<double> coeff = {}) {
    hilbertSpace = hs;
    if (!coeff.empty()) {
      coefficients = coeff;
    } else {
      coefficients = orthogonalBasisProjection(initWaveFunc);
    }

    // find the normalize value for this wavefunction
    auto ptr = [=](double x) {
      std::complex<double> sum;
      for (int n = 0; n < hilbertSpace->dim; n++)
        sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, 0);
      return pow(abs(sum), 2);
    };
    gsl_function_pp<decltype(ptr)> Fp(ptr);
    gsl_function *normFunc = static_cast<gsl_function *>(&Fp);
    gsl_integration_cquad_workspace *normW = gsl_integration_cquad_workspace_alloc(1000);
    gsl_integration_cquad(normFunc, -1, 1, 1.49e-08, 1.49e-08, normW, &integrationResult, NULL,
                          NULL);
    normF = sqrt(integrationResult);
    gsl_integration_cquad_workspace_free(normW);
  }

  std::vector<double> orthogonalBasisProjection(double waveFunc(double)) {
    std::vector<double> tempCoeff;
    tempCoeff.resize(hilbertSpace->dim);

    gsl_integration_cquad_workspace *orthoW = gsl_integration_cquad_workspace_alloc(1000);
    for (int i = 0; i < tempCoeff.size(); i++) {
      auto ptr = [=](double x) {
        return (std::conj(hilbertSpace->eigenbasis(i, x)) * waveFunc(x)).real();
      };
      gsl_function_pp<decltype(ptr)> Fp(ptr);
      gsl_function *orthoFunc = static_cast<gsl_function *>(&Fp);
      gsl_integration_cquad(orthoFunc, -1, 1, 1.49e-08, 1.49e-08, orthoW, &integrationResult, NULL,
                            NULL);
      tempCoeff[i] = integrationResult;
    }
    gsl_integration_cquad_workspace_free(orthoW);
    std::cout << tempCoeff[1] << std::endl;
    return tempCoeff;
  }

  std::complex<double> phaseFactor(double n, double t) {
    return exp(-1i * hilbertSpace->eigenvalues(n) * t);
  }

  std::complex<double> normalize(std::complex<double> waveFunc(double)) {
    auto ptr = [=](double x) { return pow(abs(waveFunc(x)), 2); };
    gsl_function_pp<decltype(ptr)> Fp(ptr);
    gsl_function *normFunc = static_cast<gsl_function *>(&Fp);
    gsl_integration_cquad_workspace *normW = gsl_integration_cquad_workspace_alloc(1000);
    gsl_integration_cquad(normFunc, -INFINITY, INFINITY, 1.49e-08, 1.49e-08, normW,
                          &integrationResult, NULL, NULL);
    gsl_integration_cquad_workspace_free(normW);
    return sqrt(integrationResult);
  }

  std::complex<double> evaluate(double x, double t) {
    std::complex<double> sum = (0, 0);
    for (int n = 0; n < hilbertSpace->dim; n++)
      sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, t);
    // std::cout << coefficients[0] << std::endl;
    return sum / normF;
  }

  double integrationResult;
  double normF;
  HilbertSpace *hilbertSpace;
  std::vector<double> coefficients;
};