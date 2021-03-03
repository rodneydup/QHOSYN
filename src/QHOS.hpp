#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <gsl/gsl_integration.h>

#include <cmath>
#include <complex>
#include <iostream>

// for master branch
// #include "al/core.hpp"

// for devel branch
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
// QHOS
#include "utility.hpp"

using namespace al;

class QHOS : public App {
 public:
  /**
   * @brief Initilialize the synth interface.
   */
  virtual void onInit() override;

  /**
   * @brief Run once on starup.
   */
  virtual void onCreate() override;

  /**
   * @brief Audio rate processing of synth.
   */
  virtual void onSound(al::AudioIOData &io) override;

  virtual void onAnimate(double dt) override;

  /**
   * @brief Draw rate processing of synth interface.
   */
  virtual void onDraw(al::Graphics &g) override;

  // Variables
  Mesh plot;

  std::vector<double> posValues = linspace(-5, 5, 300);
  std::vector<double> amplitudeValues;
  HilbertSpace hs{2};

  WaveFunction psi{&hs, [](double x) {
                     if (abs(x) < 5)
                       return 1.0;
                     else
                       return 0.0;
                   }};
};

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

  // pi^1/4 / sqrt(2^n * n!)  * e^(-x^2/2) * Hn
  double qhoBasis(int n, int x) {
    return 1.0 / sqrt(pow(2, n) * factorial(n)) * pow(M_PI, 1 / 4) * exp(pow(-x, 2) / 2) *
           std::hermite(n, x);
  };

  void calcEigenBasis() {
    for (int i = 0; i < dim; i++)
      qhoBasisApprox.push_back(std::make_unique<SampleFunction>(
        [=](int x) -> double { return qhoBasis(i, x); }, -15, 15, 2000));
  }

  double eigenbasis(int n, int x) { return qhoBasisApprox[n]->lookup(x); }

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
    if (!coeff.empty())
      coefficients = coeff;
    else {
      coefficients = orthogonalBasisProjection(initWaveFunc);
    }
    // normF = normalize(prevaluate(x,0));
    normF = normalize([this](double x) -> double {
      double sum;
      for (int n = 0; n < hilbertSpace->dim; n++)
        sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, 0);
    });
    // normalize(prevaluate(x, 0))
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
      gsl_integration_cquad(orthoFunc, -INFINITY, INFINITY, 1.49e-08, 1.49e-08, orthoW,
                            &integrationResult, NULL, NULL);
      tempCoeff[i] = integrationResult;
    }
    gsl_integration_cquad_workspace_free(orthoW);
    return tempCoeff;
  }

  double phaseFactor(int n, int t) { return exp(sqrt(-1) * hilbertSpace->eigenvalues(n) * t); }
  template <typename T>
  double normalize(T waveFunc(double)) {
    auto ptr = [=](double x) { return pow(abs(waveFunc(x)), 2); };
    gsl_function_pp<decltype(ptr)> Fp(ptr);
    gsl_function *normFunc = static_cast<gsl_function *>(&Fp);
    gsl_integration_cquad_workspace *normW = gsl_integration_cquad_workspace_alloc(1000);
    gsl_integration_cquad(normFunc, -INFINITY, INFINITY, 1.49e-08, 1.49e-08, normW,
                          &integrationResult, NULL, NULL);
    gsl_integration_cquad_workspace_free(normW);
    return sqrt(integrationResult);
  }

  // double prevaluate(int x, int t) {
  //   double sum;
  //   for (int n = 0; n < hilbertSpace->dim; n++)
  //     sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, t);
  //   return sum;
  // }

  double evaluate(int x, int t) {
    double sum;
    for (int n = 0; n < hilbertSpace->dim; n++)
      sum += coefficients[n] * hilbertSpace->eigenbasis(n, x) * phaseFactor(n, t);
    return sum / normF;
  }

  double integrationResult;
  double normF;
  HilbertSpace *hilbertSpace;
  std::vector<double> coefficients;
  // double (*evaluate)(double x, double t);
};