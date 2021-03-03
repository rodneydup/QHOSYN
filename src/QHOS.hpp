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
};

// helper functions
int factorial(int x) {
  int result = 1;
  for (int i = 1; i <= x; i++) result = result * i;
  return result;
}

template <typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in) {
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) {
    return linspaced;
  }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }

  double delta = (end - start) / (num - 1);

  for (int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);  // I want to ensure that start and end
                             // are exactly the same as the input
  return linspaced;
}

// class that discretizes (samples) a function and stores results for lookup later.
class SampleFunction {
 public:
  template <typename T>
  SampleFunction(T f, int min, int max, int numSamples) {
    minX = min;
    maxX = max;
    rangeX = max - min;
    numSamplesX = numSamples;

    std::vector<double> domain = linspace(minX, maxX, numSamplesX);
    image.resize(domain.size());
    for (int x = 0; x < domain.size(); x++) image[x] = (f(domain[x]));
  }

  double lookup(int x) {
    x = std::min(std::max(x, minX), maxX - 1);
    x -= minX;
    x = std::round(x * (numSamplesX / rangeX));
    return image[int(x)];
  }

 private:
  int minX, maxX, rangeX, numSamplesX;
  std::vector<double> image;
};

template <typename F>
class gsl_function_pp : public gsl_function {
 public:
  gsl_function_pp(const F &func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params = this;
  }

 private:
  const F &_func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp *>(params)->_func(x);
  }
};

// The state space of the wave function
class HilbertSpace {
 public:
  HilbertSpace() {}
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
    0.5 + x;
  }

  std::vector<std::unique_ptr<SampleFunction>> qhoBasisApprox;
  int dim;  // dimensions of the hilbert space (effectively the number of superposed energy states)
  double (*v)(double);  // hamiltonian operator
};

// The QHO wavefunction
class waveFunction {
 public:
  waveFunction(HilbertSpace *hs, double initWaveFunc(double), std::vector<double> coeff) {
    hilbertSpace = hs;
    if (!coeff.empty())
      coefficients = coeff;
    else {
      coefficients = orthogonalBasisProjection(initWaveFunc);
    }
  }

  std::vector<double> orthogonalBasisProjection(double waveFunc(double)) {
    std::vector<double> tempCoeff;
    tempCoeff.resize(hilbertSpace->dim);

    gsl_function gslFunc;
    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(1000);
    for (int i = 0; i < tempCoeff.size(); i++) {
      auto ptr = [=](double x) {
        return (std::conj(hilbertSpace->eigenbasis(i, x)) * waveFunc(x)).real();
      };
      gsl_function_pp<decltype(ptr)> Fp(ptr);
      gsl_function *gslFunc = static_cast<gsl_function *>(&Fp);
      gsl_integration_cquad(gslFunc, -INFINITY, INFINITY, 1.49e-08, 1.49e-08, w, &integrationResult,
                            NULL, NULL);
      tempCoeff[i] = integrationResult;
    }
    return tempCoeff;
  }

  double integrationResult;

  HilbertSpace *hilbertSpace;
  std::vector<double> coefficients;
  double (*evaluate)(double);
};