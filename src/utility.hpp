#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <gsl/gsl_integration.h>
#include <unistd.h>

#include <vector>

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
  SampleFunction(T f, float min, float max, int numSamples) {
    minX = min;
    maxX = max;
    rangeX = max - min;
    numSamplesX = numSamples;

    std::vector<double> domain = linspace(minX, maxX, numSamplesX);
    image.resize(domain.size());
    for (int x = 0; x < domain.size(); x++) image[x] = f(domain[x]);
  }

  double lookup(float x) {
    if (isnan(x)) return 0;
    x = std::min(std::max(x, minX), maxX - 1);
    x -= minX;
    x = std::round(x * (numSamplesX / rangeX));
    return image[int(x)];
  }

 private:
  int numSamplesX;
  float minX, maxX, rangeX;
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