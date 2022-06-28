#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <gsl/gsl_integration.h>

#include <vector>
#include <iostream>

#ifdef _WIN32
#include <stdlib.h>
#define strcasecmp _stricmp
#define PATH_MAX 1024
#elif __linux__
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#else
#include <assert.h>
#include <limits.h>
#include <mach-o/dyld.h>
#include <stdlib.h>
#endif

#include "../external/fparser4.5.2/fparser.hh"

// Calculate a factorial
int factorial(int x) {
  int result = 1;
  for (int i = 1; i <= x; i++) result = result * i;
  return result;
}

// create a discretized linear space in range {start_in, end_in} with num_in points
template <typename T>
std::vector<T> linspace(T start_in, T end_in, int num_in) {
  std::vector<T> linspaced;

  T start = static_cast<T>(start_in);
  T end = static_cast<T>(end_in);
  T num = static_cast<T>(num_in);

  if (num == 0) {
    return linspaced;
  }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }

  T delta = (end - start) / (num - 1);

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

    std::vector<double> domain = linspace<double>(minX, maxX, numSamplesX);
    functionSamples.resize(domain.size());
    for (int x = 0; x < domain.size(); x++) functionSamples[x] = f(domain[x]);
  }

  double lookup(float x) {
    if (isnan(x)) return 0;
    x = std::min(std::max(x, minX), maxX - 1);
    x -= minX;
    x = std::round(x * (numSamplesX / rangeX));
    return functionSamples[int(x)];
  }

  int size() { return functionSamples.size(); }

 private:
  int numSamplesX;
  float minX, maxX, rangeX;
  std::vector<double> functionSamples;
};

// some black magic to cast a function pointer in a way that gsl likes
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

// Bohlen-pierce ratios for tuning midi input
std::vector<float> bpScale{
  1.000000,  1.080000,  1.190476,  1.285714,  1.400000,  1.530612,  1.666667,  1.800000,  1.960000,
  2.142857,  2.333333,  2.520000,  2.777778,  3.000000,  3.240000,  3.571429,  3.857143,  4.200000,
  4.591837,  5.000000,  5.400000,  5.880000,  6.428571,  7.000000,  7.560000,  8.333333,  9.000000,
  9.720000,  10.714286, 11.571429, 12.600000, 13.775510, 15.000000, 16.200000, 17.640000, 19.285714,
  21.000000, 22.680000, 25.000000, 27.000000, 29.160000, 32.130000, 34.830000, 37.800000, 41.310000,
  45.090000, 48.600000, 52.920000, 57.780000, 62.910000, 68.040000, 75.060000, 81.000000};

// circular buffer with some useful methods
template <typename T>
class RingBuffer {
 public:
  RingBuffer(unsigned maxSize) : mBufferSize(maxSize) {
    mBuffer.resize(mBufferSize);
    mWriteHead = -1;
    mReadHead = -1;
    mPrevSample = 0;
  }

  // return the size of the RingBuffer
  unsigned size() const { return mBufferSize; }

  // resize the RingBuffer
  void resize(unsigned maxSize) {
    mBufferSize = maxSize;
    mBuffer.resize(mBufferSize);
  }

  // increment the write head and write value at write head
  void push_back(T value) {
    mMutex.lock();
    mWriteHead = (mWriteHead + 1) % mBufferSize;
    mBuffer[mWriteHead] = value;
    mMutex.unlock();
  }

  // increment the write head and write value at write head
  void overwrite(int index, T value) {
    if (index >= mBufferSize) {
      std::cerr << "RingBuffer index out of range." << std::endl;
      index = index % mBufferSize;
    }
    mMutex.lock();
    mBuffer[index] = value;
    mMutex.unlock();
  }

  // return the index of the write head
  unsigned getWritehead() const { return mWriteHead; }

  // return the index of the read head
  unsigned getReadhead() const { return mReadHead; }

  // Read value at index
  T at(unsigned index) {
    if (index >= mBufferSize) {
      std::cerr << "RingBuffer index out of range." << std::endl;
      index = index % mBufferSize;
    }
    if (mMutex.try_lock()) {
      mPrevSample = mBuffer[index];
      mMutex.unlock();
    }
    return mPrevSample;
  }

  // read value at index and set value to 0.
  T consume(unsigned index) {
    if (index >= mBufferSize) {
      std::cerr << "RingBuffer index out of range." << std::endl;
      index = index % mBufferSize;
    }
    if (mMutex.try_lock()) {
      mPrevSample = mBuffer[index];
      mBuffer[index] = 0;
      mMutex.unlock();
    }
    return mPrevSample;
  }

  // increment read head, read value, and set value to 0.
  T consume() {
    mReadHead = (mReadHead + 1) % mBufferSize;
    if (mMutex.try_lock()) {
      mPrevSample = mBuffer[mReadHead];
      mBuffer[mReadHead] = 0;
      mMutex.unlock();
    }
    return mPrevSample;
  }

  // overridden [] operator to call at()
  T operator[](unsigned index) { return this->at(index); }

  // overridden () operator to call at() and increment read head
  T operator()() {
    mReadHead = (mReadHead + 1) % mBufferSize;
    return this->at(mReadHead);
  }

  // return a pointer to the buffer data
  const T *data() { return mBuffer.data(); }

  // get RMS value over a number of elements before mWriteHead
  T getRMS(unsigned lookBackLength) {
    if (lookBackLength > mBufferSize) {
      std::cerr << "Lookback length must be less than ringbuffer size. Setting lookback "
                   "to length "
                   "of ringbuffer"
                << std::endl;
      lookBackLength = mBufferSize;
    }
    int start = mWriteHead - lookBackLength;
    if (start < 0) start = mBufferSize + start;

    T val = 0.0;
    for (unsigned i = 0; i < lookBackLength; i++) {
      val += pow(mBuffer[(start + i) % mBufferSize], 2);
    }
    return sqrt(val / lookBackLength);
  }

  // send buffer values to standard output
  void print() const {
    for (auto i = mBuffer.begin(); i != mBuffer.end(); ++i) std::cout << *i << " ";
    std::cout << "\n";
  }

 private:
  std::vector<T> mBuffer;
  unsigned mBufferSize;
  int mWriteHead;
  int mReadHead;
  T mPrevSample;

  std::mutex mMutex;
};

// interpolate toward new values to smooth transitions
template <class T>  // The class is templated to allow a variety of data types
class SmoothValue {
 public:
  SmoothValue(float initialTime = 50.0f,
              std::string interpolationType =
                "log") {  // how much time, in ms, does it take to arrive (or approach target value)
    arrivalTime = initialTime;  // store the input value (default to 50ms)
    interpolation = interpolationType;
    calculateCoefficients();  // calculate a and b (it is a lowpass filter)
    z = 0.0f;
  }

  inline T process() {  // this function will be called per-sample on the data
    if (stepsTaken <= numSteps) {
      if (interpolation == "log") z = (targetValue * b) + (z * a);
      if (interpolation == "lin") z += linStep;
      stepsTaken += 1;
    }
    return z;  // return the new z value (the output sample)
  }

  void setTime(float newTime) {  // set time (in ms)
    arrivalTime = newTime;       // store the input value
    numSteps = arrivalTime * 0.001f * gam::sampleRate();
    if (interpolation == "log") calculateCoefficients();  // calculate a and b
  }

  void setTarget(T newTarget) {
    targetValue = newTarget;
    linStep = (targetValue - z) / numSteps;
    stepsTaken = 0;
  }

  // Jump to a value without ramping
  void setCurrentValue(T newValue) {
    targetValue = newValue;
    z = newValue;
    stepsTaken = numSteps;
  }

  // Change between log and lin ramp types
  void changeType(std::string newType) { interpolation = newType; }

  // Get current value without processing a step
  T getCurrentValue() { return z; }

  // Get the current target value
  T getTargetValue() { return targetValue; }

  // Get the current ramp time in milliseconds
  float getTime() { return arrivalTime; }

 private:
  void calculateCoefficients() {  // called only when 'setTime' is called (and in constructor)
    a = std::exp(-(M_PI * 2) / numSteps);  // rearranged lpf coeff calculations
    b = 1.0f - a;
  }

  T targetValue;      // what is the destination (of type T, determined by implementation)
  T currentValue;     // how close to the destination? (output value)
  float arrivalTime;  // how long to take
  float a, b;         // coefficients
  T linStep;          // size of step in linear interpolation
  T numSteps;         // number of steps in interpolation
  int stepsTaken;     // variable to keep track of steps to avoid overshooting
  T z;                // storage for previous value
  std::string interpolation;
};

// Parse a string to a function
double fRand(const double *range) {
  return (double(rand()) / (RAND_MAX / (range[1] - range[0]))) + range[0];
}
class parseStringToMathFunction {
 public:
  parseStringToMathFunction(std::string defaultVariables = "x", std::string defaultFunction = "x") {
    mVariables = defaultVariables;
    mFunction = defaultFunction;
    fp.AddFunction("rand", fRand, 2);
  }
  double Square(const double *p) { return p[0] * p[0]; }
  void setFunction(std::string newFunction) { mFunction = newFunction; }

  void setVariables(std::string newVariables) { mVariables = newVariables; }

  void parse() {
    int errorCode = fp.Parse(mFunction, mVariables);
    if (errorCode == -1)
      return;
    else {
      std::cout << "\033[1;31mFunction Parser Error \033[0m" << fp.ErrorMsg() << std::endl;
    }
  }

  // three evaluate functions, supporting up to three variables
  double evaluate(double input) {
    double inputArray[1] = {input};
    double val = fp.Eval(inputArray);
    int errorCode = fp.EvalError();
    switch (errorCode) {
      case 0:
        return val;
        break;
      case 1:
        std::cout << "\033[1;31mFunction Evaluation Error: \033[0m"
                  << "division by zero" << std::endl;
        return 0;
        break;
      case 2:
        std::cout << "\033[1;31mFunction Evaluation Error: \033[0m"
                  << "sqrt error (sqrt of a negative value)" << std::endl;
        return 0;
        break;
      case 3:
        std::cout << "\033[1;31mFunction Evaluation Error: \033[0m"
                  << "log error (logarithm of a negative value)" << std::endl;
        return 0;
        break;
      case 4:
        std::cout << "\033[1;31mFunction Evaluation Error: \033[0m"
                  << "trigonometric error (asin or acos of illegal value)" << std::endl;
        return 0;
        break;
      case 5:
        std::cout << "\033[1;31mFunction Evaluation Error: \033[0m"
                  << "maximum recursion level in eval() reached" << std::endl;
        return 0;
        break;
      default:
        break;
    }
  }

  double evaluate(double inputOne, double inputTwo) {
    double inputArray[2] = {inputOne, inputTwo};
    fp.Eval(inputArray);
  }

  double evaluate(double inputOne, double inputTwo, double inputThree) {
    double inputArray[3] = {inputOne, inputTwo, inputThree};
    fp.Eval(inputArray);
  }

 private:
  FunctionParser fp;
  std::string mVariables;
  std::string mFunction;
};

/*
 * Returns the full path to the currently running executable,
 * or an empty string in case of failure.
 */
std::string getExecutablePath() {
#if _WIN32
  char *exePath;
  if (_get_pgmptr(&exePath) != 0) exePath = "";

#elif __linux__
  char exePath[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", exePath, sizeof(exePath));
  if (len == -1 || len == sizeof(exePath)) len = 0;
  exePath[len] = '\0';
#else  // THIS MEANS YOU ARE USING A >
  char exePath[PATH_MAX];
  uint32_t len = sizeof(exePath);
  if (_NSGetExecutablePath(exePath, &len) != 0) {
    exePath[0] = '\0';  // buffer too small (!)
  } else {
    // resolve symlinks, ., .. if possible
    char *canonicalPath = realpath(exePath, NULL);
    if (canonicalPath != NULL) {
      strncpy(exePath, canonicalPath, len);
      free(canonicalPath);
    }
  }
#endif
  return std::string(exePath);
}

std::string getContentPath_OSX(std::string s) {
  char delim = '/';
  size_t counter = 0;
  size_t i = s.size() - 1;
  while (counter < 2) {
    if (s[i] == delim) counter++;
    i--;
  }
  return s.substr(0, i + 2);
}