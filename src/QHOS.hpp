#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <time.h>

#include <functional>

// for master branch
// #include "al/core.hpp"

// for devel branch
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

// QHOS
#include "Wavefunction.hpp"

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
  virtual void onSound(al::AudioIOData& io) override;

  virtual void onAnimate(double dt) override;

  /**
   * @brief Draw rate processing of synth interface.
   */
  virtual void onDraw(al::Graphics& g) override;

  virtual bool onKeyDown(al::Keyboard const& k) override;

  // Variables
  Mesh plot;
  Mesh axes;

  double simTime = 0;
  int tableReader = 0;
  bool drawGUI = 1;

  // int dims = 10;
  ParameterInt dims{"Dimensions", "", 2, "", 1, 10};
  ParameterBool coeffList{"Manual Coefficient Entry", "", 0};

  // double (*initWaveFunction)(double x) = [](double x) { return 1.0; };
  std::function<double(double x)> initWaveFunction = [&](double x) {
    if (abs(x) < (2 * M_PI))
      return sin(x);
    else
      return 0.0;
  };

  std::vector<double> posValues = linspace(-5, 5, 300);
  std::vector<double> reValues;
  std::vector<double> imValues;
  std::vector<double> probValues;
  std::vector<std::complex<double>> amplitudeValues;
  HilbertSpace hs{dims};

  WaveFunction psi{&hs, initWaveFunction};
  std::vector<double> manualCoefficients = psi.coefficients;

  // WaveFunction psi{&hs, [](double x) {
  //                    if (abs(x) < 5)
  //                      return 1.0;
  //                    else
  //                      return 0.0;
  //                  }};

  // WaveFunction psi{&hs, [](double x) {
  //                    if (abs(x) < (2 * M_PI))
  //                      return sin(x);
  //                    else
  //                      return 0.0;
  //                  }};
  // WaveFunction psi{&hs, [](double x) {
  //                    srand(std::time(0));
  //                    return double(rand()) / RAND_MAX;
  //                  }};
  bool SliderDouble(const char* label, double* v, double v_min, double v_max,
                    const char* format = NULL, float power = 1.0f) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max, format, power);
  }
};