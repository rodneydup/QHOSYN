#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

// for master branch
// #include "al/core.hpp"

// for devel branch
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
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
  virtual void onSound(al::AudioIOData &io) override;

  virtual void onAnimate(double dt) override;

  /**
   * @brief Draw rate processing of synth interface.
   */
  virtual void onDraw(al::Graphics &g) override;

  // Variables
  Mesh plot;
  Mesh axes;

  double time = 0;

  std::vector<double> posValues = linspace(-5, 5, 300);
  std::vector<double> reValues;
  std::vector<double> imValues;
  std::vector<double> probValues;
  std::vector<std::complex<double>> amplitudeValues;
  HilbertSpace hs{2};

  WaveFunction psi{&hs, [](double x) {
                     if (abs(x) < 5)
                       return 1.0;
                     else
                       return 0.0;
                   }};
};