#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <time.h>

#include <fstream>
#include <functional>
#include <random>
#include <vector>

#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

// QHOS
#include "Wavefunction.hpp"
#include "shadercode.hpp"

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

  virtual void onResize(int w, int h) override;

  /**
   * @brief Do something when a key is pressed
   */
  virtual bool onKeyDown(al::Keyboard const& k) override;

  void onMessage(osc::Message& m) override;

  // Variables
  Mesh waveFunctionPlot;  // complex valued wave function
  Mesh probabilityPlot;   // probability distribution
  Mesh axes;
  Mesh grid;

  Mesh samples;
  std::array<Vec3d, 20> samplePoints;
  int sampleCounter = 0;

  ShaderProgram lineShader;
  ShaderProgram pointShader;

  Texture lineTexture;
  Texture pointTexture;

  FBO renderTarget;
  Texture rendered;

  void updateFBO(int w, int h) {
    rendered.create2D(w, h);
    renderTarget.bind();
    renderTarget.attachTexture2D(rendered);
    renderTarget.unbind();
  }

  // some variables to keep track of states
  double simTime = 0;
  float tableReader = 0;
  bool drawGUI = 1;

  // simulation parameters
  ParameterInt dims{"Dimensions", "", 2, "", 1, 20};
  ParameterBool coeffList{"Manual Coefficient Entry", "", 0};
  ParameterMenu presetFuncs{"Function Presets"};
  Parameter simSpeed{"simulation speed", 1.0, -5.0, 5.0};

  // audio parameters
  ParameterBool audioOn{"Audio On", "", 0};
  Parameter freq{"Frequency", 200, 30, 1000};
  Parameter volume{"Volume", 0.5, 0, 1};

  ParameterMenu sourceOne{"Wavetable 1 source"};
  ParameterMenu sourceTwo{"Wavetable 2 source"};
  // drawing parameters
  ParameterBool drawGrid{"Grid", "", 1};
  ParameterBool drawAxes{"Axes", "", 1};
  ParameterBool drawWavefunction{"Wave function", "", 1};
  ParameterBool drawProbability{"Probability", "", 1};
  ParameterBool drawMeasurements{"Measurements", "", 1};

  // wave function for generating coefficients
  // std::function has some undesirable overhead, I'd rather use a function pointer or lambda, but
  // if I ever want to capture some member variable to factor into the function, (for example to
  // do something like "x = 1 for highest eigenstate only" I need std::function.
  std::function<double(double x)> initWaveFunction = [&](double x) {
    if (abs(x) < (2 * M_PI))
      return sin(x);
    else
      return 0.0;
  };
  // some arrays used for plotting
  static const int resolution = 256;
  std::vector<double> posValues = linspace(-5, 5, resolution);
  std::array<std::complex<double>, resolution> psiValues;
  // These arrays store values for copying to the audio wavetable
  std::array<double, resolution> reValues;
  std::array<double, resolution> imValues;
  std::array<double, resolution> probValues;
  std::default_random_engine generator;
  // audio wavetable
  std::array<std::array<double, resolution>, 2> wavetable;
  std::array<double, resolution>* wavetableOneSource = &reValues;
  std::array<double, resolution>* wavetableTwoSource = &reValues;

  // initialize the wave function
  WaveFunction psi{new HilbertSpace(dims), initWaveFunction};

  // copy the coefficients into our variable for manual coefficient control
  std::vector<double> manualCoefficients = psi.coefficients;

  // mutex lock to avoid audio thread grabbing a half-written wavetable
  std::mutex wavetableLock;

  // custom wrapper for double precision imgui sliders
  bool SliderDouble(const char* label, double* v, double v_min, double v_max,
                    const char* format = NULL, float power = 1.0f) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max, format, power);
  }

  int oscPort = 16447;             // osc port
  char oscAddr[10] = "127.0.0.1";  // ip address
  osc::Send oscClient;             // create an osc client
  osc::Recv oscServer;
  ParameterBool oscOn{"OSC On", "OSC", 0};
  ParameterBool oscWaveform{"Send Waveform", "OSC", 0};
  ParameterBool oscMeasurement{"Send Measurements", "OSC", 0};

  ParameterBool measurementTrigger{"Measure", "", 0};

  void resetOSC() {
    oscClient.open(oscPort, oscAddr);
    std::cout << "New OSC port Selected: \n" + oscPort << std::endl;
    std::cout << oscAddr << std::endl;
  }
};