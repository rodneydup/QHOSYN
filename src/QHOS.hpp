#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <complex.h>
#include <fftw3.h>
#include <time.h>

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "Gamma/DFT.h"
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/io/al_MIDI.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "constant-q/cq.h"

// QHOS
#include "Wavefunction.hpp"
#include "shadercode.hpp"

using namespace al;

class QHOS : public App, public MIDIMessageHandler {
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

  // Meshes for visual display
  Mesh waveFunctionPlot;  // complex valued wave function
  Mesh probabilityPlot;   // probability distribution
  Mesh axes;
  Mesh grid;
  Mesh samples;

  std::array<Vec3d, 20> measurementPoints;
  std::array<int, 20> sampleDisplayTimer;
  int measurementPointsCounter = 0;
  float automeasureTimer = 0;

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

  // GUI stuff

  ImGuiWindowFlags flags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoMove |
                           ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings |
                           ImGuiWindowFlags_AlwaysAutoResize | ImGuiCond_Once;

  // simulation parameters
  ParameterInt dims{"Dimensions", "", 2, "", 1, 15};
  ParameterBool coeffList{"Manual Coefficient Entry", "", 0};
  ParameterMenu presetFuncs{"Function Presets"};
  ParameterBool project{"Project in Orthogonal Basis", "", 1};
  Parameter simSpeed{"simulation speed", 1.0, -5.0, 5.0};

  // audio parameters
  ParameterBool audioOn{"Audio On", "", 0};
  ParameterBool panner{"Panner", "", 0};
  bool pannerTrigger[2] = {0, 0};
  Parameter wavetableFreq{"Wavetable Frequency", 200, 10, 1000};
  ParameterInt ifftCenterFreq{"Start Bin", "", 1, "Bin ", 1, fftSize / 2};
  ParameterInt ifftBandwidth{"Bandwidth", "", 1, 1, 32};
  Parameter volume{"Volume", 0.0, 0, 1};
  ParameterMenu sourceOneMenu{"Channel 1 source"};
  ParameterMenu sourceTwoMenu{"Channel 2 source"};
  SmoothValue<float> pan[2]{{20, "log"}, {20, "log"}};

  // drawing parameters
  ParameterBool drawGrid{"Grid", "", 1};
  ParameterBool drawAxes{"Axes", "", 1};
  ParameterBool drawWavefunction{"Wave function", "", 1};
  ParameterBool drawProbability{"Probability", "", 1};
  ParameterBool drawMomentum{"Momentum Distribution", "", 1};

  ParameterBool drawMeasurements{"Measurements", "", 1};

  // OSC parameters
  ParameterBool oscOn{"OSC On", "OSC", 0};
  ParameterBool oscWaveform{"Send Waveform", "OSC", 0};
  ParameterBool oscMeasurement{"Send Measurements", "OSC", 0};

  // Measurement parameters
  ParameterBool measurementTrigger{"Measure", "", 0};
  ParameterBool automeasure{"Auto-measure", "", 0};
  Parameter automeasureInterval{"Auto-measure Interval (s)", 0.1, 0.016, 2};

  // custom wrapper for double precision imgui sliders
  bool SliderDouble(const char* label, double* v, double v_min, double v_max,
                    const char* format = NULL, float power = 1.0f) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max, format, power);
  }

  // wave function for generating coefficients
  // std::function has some undesirable overhead, I'd rather use a function pointer or lambda, but
  // if I ever want to capture some member variable to factor into the function, (for example to
  // do something like "x = 1 for highest eigenstate only" I need std::function.
  std::function<double(double x)> initWaveFunction = [&](double x) { return sin(x); };
  // initialize the wave function
  WaveFunction psi{new HilbertSpace(dims), initWaveFunction};

  // copy the coefficients into our variable for manual coefficient control
  std::vector<double> manualCoefficients = psi.coefficients;

  // size of all our arrays
  static const int resolution = 256;

  // range of position values to report for visualization and sonification {-5,5}
  std::vector<double> posValues = linspace<double>(-5, 5, resolution);

  // array to store results from looking up wavefunction values
  std::array<std::complex<double>, resolution> psiValues;

  // These arrays store values for copying to the audio wavetable
  std::array<float, resolution> reValues;
  std::array<float, resolution> imValues;
  std::array<float, resolution> probValues;
  std::array<float, resolution> emptyTable;

  // // ifft
  int sourceSelect[2] = {0, 0};
  static const int fftSize = resolution * 32;
  float binSize = 44100 / fftSize;

  // // fftw inverse fft
  fftw_complex* c2rIn;
  double* c2rOut;
  fftw_plan c2r;
  static const int M = resolution * 8;
  std::array<float, M> filterKernelWindow;

  double* r2cIn;
  fftw_complex* r2cOut;
  fftw_plan r2c;

  // // fftw forward fft
  fftw_complex* c2cForwardIn;
  fftw_complex* c2cForwardOut;
  fftw_plan c2cForward;
  Mesh momentumPlot;

  std::ofstream testOut{"test.csv"};
  // std::ofstream outputpost{"post-ifft.csv"};

  static const int bufferLen = fftSize / 4;
  std::array<std::array<std::array<float, bufferLen>, 4>, 2> ifftBuffers;
  std::array<float, bufferLen> overlapAddWindow;
  int fftBufferLenDiff = fftSize - bufferLen;

  int currentIfftBuffer[2] = {0, 0};
  int currentIfftSample[2] = {0, 0};

  // audio wavetable
  // output samples
  float sample[2] = {0, 0};
  std::array<std::array<float, resolution>, 2> wavetable;
  std::array<float, resolution>* source[2] = {&reValues, &imValues};

  // // mutex lock to avoid audio thread grabbing a half-written wavetable
  // std::mutex wavetableLock;

  // RNG to use in weighted distribution (to take a measurement of wavefunction)
  // I checked that this worked by taking many results and plotting it (July 20, 2021)
  std::random_device randomDevice;
  std::mt19937 gen{randomDevice()};
  std::normal_distribution<> norm{0, 1};
  std::uniform_real_distribution<> uniform{0, 1000};

  // OSC stuff

  int oscClientPort = 16447;             // osc port
  char oscClientAddr[10] = "127.0.0.1";  // ip address
  int oscServerPort = 16448;             // osc port
  char oscServerAddr[10] = "127.0.0.1";  // ip address
  osc::Send oscClient;                   // create an osc client (broadcaster)
  osc::Recv oscServer;                   // create an osc Server (listener)

  void resetOSC() {
    oscClient.open(oscClientPort, oscClientAddr);
    std::cout << "New OSC port Selected: \n" + oscClientPort << std::endl;
    std::cout << oscClientAddr << std::endl;
  }

  // MIDI stuff

  RtMidiIn midiIn;

  void onMIDIMessage(const MIDIMessage& m) {
    printf("%s: ", MIDIByte::messageTypeString(m.status()));

    switch (m.type()) {
      case MIDIByte::NOTE_ON:
        if (!m.velocity() == 0) wavetableFreq = bpScale[m.noteNumber() - 36] * 10;
        printf("Note %u, Vel %f", m.noteNumber(), m.velocity());
        break;

      case MIDIByte::NOTE_OFF:
        printf("Note %u, Vel %f", m.noteNumber(), m.velocity());
        break;

      case MIDIByte::PITCH_BEND:
        printf("Value %f", m.pitchBend());
        break;

      // Control messages need to be parsed again...
      case MIDIByte::CONTROL_CHANGE:
        printf("%s ", MIDIByte::controlNumberString(m.controlNumber()));
        switch (m.controlNumber()) {
          case MIDIByte::MODULATION:
            printf("%f", m.controlValue());
            break;

          case MIDIByte::EXPRESSION:
            printf("%f", m.controlValue());
            break;
        }
        break;
      default:;
    }

    // If it's a channel message, print out channel number
    if (m.isChannelMessage()) {
      printf(" (MIDI chan %u)", m.channel() + 1);
    }

    printf("\n");

    // Print the raw byte values and time stamp
    printf("\tBytes = ");
    for (unsigned i = 0; i < 3; ++i) {
      printf("%3u ", (int)m.bytes[i]);
    }
    printf(", time = %g\n", m.timeStamp());
  }

  virtual void onExit() {
    // destroy inverse fftw things
    fftw_destroy_plan(c2r);
    fftw_free(c2rIn);
    fftw_free(c2rOut);

    // destroy forward fftw things
    fftw_destroy_plan(r2c);
    fftw_free(r2cIn);
    fftw_free(r2cOut);

    // destroy forward fftw things
    fftw_destroy_plan(c2cForward);
    fftw_free(c2cForwardIn);
    fftw_free(c2cForwardOut);
  };
};
