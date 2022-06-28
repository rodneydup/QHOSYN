/*
Quantum Harmonic Oscillator Synthesizer (QHOSYN)
By Rodney DuPlessis (2021)
*/

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

// C libraries
#include <complex.h>
#include <time.h>

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <array>

// Native File Dialog
#include "../external/nativefiledialog/src/include/nfd.h"

// FFTW
#include <fftw3.h>

// Allolib
#include "../external/allolib/external/Gamma/Gamma/Filter.h"
#include "Gamma/tbl.h"
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/io/al_MIDI.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/io/al_File.hpp"
#include "al_ext/soundfile/al_OutputRecorder.hpp"

// QHOSYN
#include "shadercode.hpp"
#include "wavefunction.hpp"

using namespace al;

class QHOSYN : public App, public MIDIMessageHandler {
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

  static const int MAX_AUDIO_OUTS = 2;

  std::string execDir;

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
  ParameterInt dims{"Eigenstates", "", 2, "", 1, 15};
  ParameterBool coeffList{"Manual Coefficient Entry", "", 0};
  ParameterMenu presetFunctions{"Function Presets"};
  std::string customFunction = "";
  ParameterBool project{"Project in Orthogonal Basis", "", 1};
  Parameter simSpeed{"Simulation Speed", 1.0, -5.0, 5.0};

  // audio parameters
  ParameterBool audioOn{"Audio On", "", 0};
  ParameterBool panner{"Panner", "", 0};
  bool pannerTrigger[2] = {0, 0};
  Parameter wavetableFreq{"Frequency", 200, 10, 1000};
  ParameterInt ifftBin{"Start Bin", "", 1, "Bin ", 1, fftSize / 4};
  ParameterInt ifftBandwidth{"Bandwidth", "", 1, "", 1, 16};
  Parameter volume{"Volume", 0.0, 0, 1};
  ParameterMenu sourceOneMenu{"Channel 1 source"};
  ParameterMenu sourceTwoMenu{"Channel 2 source"};
  SmoothValue<float> pan[2];

  // drawing parameters
  ParameterBool drawGrid{"Grid", "", 1};
  ParameterBool drawAxes{"Axes", "", 1};
  ParameterBool drawWavefunction{"Wave function", "", 1};
  ParameterBool drawProbability{"Probability", "", 1};
  ParameterBool drawIfft{"Inverse Fourier Waveform", "", 1};
  ParameterBool drawNoise{"Noise Waveform", "", 1};
  ParameterBool drawMeasurements{"Measurements", "", 1};

  // OSC parameters
  ParameterBool oscSenderOn{"OSC Sender On", "OSC", 0};
  ParameterBool oscReceiverOn{"OSC Receiver On", "OSC", 0};
  ParameterBool oscWaveform{"Send Waveform", "OSC", 0};
  ParameterBool oscMeasurement{"Send Measurements", "OSC", 0};

  // Measurement parameters
  ParameterBool measurementTrigger{"Measure", "", 0};
  ParameterBool automeasure{"Auto-measure", "", 0};
  Parameter automeasureInterval{"Interval (s)", 0.1, 0.016, 2};

  gam::Biquad<> filter[2];
  gam::Biquad<> antiAlias[2];
  ParameterBool filterOn{"Filter", "", 1};

  // custom wrapper for double precision imgui sliders
  bool SliderDouble(const char* label, double* v, double v_min, double v_max,
                    const char* format = NULL, float power = 1.0f) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max, format, power);
  }

  // wave function for generating coefficients
  // std::function has some undesirable overhead, I'd rather use a function pointer or
  // lambda, but if I ever want to capture some member variable to factor into the
  // function, (for example to do something like "x = 1 for highest eigenstate only" I
  // need std::function.
  std::function<double(double x)> initWaveFunction = [&](double x) { return sin(x); };

  // Create Hilbert basis with 15 dimensions
  HilbertSpace basis{15};

  // initialize the wave function
  WaveFunction psi{&basis, 2, initWaveFunction};

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
  float binSize = 48000 / fftSize;

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
  Mesh ifftPlot;
  Mesh noisePlot;

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
  std::unique_ptr<osc::Send> oscSender;     // create an osc Sender
  int oscSenderPort = 16447;                // osc port
  std::string oscSenderAddr = "127.0.0.1";  // ip address

  std::unique_ptr<al::osc::Recv> oscReceiver;             // create an osc Receiver
  int oscReceiverPort = 16448;                            // osc port
  std::string oscReceiverAddr = "127.0.0.1";              // ip address
  int previousoscReceiverPort = oscReceiverPort;          // osc port
  std::string previousoscReceiverAddr = oscReceiverAddr;  // ip address
  float oscReceiverTimeout = 0.02;

  std::string measurementArg = "/measurement";
  ImGuiInputTextCallback inputTextCallback;
  void* CallbackUserData;

  parseStringToMathFunction parser{"x", "x"};

  std::string opener = "open ";
  std::string functionParserURL =
    "http://warp.povusers.org/FunctionParser/fparser.html#functionsyntax";

  bool isOscWarningWindow = false;

  void resetOSC() {
    if (oscReceiver != nullptr) oscReceiver->stop();
    oscReceiver.reset();
    oscReceiver =
      std::make_unique<al::osc::Recv>(oscReceiverPort, oscReceiverAddr.c_str(), oscReceiverTimeout);
    if (oscReceiver->isOpen()) {
      oscReceiver->handler(oscDomain()->handler());
      oscReceiver->start();
      std::cout << "OSC Receiver Settings:" << std::endl;
      std::cout << "IP Address: " << oscReceiverAddr << std::endl;
      std::cout << "Port: " << oscReceiverPort << std::endl;
      std::cout << "Timeout: " << oscReceiverTimeout << std::endl;
      previousoscReceiverAddr = oscReceiverAddr;
      previousoscReceiverPort = oscReceiverPort;
    } else {
      std::cerr << "Could not bind to UDP socket. Is there a server already bound to that port?"
                << std::endl;
      oscReceiverAddr = previousoscReceiverAddr;
      oscReceiverPort = previousoscReceiverPort;
      oscReceiver.reset();
      oscReceiver = std::make_unique<al::osc::Recv>(oscReceiverPort, oscReceiverAddr.c_str(), 0.02);
      isOscWarningWindow = true;
    }
    oscSender->open(oscSenderPort, oscSenderAddr.c_str());
    std::cout << "OSC Sender Settings:" << std::endl;
    std::cout << "IP Address: " << oscSenderAddr << std::endl;
    std::cout << "Port: " << oscSenderPort << std::endl;
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

  ImFont* bodyFont;
  ImFont* titleFont;

  // for ImGUI variable length string input fields
  struct InputTextCallback_UserData {
    std::string* Str;
    ImGuiInputTextCallback ChainCallback;
    void* ChainCallbackUserData;
  };

  static int InputTextCallback(ImGuiInputTextCallbackData* data) {
    InputTextCallback_UserData* user_data = (InputTextCallback_UserData*)data->UserData;
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
      // Resize string callback
      // If for some reason we refuse the new length (BufTextLen) and/or capacity
      // (BufSize) we need to set them back to what we want.
      std::string* str = user_data->Str;
      IM_ASSERT(data->Buf == str->c_str());
      str->resize(data->BufTextLen);
      data->Buf = (char*)str->c_str();
    } else if (user_data->ChainCallback) {
      // Forward to user callback, if any
      data->UserData = user_data->ChainCallbackUserData;
      return user_data->ChainCallback(data);
    }
    return 0;
  };

  bool InputText(const char* label, std::string* str, ImGuiInputTextFlags flags,
                 ImGuiInputTextCallback callback, void* user_data) {
    IM_ASSERT((flags & ImGuiInputTextFlags_CallbackResize) == 0);
    flags |= ImGuiInputTextFlags_CallbackResize;

    InputTextCallback_UserData cb_user_data;
    cb_user_data.Str = str;
    cb_user_data.ChainCallback = callback;
    cb_user_data.ChainCallbackUserData = user_data;
    return ImGui::InputText(label, (char*)str->c_str(), str->capacity() + 1, flags,
                            InputTextCallback, &cb_user_data);
  }

  std::string currentAudioDeviceOut;
  std::array<unsigned int, MAX_AUDIO_OUTS> AudioChanIndexOut;
  int getLeadChannelOut() const { return AudioChanIndexOut[0]; }
  bool isPaused = false;
  double globalSamplingRate = 48000;
  const int BLOCK_SIZE = 1024;
  std::string soundOutput;
  al::OutputRecorder mRecorder;
  nfdresult_t result;
  nfdchar_t* outPath = NULL;

  int getSampleRateIndex() {
    unsigned s_r = (unsigned)globalSamplingRate;
    switch (s_r) {
      case 44100:
        return 0;
      case 48000:
        return 1;
      case 88200:
        return 2;
      case 96000:
        return 3;
      default:
        return 0;
    }
  }
  
  void setSoundOutputPath(std::string sound_output_path) {
    soundOutput = al::File::conformPathToOS(sound_output_path);
  }

  void setAudioSettings(float sample_rate) {
    globalSamplingRate = sample_rate;
    configureAudio(sample_rate, BLOCK_SIZE, MAX_AUDIO_OUTS, 0);
  }

  void setOutChannels(int lead_channel, int max_possible_channels) {
    AudioChanIndexOut[0] = lead_channel;
    if (max_possible_channels == 1) {
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel;
      }
    } else {
      // assert(lead_channel + (consts::MAX_AUDIO_OUTS) < max_possible_channels);
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel + i;
      }
    }
  }

  void drawAudioIO(AudioIO* io) {
    struct AudioIOState {
      int currentSr = 1;
      int currentBufSize = 3;
      int currentDeviceOut = 0;
      int currentDeviceIn = 0;
      int currentOut = 1;
      int currentIn = 1;
      int currentMaxOut;
      int currentMaxIn;
      std::vector<std::string> devices;
    };

    auto updateOutDevices = [&](AudioIOState& state) {
      state.devices.clear();
      int numDevices = AudioDevice::numDevices();
      int dev_out_index = 0;
      for (int i = 0; i < numDevices; i++) {
        if (!AudioDevice(i).hasOutput()) continue;

        state.devices.push_back(AudioDevice(i).name());
        if (currentAudioDeviceOut == AudioDevice(i).name()) {
          state.currentDeviceOut = dev_out_index;
          state.currentOut = getLeadChannelOut() + 1;
          state.currentMaxOut = AudioDevice(i).channelsOutMax();
        }
        dev_out_index++;
      }
    };

    static std::map<AudioIO*, AudioIOState> stateMap;
    if (stateMap.find(io) == stateMap.end()) {
      stateMap[io] = AudioIOState();
      updateOutDevices(stateMap[io]);
    }
    AudioIOState& state = stateMap[io];
    ImGui::PushID(std::to_string((unsigned long)io).c_str());

    if (io->isOpen()) {
      std::string text;
      text += "Output Device: " + state.devices.at(state.currentDeviceOut);
      // text += "\nInput Device: " + state.devices.at(state.currentDeviceIn);
      text += "\nSampling Rate: " + std::to_string(int(io->fps()));
      text += "\nBuffer Size: " + std::to_string(io->framesPerBuffer());
      text += "\nOutput Channels: " + std::to_string(state.currentOut) + ", " +
              std::to_string(state.currentOut + 1);
      // text += "\nInput Channels: " + std::to_string(state.currentIn) + ", " +
      //         std::to_string(state.currentIn + 1);
      ImGui::Text("%s", text.c_str());
      if (ImGui::Button("Stop")) {
        isPaused = true;
        io->stop();
        io->close();
        state.currentSr = getSampleRateIndex();
      }
    } else {
      if (ImGui::Button("Update Devices")) {
        updateOutDevices(state);
      }

      ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
      if (ImGui::Combo("Output Device", &state.currentDeviceOut, ParameterGUI::vector_getter,
                       static_cast<void*>(&state.devices), state.devices.size())) {
        state.currentMaxOut =
          AudioDevice(state.devices.at(state.currentDeviceOut), AudioDevice::OUTPUT)
            .channelsOutMax();
      }
      std::string chan_label_out =
        "Select Outs: (Up to " + std::to_string(state.currentMaxOut) + " )";
      ImGui::Text(chan_label_out.c_str(), "%s");
      // ImGui::SameLine();
      // ImGui::Checkbox("Mono/Stereo", &isStereo);
      // ImGui::Indent(25 * fontScale);
      // ImGui::PushItemWidth(50 * fontScale);
      ImGui::DragInt("Chan 1 out", &state.currentOut, 1.0f, 0, state.currentMaxOut - 1, "%d",
                     1 << 4);

      if (state.currentOut > state.currentMaxOut - 1) state.currentOut = state.currentMaxOut - 1;
      if (state.currentOut < 1) state.currentOut = 1;

      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        ImGui::SameLine();
        int temp = state.currentOut + i;
        std::string channel = "Chan " + std::to_string(i + 1);
        ImGui::DragInt(channel.c_str(), &temp, 1.0f, 0, state.currentMaxOut, "%d", 1 << 4);
      }
      ImGui::PopStyleVar();

      // ImGui::Unindent(25 * fontScale);
      ImGui::PopItemWidth();

      std::vector<std::string> samplingRates{"44100", "48000", "88200", "96000"};
      ImGui::Combo("Sampling Rate", &state.currentSr, ParameterGUI::vector_getter,
                   static_cast<void*>(&samplingRates), samplingRates.size());
      if (ImGui::Button("Start")) {
        globalSamplingRate = std::stof(samplingRates[state.currentSr]);
        io->framesPerSecond(globalSamplingRate);
        io->framesPerBuffer(BLOCK_SIZE);
        io->deviceOut(AudioDevice(state.devices.at(state.currentDeviceOut), AudioDevice::OUTPUT));
        currentAudioDeviceOut = state.devices.at(state.currentDeviceOut);
        setOutChannels(state.currentOut - 1, state.currentMaxOut);
        io->open();
        io->start();
        isPaused = false;
      }
      ImGui::SameLine();
    }
    ImGui::PopID();
  }

  void drawRecorderWidget(al::OutputRecorder* recorder, double frameRate, uint32_t numChannels,
                          std::string directory, uint32_t bufferSize) {
    struct SoundfileRecorderState {
      bool recordButton;
      bool overwriteButton;
    };
    static std::map<SoundFileBufferedRecord*, SoundfileRecorderState> stateMap;
    if (stateMap.find(recorder) == stateMap.end()) {
      stateMap[recorder] = SoundfileRecorderState{0, false};
    }
    SoundfileRecorderState& state = stateMap[recorder];
    ImGui::PushID(std::to_string((unsigned long)recorder).c_str());
    ImGui::Text("Output File Name:");
    static char buf1[64] = "test.wav";
    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 10.0f);
    ImGui::InputText("##Record Name", buf1, 63);
    ImGui::PopItemWidth();

    if (state.recordButton) {
      ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.9, 0.3, 0.3, 1.0));
      ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.8, 0.5, 0.5, 1.0));
      ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0, 1.0, 1.0, 1.0));
    }
    std::string buttonText = state.recordButton ? "Stop" : "Record";
    bool recordButtonClicked = ImGui::Button(buttonText.c_str());
    if (state.recordButton) {
      ImGui::PopStyleColor();
      ImGui::PopStyleColor();
      ImGui::PopStyleColor();
    }
    if (recordButtonClicked) {
      state.recordButton = !state.recordButton;
      if (state.recordButton) {
        uint32_t ringBufferSize;
        if (bufferSize == 0) {
          ringBufferSize = 8192;
        } else {
          ringBufferSize = bufferSize * numChannels * 4;
        }
        std::string filename = buf1;
        if (!state.overwriteButton) {
          int counter = 1;
          while (File::exists(directory + filename) && counter < 9999) {
            filename = buf1;
            int lastDot = filename.find_last_of(".");
            filename = filename.substr(0, lastDot) + "_" + std::to_string(counter++) +
                       filename.substr(lastDot);
          }
        }
        if (!recorder->start(directory + filename, frameRate, numChannels, ringBufferSize,
                             gam::SoundFile::WAV, gam::SoundFile::FLOAT)) {
          std::cerr << "Error opening file for record" << std::endl;
        }
      } else {
        recorder->close();
      }
    }
    ImGui::SameLine();
    ImGui::Checkbox("Overwrite", &state.overwriteButton);
    ImGui::Text("Writing to:");
    ImGui::TextWrapped("%s", directory.c_str());

    ImGui::PopID();
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
