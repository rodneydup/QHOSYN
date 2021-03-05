// Quantum Harmonic Oscillator Sonifier
// By Rodney DuPlessis
// To do:
// - make new wavefunction generation asynchronous so it doesn't freeze things up
// - find more interesting-sounding situations
// - add arbitrary function input option

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOS.hpp"

void QHOS::onInit() {
  imguiInit();
  dims.registerChangeCallback([&](int x) {
    manualCoefficients.resize(dims);
    psi.newHilbertSpace(new HilbertSpace(x), initWaveFunction);
    manualCoefficients = psi.coefficients;
  });
  coeffList.registerChangeCallback([&](bool x) {
    if (x) {
      std::cout << "true" << std::endl;
      psi.newHilbertSpace(new HilbertSpace(dims), NULL, manualCoefficients);
    } else {
      psi.newHilbertSpace(new HilbertSpace(dims), initWaveFunction);
    }
  });
  presetFuncs.setElements({"psi = sin(x) if |x| < 2*pi; 0 otherwise",
                           "psi = 1 if |x|<5; 0 otherwise", "psi = random(0 to 1)"});
  presetFuncs.registerChangeCallback([&](int x) {
    switch (x) {
      case 0:
        initWaveFunction = [](double x) {
          if (abs(x) < (2 * M_PI))
            return sin(x);
          else
            return 0.0;
        };
        break;
      case 1:
        initWaveFunction = [](double x) {
          if (abs(x) < 5)
            return 1.0;
          else
            return 0.0;
        };
        break;
      case 2:
        initWaveFunction = [](double x) {
          srand(std::time(0));
          return double(rand()) / RAND_MAX;
        };
        break;

      default:
        break;
    }
    psi.newHilbertSpace(new HilbertSpace(dims), initWaveFunction);
  });
  sourceOne.setElements({"Real Values", "Imaginary Values", "Probability Values"});
  sourceOne.registerChangeCallback([&](int x) {
    switch (x) {
      case 0:
        wavetableOneSource = &reValues;
        break;
      case 1:
        wavetableOneSource = &imValues;
        break;
      case 2:
        wavetableOneSource = &probValues;
        break;
      default:
        break;
    }
  });

  sourceTwo.setElements({"Real Values", "Imaginary Values", "Probability Values"});
  sourceTwo.registerChangeCallback([&](int x) {
    switch (x) {
      case 0:
        wavetableTwoSource = &reValues;
        break;
      case 1:
        wavetableTwoSource = &imValues;
        break;
      case 2:
        wavetableTwoSource = &probValues;
        break;
      default:
        break;
    }
  });

  std::cout << "onInit() - All domains have been initialized " << std::endl;
}

void QHOS::onCreate() {
  navControl().useMouse(false);
  nav().pos(Vec3f(1, 1, 12));
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  nav().setHome();
  waveFunctionPlot.primitive(Mesh::LINE_STRIP);
  probabilityPlot.primitive(Mesh::LINE_STRIP);
  axes.primitive(Mesh::LINES);
  axes.vertex(-5, 0, 0);
  axes.vertex(5, 0, 0);
  axes.vertex(0, 5, 0);
  axes.vertex(0, -5, 0);
  axes.vertex(0, 0, 5);
  axes.vertex(0, 0, -5);

  psiValues.resize(resolution);
  imValues.resize(resolution);
  reValues.resize(resolution);
  probValues.resize(resolution);
  wavetable.resize(2);
  wavetable[0].resize(resolution);
  wavetable[1].resize(resolution);
  std::cout << "onCreate() - Graphics context now available" << std::endl;
}

void QHOS::onAnimate(double dt) {
  simTime += dt * simSpeed;
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  waveFunctionPlot.reset();
  probabilityPlot.reset();
  for (int i = 0; i < resolution; i++) {
    psiValues[i] = psi.evaluate(posValues[i], simTime);
    reValues[i] = psiValues[i].real();
    imValues[i] = psiValues[i].imag();
    probValues[i] = pow(abs(psiValues[i]), 2);
    waveFunctionPlot.vertex(posValues[i], reValues[i], imValues[i]);
    probabilityPlot.vertex(posValues[i], probValues[i], 0);
  }
}

void QHOS::onDraw(Graphics& g) {
  g.clear();

  g.lighting(false);
  g.color(HSV(0.0, 0.0, 1));
  g.draw(axes);
  g.color(HSV(1.0, 1.0, 1));
  g.draw(waveFunctionPlot);
  g.color(HSV(0.6, 1.0, 1));
  g.draw(probabilityPlot);

  if (drawGUI) {
    imguiBeginFrame();

    ParameterGUI::beginPanel("Simulation");
    ImGui::Text("Simulation Time: %.2f", simTime);
    ParameterGUI::drawParameter(&simSpeed);
    ParameterGUI::drawParameterInt(&dims, "");
    ParameterGUI::drawParameterBool(&coeffList, "");
    // Option to manually change each coefficient via slider
    if (coeffList) {
      for (int i = 0; i < manualCoefficients.size(); i++) {
        std::string coeffName = "C_n " + std::to_string(i);
        std::cout << coeffName << " = " << psi.coefficients[i] << std::endl;
        if (SliderDouble((coeffName).c_str(), &manualCoefficients[i], -10, 10)) {
          psi.newHilbertSpace(new HilbertSpace(dims), NULL, manualCoefficients);
        }
      }
    } else {
      ParameterGUI::drawMenu(&presetFuncs);
      if (ImGui::CollapsingHeader("Coefficient Values"))
        for (int i = 0; i < psi.coefficients.size(); i++) {
          ImGui::Text("coefficient %i: %.10f ", i, psi.coefficients[i]);
        }
    }
    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("Audio");
    ParameterGUI::drawMenu(&sourceOne);
    ParameterGUI::drawMenu(&sourceTwo);

    ParameterGUI::endPanel();

    imguiEndFrame();

    imguiDraw();
  }
}

void QHOS::onSound(al::AudioIOData& io) {
  // This is the sample loop
  while (io()) {
    tableReader++;
    if (tableReader == resolution) {
      tableReader = 0;
      wavetable[0] = *wavetableOneSource;
      wavetable[1] = *wavetableTwoSource;
    }

    float out1 = wavetable[0][tableReader];
    float out2 = wavetable[1][tableReader];

    io.out(0) = out1;
    io.out(1) = out2;
  }
}

bool QHOS::onKeyDown(Keyboard const& k) {
  switch (k.key()) {
    case 'g':
      drawGUI = 1 - drawGUI;
      break;
    case 'r':
      nav().home();
      break;
    default:
      break;
  }
  return true;
}

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();

  QHOS app;
  app.configureAudio(dev, dev.defaultSampleRate(), 512, dev.channelsOutMax(), dev.channelsInMax());
  app.start();
  return 0;
}