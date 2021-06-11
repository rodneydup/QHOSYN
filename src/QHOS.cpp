
/*
Quantum Harmonic Oscillator Sonifier
By Rodney DuPlessis
To do:
- make new wavefunction generation asynchronous so it doesn't freeze things up
- find more interesting-sounding situations
- add arbitrary function input option


*/
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOS.hpp"

void QHOS::onInit() {
  imguiInit();
  dims.registerChangeCallback([&](int x) {
    manualCoefficients.resize(dims);
    psi.newHilbertSpace(new HilbertSpace(x), initWaveFunction);
    manualCoefficients = psi.coefficients;
    for (int i = 0; i < manualCoefficients.size(); i++) {
      std::string coeffName = "C_n " + std::to_string(i);
      std::cout << coeffName << " = " << psi.coefficients[i] << std::endl;
    }
  });
  coeffList.registerChangeCallback([&](bool x) {
    if (x) {
      psi.newHilbertSpace(new HilbertSpace(dims), NULL, manualCoefficients);
    } else {
      psi.newHilbertSpace(new HilbertSpace(dims), initWaveFunction);
    }
    for (int i = 0; i < manualCoefficients.size(); i++) {
      std::string coeffName = "C_n " + std::to_string(i);
      std::cout << coeffName << " = " << psi.coefficients[i] << std::endl;
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

  oscClient.open(oscPort, oscAddr);

  std::cout << "onInit() - All domains have been initialized " << std::endl;
}

void QHOS::onCreate() {
  navControl().useMouse(false);
  nav().pos(Vec3f(1, 1, 12));
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  nav().setHome();

  // port, address, timeout
  // "" as address for localhost
  oscServer.open(16448, "localhost", 0.05);

  // Register ourself (osc::PacketHandler) with the server so onMessage
  // gets called.
  oscServer.handler(oscDomain()->handler());

  // Start a thread to handle incoming packets
  oscServer.start();

  waveFunctionPlot.primitive(Mesh::LINE_STRIP);
  probabilityPlot.primitive(Mesh::LINE_STRIP);
  axes.primitive(Mesh::LINES);
  axes.vertex(-5, 0, 0);
  axes.vertex(5, 0, 0);
  axes.vertex(0, 5, 0);
  axes.vertex(0, -5, 0);
  axes.vertex(0, 0, 5);
  axes.vertex(0, 0, -5);
  axes.generateNormals();

  grid.primitive(Mesh::LINES);
  for (int i = -5; i <= 5; i++) {
    // back grid
    grid.vertex(i, -5, -5);
    grid.vertex(i, 5, -5);
    grid.vertex(-5, i, -5);
    grid.vertex(5, i, -5);
    // left grid
    grid.vertex(-5, -5, i);
    grid.vertex(-5, 5, i);
    grid.vertex(-5, i, -5);
    grid.vertex(-5, i, 5);
    // bottom grid
    grid.vertex(-5, -5, i);
    grid.vertex(5, -5, i);
    grid.vertex(i, -5, -5);
    grid.vertex(i, -5, 5);

    for (int i = 0; i < 12; i++) {
      grid.texCoord(0.8, 0.0);
      grid.color(1.0, 1.0, 1.0, 0.2);  // add color for each vertex
    }
  }

  samples.primitive(Mesh::POINTS);

  // compile and link the shaders
  //
  pointShader.compile(shader::pointVertex, shader::pointFragment, shader::pointGeometry);
  lineShader.compile(shader::lineVertex, shader::lineFragment, shader::lineGeometry);

  // use a texture to control the alpha channel of each particle
  //
  pointTexture.create2D(256, 256, Texture::R8, Texture::RED, Texture::SHORT);
  int Nx = pointTexture.width();
  int Ny = pointTexture.height();
  std::vector<short> alpha;
  alpha.resize(Nx * Ny);
  for (int j = 0; j < Ny; ++j) {
    float y = float(j) / (Ny - 1) * 2 - 1;
    for (int i = 0; i < Nx; ++i) {
      float x = float(i) / (Nx - 1) * 2 - 1;
      float m = exp(-13 * (x * x + y * y));
      m *= pow(2, 15) - 1;  // scale by the largest positive short int
      alpha[j * Nx + i] = m;
    }
  }
  pointTexture.submit(&alpha[0]);

  lineTexture.create1D(256, Texture::R8, Texture::RED, Texture::SHORT);
  std::vector<short> beta;
  beta.resize(lineTexture.width());
  for (int i = 0; i < beta.size(); ++i) {
    beta[i] = alpha[128 * beta.size() + i];
  }
  lineTexture.submit(&beta[0]);

  for (int i = 0; i < 6; i++) {
    axes.texCoord(0.8, 0.0);
    axes.color(1.0, 1.0, 1.0, 1.0);  // add color for each vertex
  }

  for (int i = 0; i < resolution; i++) {
    waveFunctionPlot.texCoord(1.0, 0.0);
    probabilityPlot.texCoord(1.0, 0.0);
    waveFunctionPlot.color(1.0, 0.0, 0.0, 0.9);  // add color for each vertex
    probabilityPlot.color(0.0, 0.0, 1.0, 0.9);   // add color for each vertex
  }
  for (int i = 0; i < 20; i++) {
    samples.color(1.0, 1.0, 1.0, 1.0);
    samples.texCoord(1.0, 0.0);
  }

  std::cout << "onCreate() - Graphics context now available" << std::endl;
}

void QHOS::onAnimate(double dt) {
  simTime += dt * simSpeed;
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  waveFunctionPlot.vertices().clear();
  probabilityPlot.vertices().clear();
  samples.vertices().clear();
  wavetableLock.lock();
  for (int i = 0; i < resolution; i++) {
    psiValues[i] = psi.evaluate(posValues[i], simTime);
    reValues[i] = psiValues[i].real();
    imValues[i] = psiValues[i].imag();
    probValues[i] = pow(abs(psiValues[i]), 2);
    waveFunctionPlot.vertex(posValues[i], reValues[i], imValues[i]);
    probabilityPlot.vertex(posValues[i], probValues[i], 0);
  }
  wavetableLock.unlock();
  std::discrete_distribution<int> measurement(probValues.begin(), probValues.end());
  int measurementPoint = measurement(generator);
  samplePoints[sampleCounter] = Vec3d(posValues[measurementPoint], 0, 0);
  for (int i = 0; i < 20; i++) {
    samples.vertex(samplePoints[i]);
  }
  if (oscOn) {
    if (oscWaveform) {
      osc::Packet p;
      p.beginMessage("/realValues/firstHalf");
      for (int i = 0; i < 128; i++) p << (float)reValues[i];
      p.endMessage();
      oscClient.send(p);
      p.clear();
      p.beginMessage("/realValues/secondHalf");
      for (int i = 128; i < 256; i++) p << (float)reValues[i];
      p.endMessage();
      oscClient.send(p);

      p.clear();
      p.beginMessage("/imaginaryValues/firstHalf");
      for (int i = 0; i < 128; i++) p << (float)imValues[i];
      p.endMessage();
      oscClient.send(p);
      p.clear();
      p.beginMessage("/imaginaryValues/secondHalf");
      for (int i = 128; i < 256; i++) p << (float)imValues[i];
      p.endMessage();
      oscClient.send(p);
    }
    if (oscMeasurement) {
      if (measurementTrigger) {
        oscClient.send("/measurement", (float)posValues[measurementPoint]);
        measurementTrigger = 0;
      }
    }
  }
  sampleCounter = (sampleCounter + 1) % 20;
}

void QHOS::onDraw(Graphics& g) {
  g.clear();

  gl::blending(true);  // needed for transparency
  gl::blendTrans();    // needed for transparency

  lineTexture.bind();  // texture binding
  g.meshColor();
  g.shader(lineShader);  // run shader
  if (drawAxes) g.draw(axes);
  if (drawGrid) g.draw(grid);
  if (drawWavefunction) g.draw(waveFunctionPlot);
  if (drawProbability) g.draw(probabilityPlot);
  lineTexture.unbind();

  pointTexture.bind();
  g.shader(pointShader);
  if (drawMeasurements) g.draw(samples);
  pointTexture.unbind();

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
    ParameterGUI::drawParameterBool(&audioOn);
    ParameterGUI::drawParameter(&freq);
    ParameterGUI::drawParameter(&volume);
    ParameterGUI::drawMenu(&sourceOne);
    ParameterGUI::drawMenu(&sourceTwo);

    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("Draw");
    ParameterGUI::drawParameterBool(&drawGrid);
    ParameterGUI::drawParameterBool(&drawAxes);
    ParameterGUI::drawParameterBool(&drawWavefunction);
    ParameterGUI::drawParameterBool(&drawProbability);
    ParameterGUI::drawParameterBool(&drawMeasurements);

    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("OSC");
    ParameterGUI::drawParameterBool(&oscOn);
    ParameterGUI::drawParameterBool(&oscWaveform);
    ParameterGUI::drawParameterBool(&oscMeasurement);
    ParameterGUI::drawParameterBool(&measurementTrigger);

    ParameterGUI::endPanel();

    imguiEndFrame();

    imguiDraw();
  }
}

void QHOS::onSound(al::AudioIOData& io) {
  // This is the sample loop
  while (io()) {
    if (audioOn) {
      // increment table reader
      tableReader += freq / (io.fps() / resolution);
      // if table reader gets to the end of the table
      if (tableReader >= resolution) {
        // loop
        tableReader -= resolution;
        // update table
        if (wavetableLock.try_lock()) {
          wavetable[0] = *wavetableOneSource;
          wavetable[1] = *wavetableTwoSource;
          wavetableLock.unlock();
        }
      }
      // output sample
      float sample[2];
      // linear interpolation for table reads between indices
      for (int i = 0; i < 2; i++) {
        int j = floor(tableReader);
        float x0 = wavetable[i][j];
        float x1 = wavetable[i][(j == (wavetable[i].size() - 1)) ? 0 : j + 1];  // looping semantics
        float t = tableReader - j;
        sample[i] = (x1 * t + x0 * (1 - t)) / 2;  // (divided by 2 because it's loud)
      }
      // output
      io.out(0) = sample[0] * volume;
      io.out(1) = sample[1] * volume;
    }
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

void QHOS::onResize(int w, int h) {
  updateFBO(width(), height());
  //
}

void QHOS::onMessage(osc::Message& m) {
  // Check that the address and tags match what we expect
  if (m.addressPattern() == "/measure" && m.typeTags() == "si") {
    // Extract the data out of the packet
    std::string str;
    int val;
    m >> str >> val;

    if (val == 1) measurementTrigger = 1;
  }
}

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();

  QHOS app;
  app.configureAudio(dev, dev.defaultSampleRate(), 512, dev.channelsOutMax(), dev.channelsInMax());
  app.start();
  return 0;
}