
/*
Quantum Harmonic Oscillator Synthesizer (QHOSYN)
By Rodney DuPlessis (2021)
*/

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOSYN.hpp"

void QHOSYN::onInit() {
  title("QHOSYN");
  audioIO().setStreamName("QHOSYN");

  execDir = al::File::directory(getExecutablePath());

#ifdef __APPLE__
  execDir = getContentPath_OSX(execDir);
#endif

  audioIO().deviceOut(AudioDevice::defaultOutput());
  setAudioSettings(1);

  opener = "xdg-open ";
#ifdef __APPLE__
  opener = "open ";
#endif

#ifdef __linux__
  opener = "xdg-open ";
#endif

#ifdef _WIN32
  opener = "start ";
#endif

  audioIO().append(mRecorder);
  srand(std::time(0));
  emptyTable.fill(0);

  imguiInit();
  ImGui::GetIO().IniFilename = NULL;

  // Set up parameters
  dims.registerChangeCallback([&](int x) {
    manualCoefficients.resize(x);
    if (coeffList)
      psi.newHilbertSpace(&basis, x, NULL, manualCoefficients);
    else {
      psi.newHilbertSpace(&basis, x, initWaveFunction, {}, project);
    }
    manualCoefficients = psi.coefficients;

    for (int i = 0; i < manualCoefficients.size(); i++) {
      std::string coeffName = "C_n " + std::to_string(i);
      std::cout << coeffName << " = " << psi.coefficients[i] << std::endl;
    }
  });

  coeffList.registerChangeCallback([&](bool x) {
    if (x) {
      psi.newHilbertSpace(&basis, dims, NULL, manualCoefficients);
    } else {
      bool project = 1;
      if (presetFunctions.get() == 1) project = 0;
      psi.newHilbertSpace(&basis, dims, initWaveFunction, {}, project);
    }
    for (int i = 0; i < manualCoefficients.size(); i++) {
      std::string coeffName = "C_n " + std::to_string(i);
      std::cout << coeffName << " = " << psi.coefficients[i] << std::endl;
    }
  });

  presetFunctions.setElements({"ϕ(x) = sin(x)", "ϕ(x) = if(x=MaxEigenstate, 1, 0)",
                               "ϕ(x) = rand(0, 1)", "ϕ(x) = x - 2", "ϕ(x) = 1/(x+1)", "Custom"});
  presetFunctions.registerChangeCallback([&](int choice) {
    switch (choice) {
      case 0:
        initWaveFunction = [](double x) { return sin(x); };
        break;
      case 1:
        initWaveFunction = [&](double x) {
          if (abs(x) == 10)
            return 1.0;
          else
            return 0.0;
        };
        break;
      case 2:
        initWaveFunction = [](double x) { return double(rand()) / (RAND_MAX / 2) - 1; };
        break;
      case 3:
        initWaveFunction = [](double x) { return x - 2; };
        break;
      case 4:
        initWaveFunction = [](double x) { return 1.0 / (x + 1); };
        break;
      case 5:
        break;
      default:
        break;
    }
    if (presetFunctions.get() != 5)
      psi.newHilbertSpace(&basis, dims, initWaveFunction, {}, project);
  });

  project.registerChangeCallback([&](bool on) {
    if (on)
      psi.newHilbertSpace(&basis, dims, initWaveFunction, {}, 1);
    else
      psi.newHilbertSpace(&basis, dims, initWaveFunction, {}, 0);
  });

  sourceOneMenu.setElements({"none", "Real Wavetable", "Imaginary Wavetable",
                             "Probability Wavetable", "Probability Noise Band",
                             "Inverse Fourier Transform"});
  sourceOneMenu.registerChangeCallback([&](int x) {
    sourceSelect[0] = 0;
    switch (x) {
      case 0:
        source[0] = &emptyTable;
        break;
      case 1:
        source[0] = &reValues;
        break;
      case 2:
        source[0] = &imValues;
        break;
      case 3:
        source[0] = &probValues;
        break;
      case 4:
        sourceSelect[0] = 4;
        break;
      case 5:
        sourceSelect[0] = 5;
        break;
      default:
        break;
    }
  });
  sourceOneMenu.set(1);

  sourceTwoMenu.setElements({"none", "Real Wavetable", "Imaginary Wavetable",
                             "Probability Wavetable", "Probability Noise Band",
                             "Inverse Fourier Transform"});
  sourceTwoMenu.registerChangeCallback([&](int x) {
    sourceSelect[1] = 0;
    switch (x) {
      case 0:
        source[1] = &emptyTable;
        break;
      case 1:
        source[1] = &reValues;
        break;
      case 2:
        source[1] = &imValues;
        break;
      case 3:
        source[1] = &probValues;
        break;
      case 4:
        sourceSelect[1] = 4;
        break;
      case 5:
        sourceSelect[1] = 5;
        break;
      default:
        break;
    }
  });
  sourceTwoMenu.set(2);

  // Check for connected MIDI devices
  if (midiIn.getPortCount() > 0) {
    // Bind ourself to the RtMidiIn object, to have the onMidiMessage()
    // callback called whenever a MIDI message is received
    MIDIMessageHandler::bindTo(midiIn);

    // Open the last device found
    unsigned int port = midiIn.getPortCount() - 1;
    midiIn.openPort(port);
    printf("Opened port to %s\n", midiIn.getPortName(port).c_str());
  } else {
    printf("Warning: No MIDI devices found. Connect a MIDI device and restart the application for MIDI support.\n");
  }

  std::cout << "onInit() - All domains have been initialized " << std::endl;
}

void QHOSYN::onCreate() {
  // Camera setup
  navControl().useMouse(false);
  nav().pos(Vec3f(1, 1, 12));
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  nav().setHome();

  // OSC setup
  oscSender = std::make_unique<osc::Send>();
  resetOSC();

  ImFontConfig fontConfig;
  fontConfig.OversampleH = 4;
  fontConfig.OversampleV = 4;

  ImGuiIO& io = ImGui::GetIO();

  // build a range to inlcude special characters in font
  ImVector<ImWchar> ranges;
  ImFontGlyphRangesBuilder builder;
  builder.AddRanges(io.Fonts->GetGlyphRangesDefault());
  static const ImWchar icons_ranges[] = {0x0370, 0x03FF, 0};  // add greek characters
  builder.AddRanges(icons_ranges);
  builder.BuildRanges(
    &ranges);  // Build the final result (ordered ranges with all the unique characters submitted)

#ifdef __APPLE__
  bodyFont = io.Fonts->AddFontFromFileTTF((execDir + "Resources/fonts/Roboto-Regular.ttf").c_str(),
                                          16.0f, &fontConfig, ranges.Data);
  titleFont = io.Fonts->AddFontFromFileTTF((execDir + "Resources/fonts/Roboto-Regular.ttf").c_str(),
                                           20.0f, &fontConfig, ranges.Data);
#endif

#ifdef __linux__
  bodyFont = io.Fonts->AddFontFromFileTTF("/usr/share/qhosyn/fonts/Roboto-Regular.ttf", 16.0f,
                                          &fontConfig, ranges.Data);
  titleFont = io.Fonts->AddFontFromFileTTF("/usr/share/qhosyn/fonts/Roboto-Regular.ttf", 20.0f,
                                           &fontConfig, ranges.Data);
#endif

#ifdef _WIN32
  bodyFont = io.Fonts->AddFontFromFileTTF((execDir + "Resources/fonts/Roboto-Regular.ttf").c_str(),
                                          16.0f, &fontConfig, ranges.Data);
  titleFont = io.Fonts->AddFontFromFileTTF((execDir + "Resources/fonts/Roboto-Regular.ttf").c_str(),
                                           20.0f, &fontConfig, ranges.Data);
#endif
  io.Fonts->Build();  // Build the atlas while 'ranges' is still in scope and not deleted.

  // set up drawing meshes
  waveFunctionPlot.primitive(Mesh::LINE_STRIP);
  probabilityPlot.primitive(Mesh::LINE_STRIP);
  noisePlot.primitive(Mesh::LINE_STRIP);
  ifftPlot.primitive(Mesh::LINE_STRIP);
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
      grid.color(1.0, 1.0, 1.0, 0.2);
    }
  }
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
  for (int i = 0; i < bufferLen; i++) {
    ifftPlot.texCoord(1.0, 0.0);
    noisePlot.texCoord(1.0, 0.0);
    ifftPlot.color(0.0, 1.0, 1.0, 0.9);
    noisePlot.color(1.0, 1.0, 0.0, 0.9);
  }
  for (int i = 0; i < 20; i++) {
    samples.texCoord(1.0, 0.0);
    samples.color(1.0, 1.0, 1.0, 1.0);
  }

  samples.primitive(Mesh::POINTS);

  // compile and link shaders
  pointShader.compile(shader::pointVertex, shader::pointFragment, shader::pointGeometry);
  lineShader.compile(shader::lineVertex, shader::lineFragment, shader::lineGeometry);

  // Textures for fuzzy lines and points (thanks Karl Yerkes)
  pointTexture.create2D(512, 512, Texture::R8, Texture::RED, Texture::SHORT);
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

  lineTexture.create1D(512, Texture::R8, Texture::RED, Texture::SHORT);
  std::vector<short> beta;
  beta.resize(lineTexture.width());
  for (int i = 0; i < beta.size(); ++i) {
    beta[i] = alpha[256 * beta.size() + i];
  }
  lineTexture.submit(&beta[0]);

  // fill windowing arrays
  gam::tbl::hann(&overlapAddWindow[0], overlapAddWindow.size());
  gam::tbl::hamming(&filterKernelWindow[0], filterKernelWindow.size());

  // inverse fftw
  c2rIn = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ((fftSize / 2) + 1));
  c2rOut = (double*)fftw_malloc(sizeof(double) * fftSize);
  c2r = fftw_plan_dft_c2r_1d(fftSize, c2rIn, c2rOut, FFTW_ESTIMATE);

  // forward fftw
  r2cIn = (double*)fftw_malloc(sizeof(double) * fftSize);
  r2cOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ((fftSize / 2) + 1));
  r2c = fftw_plan_dft_r2c_1d(fftSize, r2cIn, r2cOut, FFTW_ESTIMATE);

  // zero the buffers
  for (int i = 0; i < ifftBuffers.size(); i++)
    for (int j = 0; j < ifftBuffers[i].size(); j++)
      for (int k = 0; k < ifftBuffers[i][j].size(); k++) ifftBuffers[i][j][k] = 0;

  // FFT bin size
  binSize = audioIO().framesPerSecond() / fftSize;

  // Set filters and smoothvalues
  for (int chan = 0; chan < 2; chan++) {
    pan[chan].changeType("lin");
    pan[chan].setTime(16.0f);
    pan[chan].setTarget(chan);
    filter[chan].type(gam::RESONANT);
    filter[chan].res(2);
    antiAlias[chan].type(gam::LOW_PASS);
    antiAlias[chan].res(0.5);
    antiAlias[chan].freq(15000);
  }

  std::cout << "onCreate() - Graphics context now available" << std::endl;
}

void QHOSYN::onAnimate(double dt) {
  simTime += dt * simSpeed;
  if (automeasure) automeasureTimer += dt;
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  waveFunctionPlot.vertices().clear();
  probabilityPlot.vertices().clear();
  samples.vertices().clear();
  ifftPlot.vertices().clear();
  noisePlot.vertices().clear();

  // Update wavefunction
  for (int i = 0; i < resolution; i++) {
    psiValues[i] = psi.evaluate(posValues[i], simTime);
    reValues[i] = psiValues[i].real();
    imValues[i] = psiValues[i].imag();
    probValues[i] = pow(abs(psiValues[i]), 2);
    waveFunctionPlot.vertex(posValues[i], reValues[i], imValues[i]);
    probabilityPlot.vertex(posValues[i], probValues[i], 0);
  }

  // Fourier methods, skipped if not active for efficiency
  for (int chan = 0; chan < 2; chan++)
    if (sourceSelect[chan] >= 4) {
      for (int i = 0; i < fftSize / 2 + 1; i++) {
        c2rIn[i][0] = 0;
        c2rIn[i][1] = 0;
      }
      // Noise band synthesis
      if (sourceSelect[chan] == 4) {
        for (int i = ifftBin; i < ifftBin + ((resolution - 1) * ifftBandwidth); i++) {
          if (i < fftSize / 2) {
            double phase = uniform(gen);
            float t = (i - ifftBin) / (float)ifftBandwidth;
            float temp;
            float value = ((probValues[floor(t)] * (1.0f - modf(t, &temp))) +
                           (probValues[floor(t) + 1] * modf(t, &temp))) /
                          2.0f;
            c2rIn[i][0] = cos(phase) * value / (resolution / 8) / dims;
            c2rIn[i][1] = sin(phase) * value / (resolution / 8) / dims;
          } else {
            break;
          }
        }
        fftw_execute(c2r);
        for (int i = 0; i < bufferLen; i++) {
          ifftBuffers[chan][(currentIfftBuffer[chan] + 1) % 4][i] = c2rOut[i];
        }
        if (noisePlot.vertices().size() == 0) {
          for (int i = 0; i < bufferLen; i++) {
            noisePlot.vertex(((i / (float)bufferLen) * 10.0f) - 5.0f, c2rOut[i], 0);
          }
        }
        // IFFT Synthesis
      } else if (sourceSelect[chan] == 5) {
        for (int i = ifftBin; i < ifftBin + ((resolution - 1) * ifftBandwidth); i++) {
          if (i < fftSize / 2) {
            float t = (i - ifftBin) / (float)ifftBandwidth;
            float temp;
            float reVal = ((reValues[floor(t)] * (1.0f - modf(t, &temp))) +
                           (reValues[floor(t) + 1] * modf(t, &temp))) /
                          2.0f;
            float imVal = ((imValues[floor(t)] * (1.0f - modf(t, &temp))) +
                           (imValues[floor(t) + 1] * modf(t, &temp))) /
                          2.0f;
            c2rIn[i][0] = reVal / (M / 8) / dims;
            c2rIn[i][1] = imVal / (M / 8) / dims;
          }
        }
        fftw_execute(c2r);
        for (int i = 0; i < M; i++) {
          r2cIn[i] = c2rOut[((fftSize - (M / 2)) + i) % fftSize] * filterKernelWindow[i];
        }
        for (int i = M; i < fftSize; i++) {
          r2cIn[i] = 0;
        }
        fftw_execute(r2c);

        c2rIn[0][0] = 0;
        c2rIn[0][1] = 0;
        c2rIn[fftSize / 2][0] = 0;
        c2rIn[fftSize / 2][1] = 0;
        for (int i = 1; i < fftSize / 2; i++) {
          c2rIn[i][0] = r2cOut[i][0] / M;
          c2rIn[i][1] = r2cOut[i][1] / M;
        }
        fftw_execute(c2r);

        for (int i = 0; i < bufferLen; i++) {
          ifftBuffers[chan][(currentIfftBuffer[chan] + 1) % 4][i] = c2rOut[i];
        }
        if (ifftPlot.vertices().size() == 0) {
          for (int i = 0; i < bufferLen; i++) {
            ifftPlot.vertex(((i / (float)bufferLen) * 10.0f) - 5.0f, c2rOut[i], 0);
          }
        }
      }
    }

  // Measurements
  if (automeasureTimer > automeasureInterval) {
    measurementTrigger = 1;
    automeasureTimer -= automeasureInterval;
  }
  if (measurementTrigger) {
    std::discrete_distribution<> measurement(probValues.begin(), probValues.end());
    measurementPoints[measurementPointsCounter] = Vec3d(posValues[measurement(gen)], 0, 0);
    sampleDisplayTimer[measurementPointsCounter] = 20;
    pannerTrigger[0] = 1;
    pannerTrigger[1] = 1;
  }
  for (int i = 0; i < 20; i++) {
    if (sampleDisplayTimer[i] > 0) {
      samples.vertex(measurementPoints[i]);
      sampleDisplayTimer[i] -= 1;
    }
  }

  // OSC sending
  if (oscSenderOn) {
    if (oscWaveform) {
      osc::Packet p;
      p.beginMessage("/realValues/firstHalf");
      for (int i = 0; i < 128; i++) p << (float)reValues[i];
      p.endMessage();
      oscSender->send(p);
      p.clear();
      p.beginMessage("/realValues/secondHalf");
      for (int i = 128; i < 256; i++) p << (float)reValues[i];
      p.endMessage();
      oscSender->send(p);

      p.clear();
      p.beginMessage("/imaginaryValues/firstHalf");
      for (int i = 0; i < 128; i++) p << (float)imValues[i];
      p.endMessage();
      oscSender->send(p);
      p.clear();
      p.beginMessage("/imaginaryValues/secondHalf");
      for (int i = 128; i < 256; i++) p << (float)imValues[i];
      p.endMessage();
      oscSender->send(p);
    }
    if (oscMeasurement) {
      if (measurementTrigger) {
        oscSender->send(measurementArg, (float)measurementPoints[measurementPointsCounter][0]);
      }
    }
  }
  measurementPointsCounter = (measurementPointsCounter + 1) % 20;
  measurementTrigger = 0;
}

void QHOSYN::onDraw(Graphics& g) {
  g.clear();
  // needed for transparency
  gl::blending(true);
  gl::blendMode(GL_SRC_ALPHA, GL_DST_ALPHA, GL_FUNC_ADD);
  g.lighting(true);

  // Draw lines
  lineTexture.bind();
  g.meshColor();
  g.shader(lineShader);
  if (drawAxes) g.draw(axes);
  if (drawGrid) g.draw(grid);
  if (drawWavefunction) g.draw(waveFunctionPlot);
  if (drawProbability) g.draw(probabilityPlot);
  if (drawIfft) g.draw(ifftPlot);
  if (drawNoise) g.draw(noisePlot);
  lineTexture.unbind();

  // Draw measurement points
  pointTexture.bind();
  g.shader(pointShader);
  if (drawMeasurements) g.draw(samples);
  pointTexture.unbind();

  if (drawGUI) {
    imguiBeginFrame();
    int yposition = 0;

    // Simulation Control Panel
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Simulation", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ImGui::Text("Simulation Time: %.2fs", simTime);
    ImGui::PushItemWidth(-120);
    ParameterGUI::drawParameter(&simSpeed);
    ParameterGUI::drawParameterInt(&dims, "");
    ParameterGUI::drawParameterBool(&coeffList, "");
    // Option to manually change each coefficient via slider
    if (coeffList) {
      for (int i = 0; i < manualCoefficients.size(); i++) {
        std::string coeffName = "C_n " + std::to_string(i);
        if (SliderDouble((coeffName).c_str(), &manualCoefficients[i], -10, 10)) {
          psi.newHilbertSpace(&basis, dims, NULL, manualCoefficients);
        }
      }
    } else {
      ParameterGUI::drawMenu(&presetFunctions);
      if (presetFunctions.get() == 5) {
        ImGui::Text("ϕ(x) =");
        ImGui::SameLine();
        if (InputText("##customFunction", &customFunction, ImGuiInputTextFlags_EnterReturnsTrue,
                      inputTextCallback, CallbackUserData)) {
          parser.setFunction(customFunction);
          parser.parse();
          std::cout << "got here" << std::endl;
          initWaveFunction = [this](double x) { return parser.evaluate(x); };
          psi.newHilbertSpace(&basis, dims, initWaveFunction, {}, project);
        }
        if (ImGui::IsItemHovered())
          ImGui::SetTooltip(
            "Enter a custom wavefunction such as x+2\nSupports all standard arithmetic operators "
            "such as multiply (*), divide (/), \nadd(+), subtract (-), modulo(%), and "
            "exponentiation"
            "(^).\nAlso supports rand(min,max) for uniform random number (float) generation, \n"
            "trigonometric functions such as sin(x), and many others.\nClick the Help button "
            "below for more options.");
        if (ImGui::Button("Help")) {
          system((opener + functionParserURL).c_str());
        }
      }

      ParameterGUI::drawParameterBool(&project);
      if (ImGui::CollapsingHeader("Coefficient Values"))
        for (int i = 0; i < psi.coefficients.size(); i++) {
          ImGui::Text("coefficient %i: %.10f ", i, psi.coefficients[i]);
        }
      ImGui::PopItemWidth();
    }
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    // Audio Control Panel
    ImGui::PopFont();
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Audio", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ImGui::PushItemWidth(-120);
    ParameterGUI::drawParameterBool(&audioOn);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&panner);
    ParameterGUI::drawParameter(&volume);
    ParameterGUI::drawMenu(&sourceOneMenu);
    ParameterGUI::drawMenu(&sourceTwoMenu);
    if (!sourceSelect[0] || !sourceSelect[1]) {
      ImGui::Text("Wavetable Synthesis");
      ParameterGUI::drawParameter(&wavetableFreq);
    }
    if (sourceSelect[0] || sourceSelect[1]) {
      ImGui::Text("Fourier Synthesis");
      ParameterGUI::drawParameterInt(&ifftBin, "");
      ParameterGUI::drawParameterInt(&ifftBandwidth, "");
    }
    if (ImGui::CollapsingHeader("Recorder")) {
      drawRecorderWidget(&mRecorder, audioIO().framesPerSecond(),
                         audioIO().channelsOut() <= 1 ? 1 : 2, soundOutput, BLOCK_SIZE);
      if (ImGui::Button("Change Output Path")) {
        result = NFD_PickFolder(NULL, &outPath);

        if (result == NFD_OKAY) {
          std::string temp = outPath;
          setSoundOutputPath(outPath);
        }
      }
    }
    if (ImGui::CollapsingHeader("Audio IO Settings")) {
      drawAudioIO(&audioIO());
    }
    ImGui::PopItemWidth();
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    // Draw Control Panel
    ImGui::PopFont();
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Draw", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&drawGrid);
    ParameterGUI::drawParameterBool(&drawAxes);
    ParameterGUI::drawParameterBool(&drawWavefunction);
    ParameterGUI::drawParameterBool(&drawProbability);
    ParameterGUI::drawParameterBool(&drawIfft);
    ParameterGUI::drawParameterBool(&drawNoise);
    ParameterGUI::drawParameterBool(&drawMeasurements);
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    // OSC Control Panel
    ImGui::PopFont();
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("OSC", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&oscSenderOn);
    if (oscSenderOn) {
      ImGui::PushItemWidth(200);
      if (InputText("Client IP", &oscSenderAddr, ImGuiInputTextFlags_EnterReturnsTrue,
                    inputTextCallback, CallbackUserData))
        resetOSC();
      if (ImGui::InputInt("Client Port", &oscSenderPort, 1, 1,
                          ImGuiInputTextFlags_EnterReturnsTrue))
        resetOSC();
      ParameterGUI::drawParameterBool(&oscWaveform);
      ParameterGUI::drawParameterBool(&oscMeasurement);
      InputText("Measure Argument", &measurementArg, ImGuiInputTextFlags_EnterReturnsTrue,
                inputTextCallback, CallbackUserData);
      ImGui::PopItemWidth();
    }
    ParameterGUI::drawParameterBool(&oscReceiverOn);
    if (oscReceiverOn) {
      ImGui::PushItemWidth(200);
      if (InputText("Server IP", &oscReceiverAddr, ImGuiInputTextFlags_EnterReturnsTrue,
                    inputTextCallback, CallbackUserData))
        resetOSC();
      if (ImGui::InputInt("Server Port", &oscReceiverPort, 1, 1,
                          ImGuiInputTextFlags_EnterReturnsTrue))
        resetOSC();
      ImGui::PopItemWidth();
    }
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    // Measurement Control Panel
    ImGui::PopFont();
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Measurement", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    measurementTrigger = ImGui::Button("Measure", ImVec2(100, 20));
    ParameterGUI::drawParameterBool(&automeasure);
    ImGui::PushItemWidth(-120);
    ParameterGUI::drawParameter(&automeasureInterval);
    ImGui::PopItemWidth();
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    // OSC warning window if OSC port fails to connect.
    if (isOscWarningWindow) {
      ImGui::OpenPopup("OSC Error");
    }
    bool isOSCWarningOpen = true;
    if (ImGui::BeginPopupModal("OSC Error", &isOSCWarningOpen)) {
      isOscWarningWindow = false;
      ImGui::TextColored(ImVec4(0.9, 0.2, 0.2, 0.9), "Warning");
      ImGui::Text(
        "Could not bind to UDP socket. Is there a server already bound to that "
        "port?");
      ImGui::EndPopup();
    }
    ImGui::PopFont();

    imguiEndFrame();

    imguiDraw();
  }
}

void QHOSYN::onSound(al::AudioIOData& io) {
  // This is the sample loop
  while (io()) {
    if (audioOn) {
      // Wavetable Synthesis
      if (sourceSelect[0] <= 3 || sourceSelect[1] <= 3) {
        filter[0].freq(wavetableFreq * dims);
        filter[1].freq(wavetableFreq * dims);

        // increment table reader
        tableReader += wavetableFreq / (io.fps() / resolution);
        // if table reader gets to the end of the table
        if (tableReader >= resolution) {
          // loop
          tableReader -= resolution;
          for (int chan = 0; chan < 2; chan++) {
            // new panning
            if (panner) {
              if (pannerTrigger[chan]) {
                if (sourceSelect[chan] <= 3) {
                  pan[chan].setCurrentValue((measurementPoints[measurementPointsCounter].x / 10) +
                                            0.5);
                  pannerTrigger[chan] = false;
                }
              }
            } else {
              if (sourceSelect[chan] <= 3) pan[chan].setTarget(chan);
            }

            // update table
            wavetable[chan] = *source[chan];
          }
        }
        // linear interpolation for table reads between indices
        for (int chan = 0; chan < 2; chan++) {
          int j = floor(tableReader);
          float x0 = wavetable[chan][j];
          float x1 =
            wavetable[chan]
                     [(j == (wavetable[chan].size() - 1)) ? 0 : j + 1];  // wrap at end of table
          float t = tableReader - j;
          if (sourceSelect[chan] <= 3)
            if (filterOn)
              sample[chan] =
                filter[chan]((x1 * t) + (x0 * (1 - t))) / 2;  // (divided by 2 because it's loud)
            else
              sample[chan] = ((x1 * t) + (x0 * (1 - t))) / 2.0f;
        }
      }
      // Fourier Synthesis
      for (int chan = 0; chan < 2; chan++) {
        if (sourceSelect[chan] >= 4) {
          if (currentIfftSample[chan] >= bufferLen / 2) {
            currentIfftSample[chan] = 0;
            currentIfftBuffer[chan] = (currentIfftBuffer[chan] + 1) % 4;
          }

          float ifftsample = ifftBuffers[chan][currentIfftBuffer[chan]][currentIfftSample[chan]] *
                               overlapAddWindow[currentIfftSample[chan]] +
                             ifftBuffers[chan][(currentIfftBuffer[chan] + 3) % 4]
                                        [currentIfftSample[chan] + bufferLen / 2] *
                               overlapAddWindow[currentIfftSample[chan] + bufferLen / 2];

          currentIfftSample[chan]++;
          sample[chan] = ifftsample;

          if (panner) {
            if (pannerTrigger[chan])
              if (ifftsample < 0.000001 || ifftsample > -0.000001) {
                pan[chan].setTarget((measurementPoints[measurementPointsCounter].x / 10) + 0.5);
                pannerTrigger[chan] = false;
              }
          } else {
            pan[chan].setTarget(chan);
          }
        }
      }
      for (int chan = 0; chan < 2; chan++) pan[chan].process();
      // output
      // std::cout << pan[0].getCurrentValue() << std::endl;
      io.out(0) = antiAlias[0]((sample[0] * (1.0f - pan[0].getCurrentValue())) +
                               (sample[1] * (1.0f - pan[1].getCurrentValue()))) *
                  volume;
      io.out(1) = antiAlias[1]((sample[0] * pan[0].getCurrentValue()) +
                               (sample[1] * pan[1].getCurrentValue())) *
                  volume;
    }
  }
}

bool QHOSYN::onKeyDown(Keyboard const& k) {
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

void QHOSYN::onResize(int w, int h) {
  updateFBO(width(), height());
  //
}

void QHOSYN::onMessage(osc::Message& m) {
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

  QHOSYN app;
  app.start();
  return 0;
}