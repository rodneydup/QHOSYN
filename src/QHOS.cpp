#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOS.hpp"

void QHOS::onInit() {
  imguiInit();

  std::cout << "onInit() - All domains have been initialized " << std::endl;
  // hilbertSpace hs(2, [](double x) -> double { return 1 / 2 * pow(x, 2); });
}

void QHOS::onCreate() {
  std::cout << "onCreate() - Graphics context now available" << std::endl;
  nav().pos(Vec3f(1, 1, 8));
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  nav().setHome();
  plot.primitive(Mesh::LINE_STRIP);
  axes.primitive(Mesh::LINES);
  axes.vertex(-5, 0, 0);
  axes.vertex(5, 0, 0);
  axes.vertex(0, 5, 0);
  axes.vertex(0, -5, 0);
  axes.vertex(0, 0, 5);
  axes.vertex(0, 0, -5);

  amplitudeValues.resize(posValues.size());
  imValues.resize(posValues.size());
  reValues.resize(posValues.size());
  probValues.resize(posValues.size());
}

void QHOS::onAnimate(double dt) {
  time += dt;
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  plot.reset();

  for (int i = 0; i < posValues.size(); i++) amplitudeValues[i] = psi.evaluate(posValues[i], time);
  for (int i = 0; i < posValues.size(); i++) {
    plot.vertex(posValues[i], amplitudeValues[i].real(), amplitudeValues[i].imag());
    // std::cout << amplitudeValues[i].real() << std::endl;
  }
}

void QHOS::onDraw(Graphics& g) {
  g.clear();

  g.lighting(false);
  g.color(HSV(0.0, 0.0, 1));
  g.draw(axes);
  g.color(HSV(1.0, 1.0, 1));
  g.draw(plot);

  if (drawGUI) {
    imguiBeginFrame();

    ParameterGUI::beginPanel("Simulation");
    if (ImGui::SliderInt("Dimensions", &dims, 2, 10)) {
      psi.newHilbertSpace(new HilbertSpace(dims), [](double x) {
        if (abs(x) < (2 * M_PI))
          return sin(x);
        else
          return 0.0;
      });
    }

    imguiEndFrame();

    imguiDraw();
  }
}

void QHOS::onSound(al::AudioIOData& io) {
  // This is the sample loop
  while (io()) {
    tableReader++;
    tableReader == amplitudeValues.size() ? tableReader = 0 : NULL;

    float out1 = amplitudeValues[tableReader].real();
    float out2 = amplitudeValues[tableReader].real();

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