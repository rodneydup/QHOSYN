#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOS.hpp"

void QHOS::onInit() {
  std::cout << "onInit() - All domains have been initialized " << std::endl;
  // hilbertSpace hs(2, [](double x) -> double { return 1 / 2 * pow(x, 2); });
}

void QHOS::onCreate() {
  std::cout << "onCreate() - Graphics context now available" << std::endl;
  nav().pos(Vec3f(1, 1, -8));
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
  // plot.reset();

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
}

void QHOS::onSound(al::AudioIOData& io) {
  // This is the sample loop
  while (io()) {
    // float in = io.in(0);

    float out1 = 0;
    float out2 = 0;

    io.out(0) = out1;
    io.out(1) = out2;
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