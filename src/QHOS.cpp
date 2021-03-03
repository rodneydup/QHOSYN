#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "QHOS.hpp"

void QHOS::onInit() {
  std::cout << "onInit() - All domains have been initialized " << std::endl;
  // hilbertSpace hs(2, [](double x) -> double { return 1 / 2 * pow(x, 2); });
}

void QHOS::onCreate() {
  std::cout << "onCreate() - Graphics context now available" << std::endl;
  nav().pos(Vec3f(0, 0, -2));
  plot.primitive(Mesh::LINE_STRIP);
  plot.vertex(-1, 0, 0);
  plot.vertex(1, 0, 0);

  amplitudeValues.resize(posValues.size());
}

void QHOS::onAnimate(double dt) {
  nav().faceToward(Vec3f(0.0, 0.0, 0.0));
  plot.reset();

  for (int i = 0; i < posValues.size(); i++) amplitudeValues[i] = psi.evaluate(i, i / 50);

  std::vector<double> reValues;
  std::vector<double> imValues;
}

void QHOS::onDraw(Graphics& g) {
  g.clear();

  g.lighting(false);
  g.color(HSV(0.0, 0.5, 1));
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