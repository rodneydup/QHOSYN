#include "QHOS.hpp"

void QHOS::onInit() {
  std::cout << "onInit() - All domains have been initialized " << std::endl;
  // hilbertSpace hs(2, [](double x) -> double { return 1 / 2 * pow(x, 2); });
  // HilbertSpace hs(2);
}

void QHOS::onCreate() {
  std::cout << "onCreate() - Graphics context now available" << std::endl;
  addWireBox(plot);
}

void QHOS::onAnimate(double dt) {}

void QHOS::onDraw(Graphics& g) {
  g.clear();
  // Note: we don't need to do all the normal graphics setup as this
  // is handled by the App
  // We can just draw our geometry immediately!

  g.polygonMode(Graphics::LINE);  // wireframe mode
  g.pushMatrix();
  g.color(1);
  g.draw(plot);
  g.popMatrix();
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