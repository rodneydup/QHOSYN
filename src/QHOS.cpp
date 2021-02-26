#include <iostream>

// for master branch
// #include "al/core.hpp"

// for devel branch
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"

using namespace al;

struct MyApp : App {
  void onCreate() override {}

  void onAnimate(double dt) override {}

  void onDraw(Graphics& g) override {}
};

int main() {
  MyApp app;
  app.dimensions(600, 400);
  app.start();
}