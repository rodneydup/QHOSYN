#!/bin/bash

rm -rf build
rm -rf bin

# currently not purging bin directory in case you need dynamic libraries or other config files with the binary
# rm -rf bin

(
  cd allolib
  ./distclean
)
