#!/bin/bash

rm -rf "deployment/OSX/"
mkdir -p "deployment/OSX/"

cp -r bin/QHOSYN.app deployment/OSX
cd deployment
hdiutil create -volname QHOSYN.app -srcfolder ../deployment/OSX -ov QHOSYN.dmg
mv QHOSYN.dmg OSX/
cd OSX/

echo "QHOSYN.dmg created in: $(pwd)"
