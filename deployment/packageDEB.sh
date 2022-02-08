#!/bin/bash
result=${PWD##*/}
if [ $result == "deployment" ]; then
  cd ..
fi

mkdir -p "deployment/Linux/"

if [ $# -eq 0 ]; then
  echo "Error: No version number provided. Version number required in format <MajorVersion>.<MinorVersion>"
  exit 1
fi

VERSION="$1""-1"
RELEASENAME="QHOSYN-""$VERSION""-amd64"
BUILDLOCATION=$(cd deployment && pwd)

echo "Packaging $RELEASENAME..."

# make directory structure

BUILDDIR="deployment/Linux/$RELEASENAME"

mkdir -p "$BUILDDIR/DEBIAN"
mkdir -p "$BUILDDIR/usr/bin"
mkdir -p "$BUILDDIR/usr/share/applications/"
mkdir -p "$BUILDDIR/usr/share/doc/qhosyn/"
mkdir -p "$BUILDDIR/usr/share/pixmaps/"
mkdir -p "$BUILDDIR/usr/share/qhosyn/"

# copy necessary files over

cd "$Dir"

cp -r "deployment/externalResources/fonts/" "$BUILDDIR/usr/share/qhosyn/"

objcopy --strip-debug --strip-unneeded bin/QHOSYN "$BUILDDIR/usr/bin/qhosyn"

cp "deployment/icons/QHOSYN.png" "$BUILDDIR/usr/share/pixmaps/QHOSYN.png"

# Make .desktop file
echo "[Desktop Entry]
Name=QHOSYN
Comment=Launch QHOSYN
Exec=QHOSYN
Icon=/usr/share/pixmaps/QHOSYN.png
Terminal=false
Type=Application
Categories=Audio;Music;Science;
Name[en_US]=QHOSYN" >>"$BUILDDIR/usr/share/applications/QHOSYN.desktop"

# make Debian control file
echo "Package: QHOSYN
Architecture: amd64
Section: sound
Priority: optional
Version:$VERSION
Maintainer:Rodney DuPlessis <rodney@rodneyduplessis.com>
Depends:libsndfile1, libc6
Homepage: https://github.com/rodneydup/QHOSYN
Description: This package provides QHOSYN, a synthesizer that sonifies a quantum harmonic oscillator." >>"$BUILDDIR/DEBIAN/control"

# make copyright file
echo "Format: https://www.debian.org/doc/packaging-manuals/copyright-format/1.0/
Upstream-Name: QHOSYN
Source: https://github.com/rodneydup/QHOSYN

Files: *
Copyright: 2020 Rodney DuPlessis <rodney@rodneyduplessis.com>
License: GPL-3+
 This program is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later
 version.
 .
 This program is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the GNU General Public License for more
 details.
 .
 You should have received a copy of the GNU General Public
 License along with this package; if not, write to the Free
 Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 Boston, MA  02110-1301 USA
 .
 On Debian systems, the full text of the GNU General Public
 License version 3 can be found in the file
 '/usr/share/common-licenses/GPL-3'" >>"$BUILDDIR/usr/share/doc/QHOSYN/copyright"
DATE="$(date +'%a, %d %b %Y %H:%M:%S %Z')"
echo "QHOSYN ($VERSION) stable; urgency=high
  * Initial Release
 -- Rodney DuPlessis <rodney@rodneyduplessis.com>  $DATE" >>"$BUILDDIR/usr/share/doc/QHOSYN/changelog.Debian"
gzip -9 -n "$BUILDDIR/usr/share/doc/QHOSYN/changelog.Debian"
echo "Packaging .deb at $BUILDLOCATION..."

# package .deb
cd deployment/Linux
fakeroot dpkg -b "$RELEASENAME" "$RELEASENAME.deb"

echo "Packaging Complete!"
