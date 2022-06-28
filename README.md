# QHOSYN

QHOSYN (Pronounced "cosine") is the Quantum Harmonic Oscillator Synthesizer.

![](./QHOSYN_demo.gif)

QHOSYN is a synthesizer based on the fundamental model of a quantum harmonic oscillator. 

The software employs an accurate quantum simulation of an observable in a superposition of eigenstates up to the n=15 energy level to generate sound from  a quantum wavefunction based on the time-dependent Schrödinger equation.

It synthesizes sound from the simulation in 4 main ways. The first (Wavetable) treats the wavefunction as an evolving wavetable source for wavetable synthesis. The second (Probability Noise Band Synthesis) takes the probability curve of the wavefunction as a spectral curve to create colored noise. The third (IFFT) applies an inverse fourier transform to the wavefuction, treating it as a dynamic spectrum. The fourth (Panner) can be used to generate evolving stochastic sound and control signals from the probability distribution of the wavefunction. Through sonic exploration, the user can generate "quantum sounds" and gain an aesthetic understanding of quantum space.

QHOSYN was designed as part of my PhD research in order to facilitate the composition of my piece _Psi_. It is still in a rough state that was useful for me in my research, but may be confusing and buggy for anyone else to use. Version 1, which will be more user-friendly, will be released soon.

# Building

## Building
### Debian

#### First, check the releases page to see if the pre-compiled installers work for you: https://github.com/rodneydup/QHOSYN/releases/latest

- This project uses cmake to build so if you don't have cmake then install it (Minimum version: 3.13) and make sure your c and c++ compilers are defined in your environment.

- Run the following in a terminal to install the necessary libraries for building:

`sudo apt install libgtk-3-dev libasound2-dev libsndfile1-dev libfftw3-dev libjack-dev`
 
- git clone the repository and run some the scripts that automate the configure/build/install process.

`git clone https://github.com/rodneydup/QHOSYN.git`

`cd QHOSYN`

`./configure.sh`

`./build-run.sh`

### Arch Linux

- Must have cmake (3.13 or later), git, and base-devel packages installed.

`git clone https://github.com/rodneydup/QHOSYN.git`

`cd QHOSYN`

`./configure.sh`

`./build.sh`

### OS X

#### First, check the releases page to see if the pre-compiled installers work for you: https://github.com/rodneydup/QHOSYN/releases/latest

You must have cmake installed (version 3.10 or later), and Xcode (hopefully we can get rid of this dependency soon).

You'll also need a few libraries installed: libsndfile flac libogg libvorbis boost (If you have homebrew, you can get them all with `brew install libsndfile flac libogg libvorbis boost`)

- First, clone the repo:

`git clone https://github.com/rodneydup/QHOSYN.git`

- cd into QHOSYN

`cd QHOSYN`

- run configure script:

`./configure.sh`

- run build script:

`./build-run.sh`

- QHOSYN is located in `QHOSYN/bin/`

### Windows (Visual Studio)

#### First, check the releases page to see if the pre-compiled installers work for you: https://github.com/rodneydup/QHOSYN/releases/latest

For installation through Visual Studio (taken from the allolib repo: https://github.com/AlloSphere-Research-Group/allolib):

    Install Visual Studio 2017 Community Edition from https://visualstudio.microsoft.com/downloads/

    During installation options:

    a. Install "Desktop development with C++" workload

    b. Install Individual components: C++/CLI support, Git for Windows, Visual C++ Tools for Cmake

Install libsndfile: Aim your browser at http://www.mega-nerd.com/libsndfile/#Download. Download and install the 64-bit version: libsndfile-1.0.xx-w64-setup.exe.
    - Install this at the default location.

#### Compile QHOSYN

- You'll have to do some work to get GSL linked... not easy. I recommend using vcpkg.

- Note that you need Visual Studio installed for the below scripts to work.

- Clone the repo using the git bash terminal:

`git clone https://github.com/rodneydup/QHOSYN.git`

- Open Visual Studio pointed at the directory `QHOSYN\`

- In Visual Studio, open the Developer Power Shell.

- When first running the sh scripts, Visual Studio will ask what program would you like to open this file with. Choose `Git Bash`.

- Configure the project:

`.\configure.sh`

- The above step will initially fail. This reason for this is that NFD is using an old version of Visual Studio (2010). To fix this problem, do the following.
     - In file explorer go to the folder `QHOSYN\external\nativefiledialog\build\vs2010\`
     - Open `NativeFileDialog.sln` in Visual Studio.
     - Double click on `Solution 'NativeFileDialog' (5 of 5 projects)` in the vertical window on the right side.
     - Agree to updating this project to the newest toolchain.
     - Run `.\configure.sh` again from the Emission Control 2 Visual Studio Project.
     - Note that the following will only have to be done once per git clone.
     - Sidenote: If you have a better, more stable way, let me know -- jack

- Compile the project
`.\build-run.sh`