# QHOSYN

QHOSYN (Pronounced "cosine") is the Quantum Harmonic Oscillator Synthesizer.

![](./QHOSYN_demo.gif)

QHOSYN is a synthesizer based on the fundamental physics model of a quantum harmonic oscillator. 

The software employs an accurate quantum simulation of an observable in a superposition of eigenstates up to the n=15 energy level to generate sound from  a quantum wavefunction based on the time-dependent Schrödinger equation.

It synthesizes sound from the simulation in 4 main ways. The first (Wavetable) treats the wavefunction as an evolving wavetable source for wavetable synthesis. The second (Probability Noise Band Synthesis) takes the probability curve of the wavefunction as a spectral curve to create colored noise. The third (IFFT) applies an inverse fourier transform to the wavefuction, treating it as a dynamic spectrum. The fourth (Panner) can be used to generate evolving stochastic sound and control signals from the probability distribution of the wavefunction. Through sonic exploration, the user can generate "quantum sounds" and gain an aesthetic understanding of quantum space.

QHOSYN was designed as part of my PhD research in order to facilitate the composition of my piece _Psi_. It is still in a rough state that was useful for me in my research, but may be confusing and buggy for anyone else to use. Version 1, which will be more user-friendly, will be released soon.

# Building

## Dependencies

terminal to run bash

git

cmake version 3.0 or higher

## How to setup
On a bash shell:

    git clone https://github.com/rodneydup/QHOSYN
    cd QHOSYN
    ./configure.sh
    ./run.sh

This will compile the project, and run the binary if compilation is successful.
