# Digital-Comm---Line-Codes
A digital communication project that aims to insinuate various line codes and their properties from an ensemble.

# Overview
This project explores the implementation of a software-defined radio system. The fundamental characteristic of a software radio is that the software defines the transmitted waveforms and demodulates the received waveforms. This approach contrasts with traditional radios, which rely on analog circuitry or a combination of analog and digital processing. The project focuses on generating, analyzing, and processing line codes used in digital communications.

# Project Objectives
### 1. Transmission and Line Coding:
- Implement transmission using different line coding techniques, including Unipolar Signaling, Polar Non-Return to Zero (NRZ), and Return to Zero (RZ).
- Convert binary data into a corresponding voltage level using MATLAB code and simulate the transmission over time.

### 2. Random Process Analysis:
- Generate an ensemble of 500 waveforms, each containing 100 bits, using the specified line codes.
- Compute and analyze the statistical mean, time mean, and autocorrelation functions.
- Investigate the stationarity and ergodicity of the random process.
- Determine the bandwidth of the transmitted signal.

# Implementation Details
The implementation involves the following steps:

### 1. Data Preparation:

- Binary data is prepared and converted into voltage levels using MATLAB.
### 2. Transmission Simulation:

- The code maps binary '1' to a positive voltage level (+A) and binary '0' to a negative voltage level (-A) in Polar NRZ signaling.
- The generated waveforms simulate the process of Digital-to-Analog Conversion (DAC) in real-time transmission.

### 3. Ensemble Generation:

- 500 waveforms are generated for each line coding technique.
- Each waveform starts from a random initial time shift to simulate real-world randomness.

### 4.Statistical Analysis:

- The project calculates the statistical mean, autocorrelation functions, and assesses whether the process is stationary or ergodic.
- Bandwidth analysis is conducted to understand the spectral characteristics of the transmitted signals.

# Project Structure
- Scripts/: MATLAB scripts for generating waveforms, analyzing statistics, and plotting results.
- Results/: Output files including plots of waveforms, statistical analyses, and bandwidth calculations.

# Acknowledgments
- Cairo University, Faculty of Engineering, Electronics and Communications Dept.
