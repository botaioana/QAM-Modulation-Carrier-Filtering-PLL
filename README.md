# QAM System Design and Simulation

## Project Overview
This project aims to simulate a full digital communication system using Quadrature Amplitude Modulation (QAM). It covers the generation of QAM constellations with Gray coding for minimizing bit errors, as well as key signal processing steps such as carrier extraction, filtering, limiting, and synchronization via a Phase-Locked Loop (PLL). This work provides the foundation for further exploration of higher-order QAM constellations, advanced filtering techniques, and improved PLL models.

## Current Features
- **QAM Constellation Generation**: The current implementation generates a square-shaped QAM constellation (16-QAM) and maps Gray-coded symbols to each point. The system can be extended to other orders of QAM, such as 64-QAM or 256-QAM.
- **Carrier Extraction**: Simulates the extraction of the carrier signal by raising the modulated signal to the fourth power, useful for recovering the carrier from a noisy modulated signal.
- **Filtering and Limiting**: Although not fully implemented, placeholders are present to simulate filtering and signal limiting for noise reduction and signal quality enhancement.
- **PLL Synchronization**: A Phase-Locked Loop (PLL) is applied to the extracted carrier signal to synchronize the system, ensuring proper recovery of the transmitted data.

## Future Work and Planned Features
- **Higher-Order QAM Support**: Extend the QAM constellation generator to support higher-order constellations such as 64-QAM, 128-QAM, and 256-QAM for more complex modulations.
- **Advanced Filtering Techniques**: Implement advanced filtering algorithms to improve signal integrity after carrier extraction, such as low-pass filters or adaptive filters.
- **Complete PLL Circuit Simulation**: Develop a more detailed and robust simulation of the PLL to improve carrier synchronization and reduce phase jitter.
- **Error Rate Analysis**: Add functionality to simulate and analyze the bit error rate (BER) under various noise and interference conditions.
- **Demodulation and Data Recovery**: Implement demodulation techniques for recovering transmitted data from the QAM signal, providing end-to-end simulation of a communication system.

## Requirements
- MATLAB (any recent version)

## Code Structure
- **`generate_gray.m`**: Generates Gray codes for the QAM constellation axes, ensuring minimal bit error between adjacent constellation points.
- **`main.m`**: Main script responsible for generating the QAM constellation, simulating carrier extraction, filtering, and PLL synchronization.
- **Future Scripts**: Planned future additions will include error rate analysis, advanced signal processing, and demodulation functions.

## How to Run
1. Clone the repository.
2. Run the `main.m` script to visualize the QAM constellation and simulate carrier extraction and PLL synchronization.
3. As the project expands, additional scripts and examples will be included for running more advanced simulations.

## Example Output
The figure below shows an example 16-QAM constellation, where each point is labeled with its corresponding Gray-coded bits.

![QAM Constellation](https://didatec-my.sharepoint.com/:i:/r/personal/bota_io_ioana_student_utcluj_ro/Documents/Materiale%20Facultate/an4%20-%20sem%207/TD/Assignment%201/Screenshot%202024-10-07%20152433.png?csf=1&web=1&e=fqsipi)
