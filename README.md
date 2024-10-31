# Adaptive-Particularly-Tunable-Fuzzy-PSO-(APT-FPSO)-algorithm
MATLAB code for the Adaptive Particularly Tunable Fuzzy PSO (APT-FPSO) algorithm, featuring adaptive fuzzy logic for dynamic parameter tuning to enhance optimization performance.

![Project Image](Flowchart.jpg)

## Overview
This repository features the Adaptive Particularly Tunable Fuzzy Particle Swarm Optimization (APT-FPSO) algorithm, an advanced variant of the standard Particle Swarm Optimization (PSO). Utilizing fuzzy logic, APT-FPSO adaptively tunes the learning coefficients for each particle at every iteration, enhancing the algorithm's balance between exploration and exploitation. The repository focuses on a single benchmark function, Griewangk, and demonstrates APT-FPSO’s effectiveness through detailed statistical analysis. This implementation is ideal for exploring the algorithm’s capabilities in complex optimization scenarios that require robust performance. For more detailed information, please read this [paper](https://ijfs.usb.ac.ir/article_5111.html).

## Features
- The measurement data for one single road test is provided.
- The road test corresponds to three hard-braking maneuvers on a road with a high friction coefficient.
- The friction coefficient of the road is measured from the maximum braking longitudinal acceleration before ABS activation.

## Usage
To run the main file for the project, use the following command in MATLAB:

```matlab
pso_simulation_fit_fun_Griewangk_fl_4_8_20.m
