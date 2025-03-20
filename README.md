# Shock Wave and Thermals Simulation in a Virtual Wind Tunnel with Maching Learning (CV)

This repository contains my scientific work related to the modeling of **shock waves** and **thermals** (mushroom-shaped structures) in a virtual aerodynamic wind tunnel. The project involves the implementation of a numerical method to solve a two-dimensional problem, simulating the behavior of matter inside the wind tunnel.

## Project Overview

### Goal
The primary goal of this project is to **reconstruct the initial parameters** of a real-world experiment conducted in the wind tunnel at the Physics Department of Moscow State University. Experimental snapshots of wave propagation and thermals development are available, but the initial conditions of the experiment are unknown. By generating similar images (graphs) using a numerical method, we aim to find the closest matches to the experimental images. This will allow us to approximate the initial parameters of the real experiment with a certain degree of accuracy.

### Methodology
1. **Numerical Algorithm**:
   - A numerical algorithm was developed in **C** (primarily written by my scientific supervisor).
   - The algorithm outputs files of type `fld.dat`, which describe the main physical quantities at each grid node.

2. **Python Visualization**:
   - I wrote a **Python** script to process these `fld.dat` files and generate wave distribution graphs (e.g., gradient density magnitude) that closely resemble the experimental images.
   - To create a diverse dataset, a large number of such images were generated.
   - An even number of pictures were also generated - a pair: an ordinary picture and a picture with a certain Gaussian noise superimposed on it. I will describe the motivation for adding noise below.

3. **Integration**:
   - The Python visualization code was integrated into the C code using a **subprocess** that calls the Python script and deletes the current `fld.dat` file after processing.
   - This sequential approach avoids memory overflow, as generating and storing 100,000+ `fld.dat` files would require over **400 GB** of storage, which exceeds the capacity of my computer.

4. **Dataset Augmentation**:
   - The dataset was augmented by stretching **40%** of the images along the X-axis to increase variability.
   - Augmentation is needed to increase the generalizing ability of the neural network, as well as so that the neural network can work with images of non-fixed size.
   - The result is about 80,000 pairs of images.
5. **A Neural Network for Noise Removal**:
During the work, it was observed that experimental images from the real wind tunnel contain characteristic noise, while numerically generated images are noise-free. To address this, a neural network was developed to denoise experimental images, improving subsequent image comparison metrics. The implementation started with a modified U-Net architecture adapted for variable-size inputs.
   - Base Architecture: Built upon the classic U-Net, which is widely used for image-to-image tasks.
   - 
---

## Repository Structure
- **C Code**: Contains the numerical algorithm for simulating shock waves and thermals.
- **Python Scripts**: Includes scripts for processing `fld.dat` files and generating visualizations.
- **Dataset**: A collection of generated images and augmented data.
- **Results**: Experimental and simulated images for comparison.

---

## Key Features
- **Efficient Memory Management**: Sequential processing of `fld.dat` files to avoid memory overflow.
- **Dataset Augmentation**: Enhanced dataset diversity through image stretching.
- **Integration of C and Python**: Seamless workflow between numerical simulation and visualization.

---
