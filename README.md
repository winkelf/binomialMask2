# binomialMask2

This repository contains tools for generating a binomialMask over images taken with Skipper CCDs used in MINOS measurement for the SENSEI experiment, designed for integration into the official SENSEI framework.

---

## Overview

The project consists of three main components:

### 1. `configMaker.py`

A Python script that generates a JSON configuration file containing all image and mask files to be used as inputs.

**Purpose:**
- Collect image file paths
- Collect corresponding mask file paths
- Generate a structured JSON configuration file
- Provide consistent input formatting for the analysis program

---

### 2. `helperFunctions.C`

A C++ source file containing auxiliary utility functions used by the main analysis program.

**Purpose:**
- Provide reusable helper routines
- Keep the main program modular and organized
- Support shared logic across the project

---

### 3. `binomialMask.C`

The main program of the repository.

**Key features:**
- Contains the primary function: `analyze`
- Implements the binomial mask analysis logic
- Designed to be incorporated into the SENSEI framework

This file represents the core computational component of the project.

---

## Intended Integration

The `analyze` function defined in `binomialMask.C` is intended to be integrated into the SENSEI framework so it can operate within the larger SENSEI analysis pipeline.

---
