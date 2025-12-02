# Combindexer

Universal processing pipeline for combinatorial indexing data with flexible barcode correction and UMI assignment.
Work in progress--in final stages of development. This version works; future changes are in the dev branch.

## Overview

Combindexer is a flexible barcoding pipeline for custom single-cell genomics methods, with a special focus on combinatorial indexing. This pipeline offers several key features not found in other, similar packages:

- **Flexible Barcode Configuration**: Define multiple barcode types with custom positions, correction distances, and orientations. Uniquely define cells by combinations of barcodes. 
- **Individual barcode whitelists**: Packages such as STARsolo allow multiple barcodes, but require a whitelist that individually contains each combination of barcodes. With combinatorial indexing methods, this can be millions of combinations! Combindexer allows individual barcode whitelists to be specified for each barcode.
- **Flexible barcode correction**: Barcode correction can account for SNPs, indels, 
- **Flexible export for compatibility with downstream gene counting programs:** Barcodes can be reconfigured in near-10X format for compatibility with STARSolo (see example), alevin-fry, etc. Alternatively, they can also be exported as a CSV. The order in which barcodes are exported is also customizable.

## Installation

To install Combindexer using `pyproject.toml`, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd combindexer
   ```

2. **(Recommended) Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate
   ```

3. **Install using `pip` (supports PEP 517/518 builds):**
   ```bash
   pip install .
   ```

   Alternatively, if you want to install in editable (development) mode:
   ```bash
   pip install -e .
   ```

4. **Verify installation:**
   ```bash
   python -m pip show combindexer
   ```

## Instructions
To follow (see example in the meantime)
