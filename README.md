# POMAAR

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

POMAAR (Polarimetric MIMO Arrays for Automotive Radars) is focused on the development of
polarimetric MIMO arrays for advanced driver-assistance systems. This interdisciplinary
research involves solving fundamental challenges in automotive radar, including the development
of new theory and calibration procedures to enhance the classification of vulnerable road users.

## Project Organization

```
├── LICENSE            <- Open-source license
├── README.md          <- The top-level README for anyone using this project.
├── data
│   ├── calibration    <- Calibration data (matrices, antenna patterns).
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
│
├── docs               <- A default mkdocs project; see www.mkdocs.org for details
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         and a short description, for example `1.0-data-exploration`.
│
├── pyproject.toml     <- Project configuration file with package metadata for 
│                         pomaar and configuration for tools like black
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
└── src                <- Source code for use in this project.
    └── pomaar             <- POMAAR code.
        ├── __init__.py             <- Makes pomaar a Python module.
        ├── config.py               <- Store useful variables and configuration.
        ├── dsp.py                  <- DSP algorithms, such as beamforming and CFAR detectors.
        ├── io.py                   <- Scripts to parse raw ADC data for DSP.
        └── plots.py                <- Code to create visualizations.
```

--------

