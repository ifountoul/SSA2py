<p align="center">
  <a href="">
    <img src="logo.jpg" width="400" alt="SSA2py logo">
  </a>
</p>

<h3 align="center">SSA2py</h3>

<p align="center">
   Unleashing Source Scanning at the Speed of Python
  <br>
  <a href="https://ssa2py.readthedocs.io/en/latest/"><strong>Explore SSA2py docs »</strong></a>
  <br>
  <br>
  <a href="https://github.com/ifountoul/SSA2py/issues">Report bug</a>
</p>


# SSA2py 1.0.0

Welcome to the public repo for SSA2py.

## What is SSA2py?

[SSA2py]() **is an open-source python project that follows the Source-Scanning Algorithm (SSA)**.
It provides interconnection with FDSN Compliant Web Services and it is adapted to run in GPU and CPU multiprocessing architectures. 
The aim of SSA2py is to provide rapid and accurate calculations of SSA method in near-realtime conditions.

The official documentation is hosted on [Read the Docs](https://ssa2py.readthedocs.io/en/latest/).

## Status

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/python-3.10-blue.svg)
![Repo Size](https://img.shields.io/github/repo-size/Sulstice/global-chem)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![Maturity badge - level 2](https://img.shields.io/badge/Maturity-Level%202%20--%20First%20Release-yellowgreen.svg)](https://github.com/tophat/getting-started/blob/master/scorecard.md)

## Table of contents

- [Quick start](#quick-start)
- [Supported Platforms](#supported-platforms)
- [Community](#community)
- [Contribution](#contribution)
- [Creators](#creators)
- [Thanks](#thanks)
- [Copyright and license](#copyright-and-license)


## Quick Start
- `git clone https://github.com/ifountoul/SSA2py.git`
- `cd SSA2py`
- `conda env create -f environment.yml`
- `conda activate SSA2PY`
-  Create the events directory (Events Dir) and the traveltimes directory (Traveltimes/Save) declared in the configuration file.
- `python3 SSA2py.py --download`
- `python3 SSA2py.py` **(Everything OK? Ready to go!)**

### Prerequisites
- Install conda
- If you install SSA2py on a brand new system install the C and C++ compilers before installing Anaconda.
- Make sure that you have conda-forge in your channels (`conda config --show channels`). You can add it by executing `conda config --add channels conda-forge`.
- To use GPU install the cudatoolkit throught anaconda. Please check the CUDA and NVIDIA driver versions.

To learn more about using SSA2py, follow our [guide here](https://ssa2py.readthedocs.io/en/latest/) or see the [examples](https://ssa2py.readthedocs.io/en/latest/applications.html).

## Supported Platforms

Our software is tested on the following platforms:
- Linux (Ubuntu 22.04 LTS)
- Windows (Windows 11 Pro - WSL)

**Note:** SSA2py is not vigorously tested on macOS at the moment. Contributions and feedback related to macOS testing are welcome.

## Community

Have questions, comments or feedback? Start a [discussion](https://github.com/ifountoul/SSA2py/discussions).

## Contribution

Found a bug? Please submit an [issue](https://github.com/ifountoul/SSA2py/issues).

## Creators

**Ioannis Fountoulakis**

:email: ifountoul@noa.gr

**Christos Evangelidis**

:email: cevan@noa.gr


## Thanks 

<a href="https://www.elidek.gr/en/homepage/">
  <img src="https://www.elidek.gr/wp-content/themes/elidek/images/elidek_logo_en.png" alt="H.F.R.I" width="310" height="90">
</a>

The research was supported by the Hellenic Foundation for Research and Innovation ([H.F.R.I.](https://www.elidek.gr/en/homepage/)) under the 
“First Call for H.F.R.I. Research Projects to support Faculty members and Researchers and the procurement of high-cost research equipment grant” (SIREN, Project Number: 910).

<a href="https://www.noa.gr/en/">
  <img src="https://www.noa.gr/wp-content/uploads/2019/12/noa_logo.svg" alt="NOA" width="110" height="110">
</a>

Thanks to [NOA](https://www.noa.gr/en/) for providing the infrastructure to develop this program!

## Copyright and license

Code released under the GNU [GENERAL PUBLIC LICENSE Version 3](https://github.com/ifountoul/SSA2py-Ghost/blob/master/LICENSE)
