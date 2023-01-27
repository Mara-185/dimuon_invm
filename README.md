# dimuon_invm

[![Documentation Status](https://readthedocs.org/projects/dimuon-invm/badge/?version=latest)](https://dimuon-invm.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://app.travis-ci.com/github/Mara-185/dimuon_invm.svg?branch=master)](https://app.travis-ci.com/github/Mara-185/dimuon_invm)

Il pacchetto propone due analisi differenti su dati di tipo NanoAOD:
* **dimuon_invm**: Seleziona le coppie di dimuoni con carica opposta al fine di riprodurre le risonanze nello spettro in massa dei dimuoni. Per ciascuna risonanza è possibile riprodurre plot di impulso trasverso, pseudorapidità e angolo azimutale e anche fare un fit del picco.
E.g. for the Y(1S,2S,3S):

![Image](Y_fit.png "icon")

* **Z_asymmetry**: Seleziona coppie di muoni che passano una particolare selezione (referenza: [link](https://arxiv.org/abs/1806.00863)) e riproduce il plot della Forward-Backward asymmetry in sei diversi range di rapidita'.
