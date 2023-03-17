# dimuon_invm

[![Documentation Status](https://readthedocs.org/projects/dimuon-invm/badge/?version=latest)](https://dimuon-invm.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://app.travis-ci.com/github/Mara-185/dimuon_invm.svg?branch=master)](https://app.travis-ci.com/github/Mara-185/dimuon_invm)

Il pacchetto propone due analisi differenti su dati di tipo NanoAOD:
* **dimuon_invm**: Seleziona le coppie di dimuoni con carica opposta al fine di riprodurre le risonanze nello spettro in massa dei dimuoni. Per ciascuna risonanza è possibile riprodurre plot di impulso trasverso, pseudorapidità e angolo azimutale e fare un fit del picco per ricavare il valore della massae della sezione d'urto. Risoluzione?
(Usa nanoAOD_DoubleMu_Parked)

E.g. for the Y(1S,2S,3S):

![Image](./dimuon_invm/Fit/Y_fit.png "icon")

* **Z_asymmetry**: Seleziona coppie di muoni che passano una particolare selezione (referenza: [link](https://arxiv.org/abs/1806.00863)) e riproduce il plot della Forward-Backward asymmetry in sei diversi range di rapidita'.

E.g. for data of Run2012B-C_SingleMu:

![Image](./Z_asymmetry/Plot/afb_y_all.png "icon")
