"""Utilities for plotting scripts."""

import logging
import ROOT

def set_logger(name):

    # Create logger and set its level
    logger = logging.getLogger(f"{name}")
    logging.basicConfig(level=logging.DEBUG)
    handler = logging.FileHandler(f"{name}.log", "w+")
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - "
     "%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger

def write_fitresults(results, filename):
    """Function used to save results from RooFit fit in a txt file"""

    text = ROOT.std.ofstream(filename)
    results.printMultiline(text, 1111, True)
    text.close()


# Dictionary with particles' mass range for the fit.
# "particle" : (lower_mass_limit, upper_mass_limit)

PARTICLES_MASS_RANGE = {
  "eta" : (0.52, 0.57),
  "rho" : (0.72, 0.84),
  "omega": (0.72, 0.84),
  "phi" : (0.96, 1.07),
  "J-psi" : (2.65, 3.55),
  "psi'" : (3.55, 3.85),
  "Z" : (82, 98),
  "Y": (8.5, 11),
  "Y1": (9.1, 9.75),
  "Y2": (9.75, 10.2),
  "Y3": (10.1, 10.6)
}

# Dictionary with intialization parameters for the fit and plot styling
# "particle" :
# (mean0,sigma0,signal,n° parameters for background,xmax(box),ymax(box),name,
#   n° total parameters )
FIT_INIT_PARAM = {
    "eta" : (ROOT.RooRealVar("mean", "mean", 0.552, 0.53, 0.56),
        ROOT.RooRealVar("sigma", "sigma", 0.07, 0.01, 0.1),
        "gaus",2, 0.39, 0.43, "#eta", 6),
    "rho" : (ROOT.RooRealVar("mean", "mean", 0.78, 0.75, 0.8),
        ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1),
        "gaus",2, 0.39, 0.8, "#rho", 6),
    "omega" : (ROOT.RooRealVar("mean", "mean", 0.78, 0.75, 0.8),
        ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1),
        "gaus",2, 0.39, 0.8, "#omega", 6),
    "phi": (ROOT.RooRealVar("mean", "mean", 1.02, 0.98, 1.06),
        ROOT.RooRealVar("sigma", "sigma", 0.004, 0.0001, 0.1),
        "gaus",2, 0.39, 0.8, "#phi", 6),
    "J-psi": (ROOT.RooRealVar("mean", "mean", 3.10, 2.9, 3.2),
        ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.),
        "Crystal ball", 3, 0.4, 0.8, "J/#psi", 9),
    "psi'" : (ROOT.RooRealVar("mean", "mean", 3.7, 3.6, 3.8),
        ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.),
        "gaus", 3, 0.9, 0.8, "#psi'", 7),
    "Z" : (ROOT.RooRealVar("mean", "mean", 91, 89, 93),
        ROOT.RooRealVar("sigma", "sigma",2, 0.01,4),
        "gaus", 2,0.42, 0.8, "Z", 6),
    "Y" : (0, 0, "3gaus", 3, 0.4, 0.8, "Y", 13),
    "Y1" : (ROOT.RooRealVar("mean1", "mean1", 9.4, 9.2, 9.7),
        ROOT.RooRealVar("sigma1", "sigma1", 0.005, 0.001, 1.)),
    "Y2" : (ROOT.RooRealVar("mean2", "mean2", 10, 9.8, 10.15),
        ROOT.RooRealVar("sigma2", "sigma2", 0.005, 0.001, 1.)),
    "Y3": (ROOT.RooRealVar("mean3", "mean3", 10.3, 10., 10.5),
        ROOT.RooRealVar("sigma3", "sigma3", 0.005, 0.001, 1.))
}


# Z analysis
# Dictionaries with mass bin values and eta bin values

# Mass range from article : 60 < M < 120:
    # 60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120;
MASS_BIN = {
    "0": (60, 70),
    "1": (70, 78),
    "2": (78, 84),
    "3": (84, 87),
    "4": (87, 89),
    "5": (89, 91),
    "6": (91, 93),
    "7": (93, 95),
    "8": (95, 98),
    "9": (98, 104),
    "10": (104, 112),
    "11": (112, 120)
}

# Eta bins of equal size for |yll| < 2.4:
    # 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4
RAPIDITY_BIN = {
    "0": (0.0, 0.4),
    "1": (0.4, 0.8),
    "2": (0.8, 1.2),
    "3": (1.2, 1.6),
    "4": (1.6, 2.0),
    "5": (2.0, 2.4),
}
