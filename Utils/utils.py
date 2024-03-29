"""Utilities for "dimuon_invm.py" and "Z_asymmetry.py" scripts."""

import logging
import ROOT

# pylint: disable=E1101
# (10/10)

def set_logger(name, log_level):
    """
    Create a new logger.

    :param name: name of the logger
    :type name: string, required
    :param log_level: level to set the logger
    :type log_level: e.g. "logging.DEBUG", required
    :return: logger

    """

    # Create logger and set its level
    logger = logging.getLogger(f"{name}")
    logging.basicConfig(level=log_level)
    handler = logging.FileHandler(f"{name}.log", "w+")
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - "
     "%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger

def write_fitresults(results, filename):
    """
    Function used to save results from RooFit fit in a txt file.

    :param results: results from a fit in RooFit
    :type results: RooFitResult, required
    :param filename: name of the txt output file
    :type filename: string, required

    """

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

# Dictionary with intialization parameters for the fit and plot styling.
# "particle" :
# (mean0,sigma0,signal,n° parameters for background,xmax(box),ymax(box),name)

FIT_INIT_PARAM = {
    "eta" : (ROOT.RooRealVar("mean", "mean", 0.552, 0.53, 0.56),
        ROOT.RooRealVar("sigma", "sigma", 0.07, 0.01, 0.1),
        "gaus",2, 0.39, 0.43, "#eta"),
    "rho" : (ROOT.RooRealVar("mean", "mean", 0.78, 0.75, 0.8),
        ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1),
        "gaus",2, 0.39, 0.8, "#rho"),
    "omega" : (ROOT.RooRealVar("mean", "mean", 0.78, 0.75, 0.8),
        ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1),
        "gaus",2, 0.39, 0.8, "#omega"),
    "phi": (ROOT.RooRealVar("mean", "mean", 1.02, 0.98, 1.06),
        ROOT.RooRealVar("sigma", "sigma", 0.004, 0.0001, 0.1),
        "gaus",2, 0.39, 0.8, "#phi"),
    "J-psi": (ROOT.RooRealVar("mean", "mean", 3.10, 2.9, 3.2),
        ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.),
        "Crystal ball", 3, 0.4, 0.8, "J/#psi"),
    "psi'" : (ROOT.RooRealVar("mean", "mean", 3.7, 3.6, 3.8),
        ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.),
        "gaus", 3, 0.9, 0.8, "#psi'"),
    "Z" : (ROOT.RooRealVar("mean", "mean", 91, 89, 93),
        ROOT.RooRealVar("sigma", "sigma",2, 0.01,4),
        "gaus", 2,0.42, 0.8, "Z"),
    "Y" : (0, 0, "3gaus", 3, 0.4, 0.8, "Y"),
    "Y1" : (ROOT.RooRealVar("mean1", "mean1", 9.4, 9.2, 9.7),
        ROOT.RooRealVar("sigma1", "sigma1", 0.005, 0.001, 1.)),
    "Y2" : (ROOT.RooRealVar("mean2", "mean2", 10, 9.8, 10.15),
        ROOT.RooRealVar("sigma2", "sigma2", 0.005, 0.001, 1.)),
    "Y3": (ROOT.RooRealVar("mean3", "mean3", 10.3, 10., 10.5),
        ROOT.RooRealVar("sigma3", "sigma3", 0.005, 0.001, 1.))
}
