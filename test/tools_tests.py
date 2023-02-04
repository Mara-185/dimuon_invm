"""Tests on \tools.cpp\" module (C++ shared library)."""

import logging
import unittest
import os
import sys
import ROOT

# Add my modules to the path
root_utils = os.path.abspath('../Utils')
sys.path.insert(0, root_utils)
#sys.path.insert(0, os.path.join(root_utils, "utils"))
import utils

#Import shared library to test
#ROOT.gInterpreter.ProcessLine('#include "tools.h"')
ROOT.gSystem.Load('../Utils/tools_cpp.so')


def create_example():
    """Create a known example to test the functions."""

    # Create a standard vector with transvers momentum values
    pt =  ROOT.std.vector("float")()
    pt.emplace_back(57.2187)
    pt.emplace_back(31.1705)

    # Create a standard vector with pseudorapidity values
    eta = ROOT.std.vector("float")()
    eta.emplace_back(2.21451)
    eta.emplace_back(1.15199)

    # Create a standard vector with azimuthal angle values
    phi = ROOT.std.vector("float")()
    phi.emplace_back(3.13081)
    phi.emplace_back(-0.552126)

    # Create a standard vector with mass values
    mass = ROOT.std.vector("float")()
    mass.emplace_back(0.1056580)
    mass.emplace_back(0.1056580)

    return pt, eta, phi, mass


class Z_asymmetryTest(unittest.TestCase):
    """Define a class with the tests."""

    def test_dilepton_vec(self):
        """Test on \"dilepton_vec\" function."""

        logger.info("Test on \"dilepton_vec\" function...")
        pt, eta, phi, mass = create_example()
        dilepton = ROOT.dilepton_vec(pt[0], eta[0], phi[0], mass[0], pt[1], \
            eta[1], phi[1], mass[1])
        self.assertAlmostEqual(dilepton[0], 34.4752, 4)
        self.assertAlmostEqual(dilepton[1], 2.87066, 4)
        self.assertAlmostEqual(dilepton[2], -2.66773, 4)
        self.assertAlmostEqual(dilepton[3], 93.9916, 4)


    def test_cos_rapidity(self):
        """Test on \"cos_rapidity\" function."""

        logger.info("Test on \"cos_rapidity\" function...")
        pt, eta, phi, mass = create_example()
        cos_rap = ROOT.cos_rapidity(pt[0], eta[0], phi[0], mass[0], pt[1], \
            eta[1], phi[1], mass[1])
        self.assertAlmostEqual(cos_rap[0], 0.48296, 5)
        self.assertAlmostEqual(cos_rap[1], 1.82757, 5)


    def test_weights(self):
        """Test on \"weights\" function."""

        logger.info("Test on \"weights\" function...")
        pt, eta, phi, mass = create_example()
        dilepton = ROOT.dilepton_vec(pt[0], eta[0], phi[0], mass[0], pt[1], \
            eta[1], phi[1], mass[1])
        cos_rap = ROOT.cos_rapidity(pt[0], eta[0], phi[0], mass[0], pt[1], \
            eta[1], phi[1], mass[1])
        weights = ROOT.weights(dilepton[0], dilepton[3], cos_rap[0])
        self.assertAlmostEqual(weights[0], 0.05956, 5)
        self.assertAlmostEqual(weights[1], 0.15429, 5)


if __name__ == "__main__":

    # Create logger
    logger = utils.set_logger("Unit test", logging.DEBUG)
    logger.info("Starting the tests...")

    # Start all tests
    unittest.main()
