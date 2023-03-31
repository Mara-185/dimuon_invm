"""Tests on \tools.cpp\" module (C++ shared library)."""

import logging
import unittest
import os
import sys
import ROOT

# pylint: disable=E1101
# (8.62/10)

# Add my modules to the path
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT_DIR)
sys.path.insert(0, os.path.abspath('../Utils'))
import utils

#Import shared library to test (Necessary for TRAVIS CI)
os.chdir("../Utils")
ROOT.gInterpreter.ProcessLine('#include "tools.h"')
ROOT.gInterpreter.ProcessLine('.L tools.cpp+')
os.chdir(ROOT_DIR)


# Create logger
logger = utils.set_logger("Unit test", logging.DEBUG)


def create_example():
    """Create a known example to test the functions."""

    # Create a standard vector with transvers momentum values
    pt = ROOT.std.vector("float")()
    pt.reserve(2)
    pt.emplace_back(57.2187)
    pt.emplace_back(31.1705)

    # Create a standard vector with pseudorapidity values
    eta = ROOT.std.vector("float")()
    eta.reserve(2)
    eta.emplace_back(2.21451)
    eta.emplace_back(1.15199)

    # Create a standard vector with azimuthal angle values
    phi = ROOT.std.vector("float")()
    phi.reserve(2)
    phi.emplace_back(3.13081)
    phi.emplace_back(-0.552126)

    # Create a standard vector with mass values
    mass = ROOT.std.vector("float")()
    mass.reserve(2)
    mass.emplace_back(0.1056580)
    mass.emplace_back(0.1056580)

    charge = ROOT.std.vector("int")()
    charge.reserve(2)
    charge.emplace_back(+1)
    charge.emplace_back(-1)

    return pt, eta, phi, mass, charge


class ZAsymmetryTest(unittest.TestCase):
    """Define a class with the tests."""

    # Test case
    pt, eta, phi, mass, charge = create_example()
    dilepton = ROOT.dilepton_vec(pt[0], eta[0], phi[0], mass[0], pt[1], \
        eta[1], phi[1], mass[1])
    cos_rap = ROOT.cos_rapidity(pt[0], eta[0], phi[0], mass[0], charge[0], pt[1], \
        eta[1], phi[1], mass[1])
    weights = ROOT.weights(dilepton[0], dilepton[3], cos_rap[0])

    # Actual tests

    def test_dilepton_vec(self):
        """Test on \"dilepton_vec\" function."""

        logger.info("Test on \"dilepton_vec\" function...")
        self.assertAlmostEqual(self.dilepton[0], 34.4752, 4)
        self.assertAlmostEqual(self.dilepton[1], 2.87066, 4)
        self.assertAlmostEqual(self.dilepton[2], -2.66773, 4)
        self.assertAlmostEqual(self.dilepton[3], 93.9916, 4)


    def test_cos_rapidity(self):
        """Test on \"cos_rapidity\" function."""

        logger.info("Test on \"cos_rapidity\" function...")
        self.assertAlmostEqual(self.cos_rap[0], -0.48296, 5)
        self.assertAlmostEqual(self.cos_rap[1], 1.82757, 5)


    def test_weights(self):
        """Test on \"weights\" function."""

        logger.info("Test on \"weights\" function...")
        self.assertAlmostEqual(self.weights[0], 0.05956, 5)
        self.assertAlmostEqual(self.weights[1], 0.15429, 5)


if __name__ == "__main__":

    # Create logger
    logger.info("Starting the tests...")

    # Start all tests
    unittest.main()
