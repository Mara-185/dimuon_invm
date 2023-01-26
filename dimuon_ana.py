import ROOT
import argparse
import logging
import os
import math
import sys
import glob
import numpy as np
import utils
from pathlib import Path

"""The script takes as argument:
    - the data file (URL) of dileptons (-f), for example: "root:
     //eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root";
    - the string of the particle's name (-p) to analyze and fit, among those in
      the dimuon spectrum, which are :
      "eta", "rho","omega", "phi", "J-psi", "psi'", "Y", "Z".
      Better put the string in quotes, because for example for the "psi'" the
      " ' " character gives some troubles.

      There are different functions for the main analysis, among which:
      - "leptons_analysis" which selects the couple of muons and electrons which
        are interesting in order to create the dimuon mass spectrum and for
        further analysis;
      - "mumu_spectrum" which plot the dimuon mass spectrum.
      - "resonance_fit" solves every resonance and returns the plot with the fit
        and a txt with the fit results.
      - "resonance_prop" creates different plots of the main characteristics of
        the particle chosen in the spectrum;

      For the extimation of the weak mixing angle studying the angular properties
      of the Z boson:
      - "weight" which calculates other useful variables for the analysis;
      - "afb" which estimate, from the varibles obtained by the previous function,
        the mean values of Afb (forward-backward asymmetry) in different bins
        of mass and pseudorapidity (in total 6 bins of pseudorapidity and 12 bins
        of mass).

    In the analysis the following version have been used:
    - Python v3.8
    - ROOT v6.24 ("source ~/root/bin/thisroot.sh" command needed before starting
        the analysis to set the environment of ROOT)
"""

# def main(infile, particle="all"):

# Creating the logger and setting its level
logger = utils.set_logger("Analysis")

def leptons_analysis(infile):
    """ It takes in input a nano-AOD data file and returns two root files, named:
    "dimuon.root" and "dielectron.root" both with four columns, each one related
    to different quantity, for example:
    ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi"].
    It returns the data cached in memory, in order to speed up further analysis.
    """

    # Enable parallel analysis
    nthreads = 8
    ROOT.ROOT.EnableImplicitMT(nthreads)

    logging.info("Start the analysis...")

    # Create an RDataFrame of useful data from the root file.
    root_file = ROOT.TFile.Open(infile, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug("The RDataFrame have been created.")

    # Filter on two muons and two electrons with opposite charge
    rdf_mu = rdf.Filter("nMuon==2","Selection of two muons").\
     Filter("Muon_charge[0]!=Muon_charge[1]",\
                "Selection of muons with opposite charge")
    logger.info("The cut on two muons with opposite charge is done.")

    rdf_e = rdf.Filter("nElectron==2","Selection of two electrons").\
     Filter("Electron_charge[0]!=Electron_charge[1]",\
                "Selection of electrons with opposite charge")
    logger.info("The cut on two electrons with opposite charge is done.")

    # Print cutflows
    logger.info("CutFlow muons:")
    rdf_mu.Report().Print()

    logger.info("CutFlow electrons:")
    rdf_e.Report().Print()

    # Save Node Graph
    ROOT.RDF.SaveGraph(rdf, "rdf(first).dot")
    del rdf


    # By "Jitting", two functions are created:
    #   - "dilepton_vec" takes in input pt (transverse momentum), eta
    #   (pseudorapidity), phi (azimuthal angle), and mass of the chosen two
    #   particles and it returns a four-vector with the dilepton invariant mass,
    #   pt, phi and eta.
    #   -"cos_rapidity" takes in input pt, eta, phi and mass of the chosen two
    #   particles and it returns a vector with the values of rapidity and cosine
    #   of the angle of the negative leptons in the Collins-Soper frame of the
    #   dilepton system.
    ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<double> dilepton_vec(float pt0, float eta0, float phi0,
        float mass0, float pt1, float eta1, float phi1, float mass1)
    {
    	ROOT::Math::PtEtaPhiMVector p0, p1;
    	p0.SetCoordinates(pt0, eta0, phi0, mass0);
    	p1.SetCoordinates(pt1, eta1, phi1, mass1);
    	ROOT::VecOps::RVec<double> P{(p1+p0).Pt(), (p1+p0).Eta(), (p1+p0).Phi(),
                (p1+p0).M()};
    	return P;
    }


    ROOT::VecOps::RVec<float> cos_rapidity(float pt0, float eta0, float phi0,
            float mass0, float pt1, float eta1, float phi1, float mass1)
    {
        float pz0, pz1, E0, E1, P0_1, P0_2, P1_1, P1_2, cos, numer, denom, mll,
            ptt, pzll, y;
        ROOT::Math::PtEtaPhiMVector p0, p1;
  		p0.SetCoordinates(pt0, eta0, phi0, mass0);
  		p1.SetCoordinates(pt1, eta1, phi1, mass1);

        pz0 = pt0*sinh(eta0);
        E0 = sqrt(pow(pz0,2)+pow(pt0,2)+pow(mass0,2));
        pz1 = pt1*sinh(eta1);
        E1 = sqrt(pow(pz1,2)+pow(pt1,2)+pow(mass1,2));
        P0_1 = (E0+pz0)/sqrt(2);
        P0_2 = (E0-pz0)/sqrt(2);
        P1_1 = (E1+pz1)/sqrt(2);
        P1_2 = (E1-pz1)/sqrt(2);
        numer = 2*((P0_1*P1_2) - (P0_2*P1_1));
        mll = pow((p1+p0).M(),2);
        ptt = pow((p1+p0).Pt(),2);
        denom = sqrt(mll*(mll+ptt));
        pzll = pz0+pz1;

        cos = (numer/denom)*(pzll/abs(pzll));
        y = (0.5)*log(((E0+E1)+(pz0+pz1))/((E0+E1)-(pz0+pz1)));

        ROOT::VecOps::RVec<float> A{cos,y};
  		return A;
  	}
    """)
    #cos =( ( 2*((P0_1*P1_2) - (P0_2*P1_1)) ) * (pzll) ) / ( (sqrt(mll*(mll+ptt))) * (abs(pzll)) );

    logger.debug("The two c++ functions have been defined.")

    # Define new columns in the RDataFrame with the variables needed for the analysis
    rdf_dimu = rdf_mu.\
        Define("Dimuon_mass", \
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[3]").\
        Define("Dimuon_pt",\
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_eta", \
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]").\
        Define("Dimuon_phi",\
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[2]").\
        Define("Dimuon_cos",\
            "cos_rapidity(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_y",\
            "cos_rapidity(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]")

    # Save Node Graph
    ROOT.RDF.SaveGraph(rdf_mu, "rdf_mu.dot")
    del rdf_mu

    rdf_diel = rdf_e.\
        Define("Dielectron_mass", "dilepton_vec(Electron_pt[0], Electron_eta[0],\
         Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
         Electron_phi[1], Electron_mass[1])[3]").\
        Define("Dielectron_pt", "dilepton_vec(Electron_pt[0], Electron_eta[0],\
         Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
         Electron_phi[1], Electron_mass[1])[0]").\
        Define("Dielectron_eta", "dilepton_vec(Electron_pt[0], Electron_eta[0],\
         Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
         Electron_phi[1], Electron_mass[1])[1]").\
        Define("Dielectron_phi", "dilepton_vec(Electron_pt[0], Electron_eta[0],\
         Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
         Electron_phi[1], Electron_mass[1])[2]").\
        Define("Dielectron_cos", "cos_rapidity(Electron_pt[0], Electron_eta[0],\
         Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
         Electron_phi[1], Electron_mass[1])[0]").\
        Define("Dielectron_y", "cos_rapidity(Electron_pt[0], Electron_eta[0], \
        Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], \
        Electron_phi[1], Electron_mass[1])[1]")
    logger.debug("The new columns of dimuons and dieletrons are defined.")

    # Save Node Graph
    ROOT.RDF.SaveGraph(rdf_e, "rdf_e.dot")
    del rdf_e

    # Lists with the names of the new columns are created in order to make a
    # snapshot of them
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi", \
        "Dimuon_cos", "Dimuon_y"]
    branchlist_e = [ "Dielectron_mass", "Dielectron_pt", "Dielectron_eta", \
        "Dielectron_phi", "Dielectron_cos", "Dielectron_y"]

    # Two snapshots are done to collect the useful physiscal quantity of the
    # dimuons and dielectrons in a single root file.
    # Data are also cached to speed up the analysis.
    logger.debug("Starting snapshots...")
    dimu_cached = rdf_dimu.Snapshot("dimuon", "dimuon.root", branchlist_mu).\
        Cache()
    diel_cached = rdf_diel.Snapshot("dielectron", "dielectron.root",
        branchlist_e).Cache()
    logger.info("The snapshots are done and data are cached.")

    # Close the file and return data cached of dimuons and dielectrons.
    root_file.Close()
    return dimu_cached, diel_cached


def mumu_spectrum(infile, mu_cached=None):
    """It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" named "dimuon.root" and plot the histogram of
    the dimuons invariant mass."""

    #Histogram of dimuon mass
    if mu_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        h = rdf.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")
    else:
        h = mu_cached.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    c = ROOT.TCanvas("dimuon spectrum", "#mu^{+}#mu^{-} invariant mass")
    c.SetLogx()
    c.SetLogy()
    h.GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
    h.GetXaxis().SetTitleSize(0.04)
    h.GetXaxis().CenterTitle()
    h.GetYaxis().SetTitle("Events")
    h.GetYaxis().SetTitleSize(0.04)
    h.GetYaxis().CenterTitle()
    h.Draw()

    # Labels
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.DrawLatex(0.165, 0.720, "#eta")
    label.DrawLatex(0.190, 0.772, "#rho,#omega")
    label.DrawLatex(0.245, 0.775, "#phi")
    label.DrawLatex(0.400, 0.850, "J/#psi")
    label.DrawLatex(0.410, 0.700, "#psi'")
    label.DrawLatex(0.485, 0.700, "Y(1, 2, 3S)")
    label.DrawLatex(0.795, 0.680, "Z")

    # Save results in pdf and png in the folder "Spectra"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_spectrum.pdf")
    c.SaveAs("dimuon_spectrum.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The files \"dimuon_spectrum.***\"(pdf, png) have been created.")


def mumu_spectrum_bump(infile, mu_cached=None):
    """It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" named "dimuon.root" and plot the histogram of
    the dimuons invariant mass with a cut for events which have pt>20 GeV."""

    if mu_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")
        del rdf
    else:
        rdf_m = mu_c.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    c = ROOT.TCanvas("dimuon spectrum", "#mu^{+}#mu^{-} invariant mass")
    c.SetLogx()
    c.SetLogy()
    h.SetTitle("Dimuon mass without bump")
    h.GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
    h.GetXaxis().SetTitleSize(0.04)
    h.GetXaxis().CenterTitle()
    h.GetYaxis().SetTitle("Events")
    h.GetYaxis().SetTitleSize(0.04)
    h.GetYaxis().CenterTitle()
    h.Draw()

    # Labels
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.DrawLatex(0.165, 0.720, "#eta")
    label.DrawLatex(0.190, 0.772, "#rho,#omega")
    label.DrawLatex(0.245, 0.775, "#phi")
    label.DrawLatex(0.400, 0.850, "J/#psi")
    label.DrawLatex(0.410, 0.700, "#psi'")
    label.DrawLatex(0.485, 0.700, "Y(1, 2, 3S)")
    label.DrawLatex(0.795, 0.680, "Z")
    label.DrawLatex(0.7, 0.85, "p_{t}> 20 GeV")

    # Save results in pdf and png in the folder "Spectra"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_spectrum_bump.pdf")
    c.SaveAs("dimuon_spectrum_bump.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The files \"dimuon_spectrum_bump.***\"(pdf, png) have been created.")


def mumu_eta(infile, mu_cached=None):
    """It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" named "dimuon.root" and plot the histogram of
    the dimuons pseudorapidity."""

    if mu_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon eta", "Dimuon eta",
            150, -4, 4), "Dimuon_eta")
        del rdf
        #del rdf_m????
    else:
        rdf_m = mu_cached.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon eta", "Dimuon eta",
            150, -4, 4), "Dimuon_eta")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    c = ROOT.TCanvas("dimuon spectrum", "#mu^{+}#mu^{-} pt")
    c.SetGrid()
    h.GetXaxis().SetTitle("#eta_{#mu^{+}#mu^{-}} [GeV]")
    h.GetXaxis().SetTitleSize(0.04)
    h.GetXaxis().CenterTitle()
    h.GetYaxis().SetTitle("Events")
    h.GetYaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetYaxis().CenterTitle()
    h.Draw()

    # Save results in pdf and png in the folder "Spectra"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_eta.pdf")
    c.SaveAs("dimuon_eta.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The files \"dimuon_eta.***\"(pdf,png) have been created.")


def resonance_fit(infile, particle, mu_cached=None):
    """It takes in input the root data file obtained with "leptons_analysis"
        and a string with the name of the resonance to fit. Possible arguments
        are: \"eta\", \"rho\",\"omega\", \"phi\", \"J-psi\", \"psi'\", \"Y\",
        \"Z\", \"all\."
        It retrieves a plot with the fit and a txt file with fit results, for
        each resonance.

         and a workspace
        ????????????????????????????????????????????????????????"""

    # An auxiliary TTree from the root data file is created.
    if mu_cached == None:
        f = ROOT.TFile.Open(infile)
        tree_name = infile.replace(".root", "")
        tree = f.Get(f"{tree_name}")
        logger.debug("The tree is created.")
    else:
        logger.debug("Devi aggiornare la versione di root a 6.26 e usare"
            "RooDataSetHelper")

    low_edge_pt = 20
    upper_edge_pt = 100
    low_edge_eta = -1.2
    upper_edge_eta = 1.2

    error_flag=1
    bkg = "Chebychev" # All background is characterized by Chebychev polynomials
    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    if particle=="all":
        for p in arguments:
            resonance_fit(infile, p)
    elif particle not in arguments:
        logging.error("Invalid argument! \nPossible arguments are:\"eta\", "
            "\"rho\", \"omega\",\"phi\", \"J-psi\", \"psi'\", \"Y\", \"Z\" or "
            "\"all\"")
        inp = input("Insert the right string or leave it empty to exit the "
            "program:")
        if inp=="":
            sys.exit(1)
        else:
            logging.info(f"The input chosen is: \"{inp}\"")
            resonance_fit(infile, inp)
            error_flag=0;

    if particle!="all" and error_flag:
        logger.info(f"Fit of particle {particle}...")
        logger.debug("No error flag")

        # Mass limits are defined in order to make a selection on data for each
        # particle. The shapes of singal and background are choosed, too.
        if particle=="Y":
            mean1 = utils.FIT_INIT_PARAM["Y1"][0]
            mean2 = utils.FIT_INIT_PARAM["Y2"][0]
            mean3 = utils.FIT_INIT_PARAM["Y3"][0]
            sigma1 = utils.FIT_INIT_PARAM["Y1"][1]
            sigma2 = utils.FIT_INIT_PARAM["Y2"][1]
            sigma3 = utils.FIT_INIT_PARAM["Y3"][1]
        else:
            mean = utils.FIT_INIT_PARAM[f"{particle}"][0]
            sigma = utils.FIT_INIT_PARAM[f"{particle}"][1]

        sig = utils.FIT_INIT_PARAM[f"{particle}"][2]
        nparam = utils.FIT_INIT_PARAM[f"{particle}"][3]
        xmax = utils.FIT_INIT_PARAM[f"{particle}"][4]
        ymax = utils.FIT_INIT_PARAM[f"{particle}"][5]
        name = utils.FIT_INIT_PARAM[f"{particle}"][6]
        ndf = utils.FIT_INIT_PARAM[f"{particle}"][7]

        # Retrieve values of mass range for each particles
        lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][0]
        upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][1]

        Dimuon_mass = ROOT.RooRealVar("Dimuon_mass", "Dimuon_mass",
            lower_mass_edge,upper_mass_edge)
        Dimuon_pt = ROOT.RooRealVar("Dimuon_pt", "Dimuon_pt", low_edge_pt,
            upper_edge_pt)
        Dimuon_eta = ROOT.RooRealVar("Dimuon_eta", "Dimuon_eta", low_edge_eta,
            upper_edge_eta)

        # Define "signal", "background" and cuts for each particle in order
        # to fit them.
        cut = ROOT.RooFormulaVar("cut on mass, pt, eta",\
            f"(Dimuon_mass>{lower_mass_edge})&&(Dimuon_mass<{upper_mass_edge})\
            &&(Dimuon_pt>{low_edge_pt})&&(Dimuon_pt<{upper_edge_pt})\
            &&(Dimuon_eta>{low_edge_eta})&&(Dimuon_eta<{upper_edge_eta})", \
            ROOT.RooArgList(Dimuon_mass, Dimuon_pt, Dimuon_eta))
        rds = ROOT.RooDataSet("rds","rds",tree,
            ROOT.RooArgSet(Dimuon_mass, Dimuon_pt, Dimuon_eta),cut)
        logger.debug("The roodataset is created from the tree with cut on mass,"
            " pt, and eta.")
        num = rds.sumEntries()

        if sig=="gaus":
            sigf = ROOT.RooGaussian(sig, "The resonance",
                Dimuon_mass, mean, sigma)
        elif sig=="Crystal ball":
            alpha = ROOT.RooRealVar("alpha", "alpha", 1.5, 0.5, 5)
            n = ROOT.RooRealVar("n", "n", 5, 0., 10)
            sigf = ROOT.RooCBShape(sig, sig, Dimuon_mass, mean, sigma, alpha, n)
        elif sig=="3gaus":
            sig1 = ROOT.RooGaussian(f"{sig}1", "The resonance1",
                Dimuon_mass, mean1, sigma1)
            sig2 = ROOT.RooGaussian(f"{sig}2", "The resonance2",
                Dimuon_mass, mean2, sigma2)
            sig3 = ROOT.RooGaussian(f"{sig}3", "The resonance3",
                Dimuon_mass, mean3, sigma3)

        # Background
        a0 = ROOT.RooRealVar("a0", "a0", -0.17, -1, 1)
        a1 = ROOT.RooRealVar("a1", "a1", 0.4, -1, 1)
        if nparam == 2:
            bkgf= ROOT.RooChebychev(bkg, bkg, Dimuon_mass,
                ROOT.RooArgList(a0, a1))
        elif nparam == 3:
            a2 = ROOT.RooRealVar("a2", "a2", -1, 1)
            bkgf= ROOT.RooChebychev(bkg, bkg, Dimuon_mass,
                ROOT.RooArgList(a0, a1, a2))

        #Extended Fit
        nbkg = ROOT.RooRealVar("nbkg", "backgrounds events", num/2, 0, num)

        if particle=="Y":
            nsig1 = ROOT.RooRealVar("nsig1", "signal events", num/6, 0, num)
            nsig2 = ROOT.RooRealVar("nsig2", "signal events", num/7, 0, num)
            nsig3 = ROOT.RooRealVar("nsig3", "signal events", num/7, 0, num)
            tot = ROOT.RooAddPdf("Model", "The total pdf",
                ROOT.RooArgList(sig1, sig2, sig3, bkgf),
                ROOT.RooArgList(nsig1, nsig2, nsig3, nbkg))
        else:
            nsig = ROOT.RooRealVar("nsig", "signal events", num/5, 0, num)
            tot = ROOT.RooAddPdf("Model", "The total pdf",
                ROOT.RooArgList(sigf, bkgf), ROOT.RooArgList(nsig, nbkg))

        # Fit properties
        opt_list = ROOT.RooLinkedList()
        opt_list.Add(ROOT.RooFit.Extended(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.Save())
        opt_list.Add(ROOT.RooFit.BatchMode(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.NumCPU(4))
        #opt_list.Add(ROOT.RooFit.Timer(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.PrintLevel(-1))

        # Fit
        results = tot.fitTo(rds,opt_list)

        #Plot and styling
        mframe = Dimuon_mass.frame()
        mframe.SetTitle(f"{name} mass")
        rds.plotOn(mframe,ROOT.RooFit.Name("Data"),ROOT.RooFit.MarkerColor(861),
            ROOT.RooFit.Name("Data"))
        tot.plotOn(mframe, ROOT.RooFit.Name("Model"),ROOT.RooFit.LineColor(798))
        tot.plotOn(mframe, ROOT.RooFit.Name(f"{bkg}"),ROOT.RooFit.Components(bkg),
            ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

        mframe.GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
        mframe.GetXaxis().SetTitleSize(0.045)
        mframe.GetXaxis().CenterTitle()
        mframe.GetYaxis().CenterTitle()

        c = ROOT.TCanvas(f"{name} Resonance", f"{name} Resonance")

        # Retrieve the Chi Square
        chi = mframe.chiSquare("Model","Data", ndf)
        chi2 = f"chi2 = {round(chi,2)}"

        if particle =="Y":
            for i in range(1,4,1):
                tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}{i}"),
                    ROOT.RooFit.Components(f"{sig}{i}"),
                    ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))

            para = ROOT.RooArgSet(mean1, sigma1, mean2, sigma2, mean3, sigma3)
            tot.paramOn(mframe,ROOT.RooFit.Parameters(para),ROOT.RooFit.Label(chi2),
                ROOT.RooFit.Format("NEU",
                ROOT.RooFit.AutoPrecision(0)),ROOT.RooFit.Layout(0.62, 0.9, 0.9))
            mframe.getAttText().SetTextSize(0.033)

            #Labels
            ROOT.gPad.SetGrid()
            mframe.Draw()
            label = ROOT.TLatex()
            label.SetNDC(True)
            label.DrawLatex(0.3, 0.8, "Y(1S)")
            label.DrawLatex(0.5, 0.55, "Y(2S)")
            label.DrawLatex(0.72, 0.4, "Y(3S)")

            # Legend
            leg = ROOT.TLegend(xmax-.3, ymax+.1, xmax-.12, ymax)
            leg.AddEntry(f"{sig}1",f"{sig}", "l")

            # Print fit results
            print(f"The mass and the width of the particles {name} obtained from"
            f"the fit are: \
            \nm(Y1S) = {mean1.getValV()}; width(Y1S) = {sigma1.getValV()}\
            \nm(Y2S) = {mean2.getValV()}; width(Y2S) = {sigma2.getValV()}\
            \nm(Y3S) = {mean3.getValV()}; width(Y3S) = {sigma3.getValV()}\
            \nchi2\ndf = {chi}\n")
        else:
            tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}"),
                ROOT.RooFit.Components(sig),
                ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))

            para = ROOT.RooArgSet(mean, sigma)
            tot.paramOn(mframe,ROOT.RooFit.Parameters(para),ROOT.RooFit.Label(chi2),
                ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(0)),
                ROOT.RooFit.Layout(xmax-.27, xmax, ymax))
            mframe.getAttText().SetTextSize(0.03)
            ROOT.gPad.SetGrid()
            mframe.Draw()

            # Legend
            leg = ROOT.TLegend(xmax-.27,ymax+.1,xmax-.1,ymax)
            leg.AddEntry(f"{sig}",f"{sig}", "l")

            # Print fit results
            print(f"The mass and the width of the particle {particle} obtained "
                f"from the fit are: \
                \nm = {mean.getValV()}; width = {sigma.getValV()}\
                \nchi2\ndf = {chi}\n")

        leg.AddEntry(f"{bkg}",f"{bkg}", "l")
        leg.AddEntry("Model","Total", "l")
        leg.SetShadowColor(0)
        leg.Draw()

        #Save results in the directory "Fit"
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Fit')))
        c.SaveAs(f"{particle}_fit.pdf")
        c.SaveAs(f"{particle}_fit.png")
        utils.write_fitresults(results, f"res_{particle}.txt")

        # Return in main directory
        os.chdir(os.path.dirname(os. getcwd()))
        f.Close()


def resonance_prop(infile, mu_cached=None, particle="all"):
    """The function creates different plots of the main properties of the
    resonances such as transverse momentum, pseudorapidity and azimuthal angle.
    It takes in input :
    - the root data file or the data cached obtained by the function "
        leptons_analysis" named "dimuon.root" which contains pt, phi, eta and
        the invariant mass of the dimuons;
    - a string with the name of the resonance; the default value is "all", in
        this case plots of all resonances' properties are created.
    """

    # Cuts on pt
    pt_cut_inf = 0
    pt_cut_sup = 120

    error_flag = 1

    # Rdataframe with useful cuts is created.
    if mu_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_cut=rdf.Filter(f"(Dimuon_pt>{pt_cut_inf})&&(Dimuon_pt<{pt_cut_sup})")
        del rdf
    else:
        rdf_cut=mu_cached.Filter(f"(Dimuon_pt>{pt_cut_inf})&&(Dimuon_pt<{pt_cut_sup})")

    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    # Some check on input arguments
    if particle=="all":
        for p in arguments:
            resonance_prop(infile,mu_cached, p)
    elif particle not in arguments:
        logging.error("Invalid argument! \nPossible arguments are:\"eta\", "
            "\"rho\", \"omega\",\"phi\", \"J-psi\", \"psi'\", \"Y\", \"Z\" or "
            "\"all\"")
        inp = input("Insert the right string or leave it empty to exit the "
            "program:")
        if inp=="":
            sys.exit(1)
        else:
            logging.info(f"The particle chosen is: \"{inp}\"")
            resonance_prop(infile,mu_cached, inp)
            error_flag=0

    if particle!="all" and error_flag:

        logging.info(f"Plot properties of particle {particle}...")

        # Retrieve values of mass range for each particle
        lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][0]
        upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][1]
        rdf_m = rdf_cut.Filter(f"(Dimuon_mass>={lower_mass_edge})\
            &&(Dimuon_mass<={upper_mass_edge})", f"{particle} cut")
        del rdf_cut

        # Styling
        name = utils.FIT_INIT_PARAM[f"{particle}"][6]
        c = ROOT.TCanvas()
        c.UseCurrentStyle()
        c.SetGrid()
        ROOT.gStyle.SetOptStat("e")

        eta_lim=3.5

        # Properties for the three Y resonances
        if particle=="Y":
            eta_lim=5
            for i in range(1, 4, 1):
                # Retrieve values of mass range and do the cuts
                lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"Y{i}"][0]
                upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"Y{i}"][1]
                rdf_m_y = rdf_m.Filter(f"(Dimuon_mass>={lower_mass_edge})\
                    &&(Dimuon_mass<={upper_mass_edge})", f"Y cut{i}S")

                n = rdf_m_y.Count().GetValue()
                nbin_y = math.floor(ROOT.TMath.Sqrt(n))

                # Create a dictonary with all histograms booked
                h_y = {}
                h_y["Transverse momentum"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Transverse momentum Y({i}S)",f"Transverse momentum Y({i}S);"
                    "p_{t} [MeV];Events", nbin_y, 0, 120), "Dimuon_pt")
                h_y["Pseudorapidity"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Pseudorapidity Y({i}S)", f"Pseudorapidity Y({i}S);#eta;"
                    "Events", nbin_y, -eta_lim, eta_lim), "Dimuon_eta")
                h_y["Azimuthal angle"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Azimuthal angle Y({i}S)",f"Azimuthal angle Y({i}S);#phi "
                    "[rad];Events", nbin_y, -3.5, 3.5), "Dimuon_phi")

            # Do and save the histogram in the directory "Properties"
                for j, (k, h) in enumerate(h_y.items()):
                    h.GetXaxis().CenterTitle()
                    h.GetYaxis().CenterTitle()
                    h.GetXaxis().SetTitleSize(0.04)
                    h.SetFillColorAlpha(9, .8)
                    if j==0:
                        os.chdir(os.path.abspath(os.path.join(os.sep,
                            f'{os.getcwd()}', 'Properties')))
                        c.Print(f"Y({i}S) properties.pdf[", "pdf")
                        logger.debug("In Properties")
                    logger.debug(f"Iteration on histograms number {j}")
                    h.Draw()
                    c.Print(f"Y({i}S) properties.pdf", f"Title:Y({i}S) {k}")
                    c.SaveAs(f"Y({i}S) {k}.png")
                c.Print(f"Y({i}S) properties.pdf]", "pdf")

                # Return in main directory
                os.chdir(os.path.dirname(os. getcwd()))

        else:
            if particle=="Z":
                eta_lim=10
            # Plots for the other particles
            n = rdf_m.Count().GetValue()
            nbin = math.floor(ROOT.TMath.Sqrt(n))

            # Create a dictonary with all histograms booked
            h_p = {}
            h_p["Transverse momentum"] = rdf_m.Histo1D(ROOT.RDF.TH1DModel(
                f"Transverse momentum {name}",f"Transverse momentum {name};"
                "p_{t} [MeV];Events", nbin, 0, 120), "Dimuon_pt")
            h_p ["Pseudorapidity"]= rdf_m.Histo1D(ROOT.RDF.TH1DModel(
                f"Pseudorapidity {name}",f"Pseudorapidity {name};#eta;Events",
                nbin, -eta_lim, eta_lim), "Dimuon_eta")
            h_p["Azimuthal angle"] = rdf_m.Histo1D(ROOT.RDF.TH1DModel(
                f"Azimuthal angle {name}",f"Azimuthal angle {name};#phi [rad];"
                "Events", nbin, -3.5, 3.5), "Dimuon_phi")

            # Do and save the histograms in the directory "Properties"
            for j, (k,h) in enumerate(h_p.items()):
                h.GetXaxis().CenterTitle()
                h.GetYaxis().CenterTitle()
                h.GetXaxis().SetTitleSize(0.04)
                h.SetFillColorAlpha(9, .8)
                if j==0:
                    os.chdir(os.path.abspath(os.path.join(os.sep,
                        f'{os.getcwd()}','Properties')))
                    c.Print(f"{particle} properties.pdf[", "pdf")
                    logger.debug("In Properties")
                h.Draw()
                c.Print(f"{particle} properties.pdf", f"Title:{particle} {k}")
                c.SaveAs(f"{particle} {k}.png")
            c.Print(f"{particle} properties.pdf]", "pdf")

            # Return in main directory
            os.chdir(os.path.dirname(os. getcwd()))

# Z Analysis

# Mass range from article : 60 < M < 120:
    # 60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120;
# Eta bins of equal size for |yll| < 2.4:
    # 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4

def mass_vs_eta(infile, rap_inf, rap_sup, bin, data_cached):

    # Cuts on pt
    pt_inf = 80
    pt_sup = 120

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        branch = t_name.capitalize()
        rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_inf})&&({branch}_y<{rap_sup})").Filter(
            f"({branch}_pt>{pt_inf})&&({branch}_pt<{pt_sup})")
        del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_inf})&&({branch}_y<{rap_sup})").Filter(
            f"({branch}_pt>{pt_inf})&&({branch}_pt<{pt_sup}")
    logger.info("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel(f"Mass in defined range rapidity",
        f"Mass in defined range rapidity;Mass [MeV];Events", bin, 60, 120),
        f"{branch}_mass")

    # Styling
    c = ROOT.TCanvas("mass", "mass")
    c.SetGrid()
    ROOT.gPad.SetLogy()
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"{branch}_mass_y({rap_inf},{rap_sup})_({pt_inf},{pt_sup}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def cos_vs_eta(infile, rap_inf, rap_sup, bin, data_cached):

    # Cuts on pt
    pt_inf = 10
    pt_sup = 100

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        branch = t_name.capitalize()
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_inf})&&({branch}_y<{rap_sup})").\
            Filter(f"({branch}_pt>{pt_inf})&&({branch}_pt<{pt_sup})")
        del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
        f"({branch}_y>{rap_inf})&&({branch}_y<{rap_sup})").\
        Filter(f"({branch}_pt>{pt_inf})&&({branch}_pt<{pt_sup})")
    logger.info("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel(f"Cos in defined range rapidity",
        f"Mass in defined range rapidity;Cos;Events",bin,-1,1),f"{branch}_cos")

    # Styling
    c = ROOT.TCanvas("cos", "cos")
    c.SetGrid()
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"{branch}_cos_y({eta_inf},{eta_sup}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


#no
def cos_vs_pt(infile, eta_inf, eta_sup, bin):

    t_name = infile.replace(".root", "")
    branch = t_name.capitalize()
    rdf = ROOT.RDataFrame(t_name,infile)
    rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&({branch}_y>{eta_inf})&&({branch}_y<{eta_sup})").\
        Filter(f"({branch}_pt>10)&&({branch}_pt<100)")
    logger.info("The RDataFrame is created.")

    h = rdf_m.Graph(f"{branch}_cos", f"{branch}_pt")

    c = ROOT.TCanvas("cos", "cos")
    c.SetGrid()
    h.Draw("AP")
    c.SaveAs(f"{branch}_cos_pt({eta_inf},{eta_sup}).png")


def weight(infile, data_cached=None):

    # Cuts on mass, we are intersted in Z resonance and sidebands
    m_inf = 60
    m_sup = 120

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        branch = t_name.capitalize()
        rdf = ROOT.RDataFrame(t_name,infile).\
            Filter(f"({branch}_mass>{m_inf})&&({branch}_mass<{m_sup})")
    else:
        rdf=data_cached.Filter(f"({branch}_mass>{m_inf})&&({branch}_mass<{m_sup})")
    logger.info("The RDataFrame is created.")

    # Using JITting to define a function in order to calculate useful quantities
    # for further analysis.
    cppcode = """
    //The function takes in input the fourvector of each muon and it returns a
    //"VecOps" vector with the invariant mass, pt, phi and eta of the dimuon.

    //Function definition

    ROOT::VecOps::RVec<float> w(float pt, float m, float c)
    {
        float A0, h, den, w_d, w_n;
        A0 = (pow(pt,2))/(pow(m, 2)+pow(pt,2));
        h = (0.5)*A0*(1-(3*pow(c,2)));
        den = 1+pow(c,2)+h;
        w_d = (0.5)*(pow(c,2)/pow(den,3));
        w_n = (0.5)*(abs(c)/pow(den,2));

        ROOT::VecOps::RVec<float> P{w_d,w_n};
    	return P;
    }
    """

    ROOT.gInterpreter.ProcessLine(cppcode)                                      #It compiles the C++ code
    logger.debug("The c++ function is defined.")

    rdf_w = rdf.\
        Define("w_d", f"w({branch}_pt,{branch}_mass, {branch}_cos)[0]").\
        Define("w_n", f"w({branch}_pt,{branch}_mass, {branch}_cos)[1]")
    logger.debug("The new columns are defined.")
    del rdf

    # A list with the names of the new columns is created in order to make a
    # snapshot of them
    branchlist = [f"{branch}_mass", f"{branch}_pt", f"{branch}_eta", f"{branch}_phi",\
        f"{branch}_cos",f"{branch}_y","w_d", "w_n"]

    # A snapshot is done to collect the useful physiscal quantity of the dileptons
    # in a single root file.
    logger.debug("Starting snapshot...")
    rdf_w.Snapshot(f"{t_name}_w", f"{t_name}_w.root", branchlist)
    logger.info("The snapshot is done.")

def afb(infile):

    # Cuts on pt (to reduce backgorund?)
    pt_inf = 10
    pt_sup = 120

    # Create RDataFrame
    t_name = infile.replace(".root", "")
    branch = t_name.replace("_w", "").capitalize()
    rdf = ROOT.RDataFrame(t_name, infile).Filter(f"({branch}_pt>{pt_inf})&&"
        f"({branch}_pt<{pt_sup})")
    logger.info("The RDataFrame has been created.")

    #
    for j in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve rapidity range bin from dictionary RAPIDITY_BIN
        y_inf = utils.RAPIDITY_BIN[f"{j}"][0]
        y_sup = utils.RAPIDITY_BIN[f"{j}"][1]

        rdf_y = rdf.Filter(f"({branch}_y>{y_inf})&&({branch}_y<{y_sup})")
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))
        with open(f"Afb[y({y_inf},{y_sup})]_({pt_inf},{pt_sup}){t_name}.txt", "w+") as outf:
            print("Afb:M:A4", file=outf)
            os.chdir(os.path.dirname(os. getcwd()))
            for i in range(0, len(utils.MASS_BIN), 1):
                # Retrive mass range bin from dictonary MASS_BIN
                m_inf = utils.MASS_BIN[f"{i}"][0]
                m_sup = utils.MASS_BIN[f"{i}"][1]

                rdf_m = rdf_y.Filter(f"({branch}_mass>{m_inf})&&({branch}_mass<{m_sup})")
                rdf_f = rdf_m.Filter(f"{branch}_cos>0")
                D_f = rdf_f.Sum("w_d")
                N_f = rdf_f.Sum("w_n")
                rdf_b = rdf_m.Filter(f"{branch}_cos<0")
                D_b = rdf_b.Sum("w_d")
                N_b = rdf_b.Sum("w_n")

                Afb = (3*(N_f.GetValue()-N_b.GetValue()))/(8*(D_f.GetValue()+D_b.GetValue()))
                M = m_inf +((m_sup-m_inf)/2)
                A4 = Afb*(8/3)

                print(Afb, M, A4,file=outf)
    del rdf_y, rdf_m, rdf_f, rdf_b

def afb_plot(infile):
    with open(infile) as ft:
        data = [_.rstrip('\n') for _ in ft]
    c = ROOT.TCanvas("Afb")
    file_name = infile.replace("txt", "pdf")
    c.Print(f"{file_name}[")
    for inf in data:
        t = ROOT.TTree("tree", "tree")
        t.ReadFile(inf)

        c.UseCurrentStyle()
        c.SetGrid()
        ROOT.gStyle.SetOptStat("e")
        t.Draw("Afb:M","")
        name = inf.replace("txt", "")
        c.Print(f"{file_name}", f"Title:{name}")
    c.Print(f"{file_name}]")#, #f"Title:{name}")



# Fit mass, pt, phi, eta from Z on DIMUONS
def fit_pt(infile):

    # Create the Ttree for the analysis
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    # Histogram of pt in mass range 60 - 120 GeV (Z boson)
    c = ROOT.TCanvas("Dimuon_pt", "Dimuon_pt")
    tree.Draw("Dimuon_pt>>h(100,0,120)", "(Dimuon_mass>60)&&(Dimuon_mass<120)")
    h_pt = ROOT.gPad.GetPrimitive("h")
    h_pt.Draw()
    results = h_pt.Fit("landau", "S")
    MPV = results.Parameter(1)
    sigma = results.Parameter(2)
    with open("pt_Z_result.txt","w+") as outf:
        print("Results of Z transverse momentum distribution fit:", file=outf)
        print(f"MPV = {MPV}", file=outf)
        print(f"Sigma = {sigma}", file=outf)
    c.SaveAs("pt_Z.png")

def fit_eta(infile):

    # Create the Ttree for the analysis
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    # Histogram of pt in mass range 60 - 120 GeV (Z boson)
    c = ROOT.TCanvas("Dimuon_eta", "Dimuon_eta")
    tree.Draw("Dimuon_eta>>h(100,-10,10)", "(Dimuon_mass>60)&&(Dimuon_mass<120)")
    h_eta = ROOT.gPad.GetPrimitive("h")
    h_eta.Draw()
    f1 = ROOT.TF1("f1", "gaus", -10, 0)
    f2 = ROOT.TF1("f2", "gaus", 0,10)
    h_eta.Fit(f1,"R")
    h_eta.Fit(f2,"R")
    par = np.zeros(6)
    f1.GetParameters(par[0:2])
    f2.GetParameters(par[3:5])

    total = ROOT.TF1("total", "gaus(0)+gaus(3)", -10,10)
    total.SetParameters(par)
    h_eta.Fit(total,"R")
    total.GetParameters(par)

    with open("eta_Z_result.txt","w+") as outf:
        print("Results of Z pseudorapidity distribution fit:", file=outf)
        print(f"Mean1 = {par[1]}", file=outf)
        print(f"Sigma1 = {par[2]}", file=outf)
        print(f"Mean1 = {par[4]}", file=outf)
        print(f"Sigma2 = {par[5]}", file=outf)
    c.SaveAs("eta_Z.png")

def fit_mass(infile):

    # Create the Ttree for the analysis
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    # Histogram of pt in mass range 60 - 120 GeV (Z boson)
    c = ROOT.TCanvas("Dimuon_mass", "Dimuon_mass")
    tree.Draw("Dimuon_mass>>h(100,60,120)", "(Dimuon_mass>60)&&(Dimuon_mass<120)")
    h_mass = ROOT.gPad.GetPrimitive("h")
    h_mass.Draw()
    f1 = ROOT.TF1("gaus", "gaus", 82,96)
    f2 = ROOT.TF1("expo","expo",60,120)
    h_mass.Fit(f1,"R")
    h_mass.Fit(f2,"R")
    par = np.zeros(5)
    f1.GetParameters(par[0:2])
    f2.GetParameters(par[3:4])

    total = ROOT.TF1("total", "gaus(0)+expo(3)", 60,120)
    total.SetParameters(par)
    h_mass.Fit(total,"R")
    total.GetParameters(par)

    with open("mass_Z_result.txt","w+") as outf:
        print("Results of Z mass distribution fit:", file=outf)
        print(f"Mean = {par[1]}", file=outf)
        print(f"Sigma = {par[2]}", file=outf)
        print(f"Const = {par[3]}", file=outf)
        print(f"Slope = {par[4]}", file=outf)
    c.SaveAs("mass_Z.png")

def fit_cos(infile):
    # Create the Ttree for the analysis
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    # Histogram of pt in mass range 60 - 120 GeV (Z boson)
    c = ROOT.TCanvas("Dimuon_mass", "Dimuon_mass")
    tree.Draw("Dimuon_cos>>h(100,-1,1)", "(Dimuon_mass>86)&&(Dimuon_mass<94)&&(Dimuon_pt>60)&&(Dimuon_pt<100)&&(Dimuon_eta>-2.4)&&(Dimuon_eta<2.4)")
    h_mass = ROOT.gPad.GetPrimitive("h")
    h_mass.Draw()
    # f1 = ROOT.TF1("gaus", "gaus", 82,96)
    # f2 = ROOT.TF1("expo","expo",60,120)
    # h_mass.Fit(f1,"R")
    # h_mass.Fit(f2,"R")
    # par = np.zeros(5)
    # f1.GetParameters(par[0:2])
    # f2.GetParameters(par[3:4])
    #
    # total = ROOT.TF1("total", "gaus(0)+expo(3)", 60,120)
    # total.SetParameters(par)
    # h_mass.Fit(total,"R")
    # total.GetParameters(par)
    #
    # with open("mass_Z_result.txt","w+") as outf:
    #     print("Results of Z mass distribution fit:", file=outf)
    #     print(f"Mean = {par[1]}", file=outf)
    #     print(f"Sigma = {par[2]}", file=outf)
    #     print(f"Const = {par[3]}", file=outf)
    #     print(f"Slope = {par[4]}", file=outf)
    c.SaveAs("cos_Z.png")

# Try to retrieve distribution of single muons
def test_begin(infile ,col):

    # Con TTree (per i tagli ci mette molto)

    # # Create the Ttree for the analysis
    # f = ROOT.TFile.Open(infile)
    # tree = f.Get(f"Events")
    # logger.debug("The tree is created.")
    # c = ROOT.TCanvas("", "")
    #tree.Draw(f"{col}>>h(100,0,120)","(nMuon==2)&&(Muon_charge[0]!=Muon_charge[1])")
    # h_mass = ROOT.gPad.GetPrimitive("h")
    # h_mass.Draw()
    # c.SaveAs(f"{col}.png")

    # Con RDataFrame (molto più veloce sui tagli)
    root_file = ROOT.TFile.Open(infile, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    rdf_cut=rdf.Filter(f"(nMuon==2)&&(Muon_charge[0]!=Muon_charge[1])")

    if col=="Muon_pt":
        h=rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"{col}", f"{col}", 200, 0, 120), f"{col}")
    elif col=="Muon_eta":
        h=rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"{col}", f"{col}", 200, -10, 10), f"{col}")
    elif col=="Muon_phi":
        h=rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"{col}", f"{col}", 200, -3.4, 3.4), f"{col}")
    elif col=="Muon_mass":
        h=rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"{col}", f"{col}", 10, 0.1, .2), f"{col}")


    c = ROOT.TCanvas("", "")
    h.Draw()
    c.SaveAs(f"{col}.png")





if __name__ == "__main__":

    # Set ROOT environment
    ROOT.gROOT.SetBatch()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)

    # Enable the multi-threading analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    # Style of the canvas
    myStyle = ROOT.TStyle('My Style', 'My Style')
    myStyle.SetCanvasColor(0)
    myStyle.SetMarkerSize(1.2)
    myStyle.SetMarkerColor(861)
    myStyle.SetMarkerStyle(20)
    myStyle.SetHistLineColor(9)
    myStyle.SetLabelFont(42)
    myStyle.SetTitleFont(22)
    myStyle.SetLabelSize(0.035)
    ROOT.gROOT.SetStyle("My Style")
    myStyle.cd()

    # Start the the timer
    timer = ROOT.TStopwatch()
    timer.Start()

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str,
        help="The path of the nanoAOD file to analyze.")
    parser.add_argument("-p", "--particle",required=False, type=str,
        help="The name of the particle to analyze. It is necessary if you want \
            to see resonances' fit and properties. Possible arguments are:\
            \"eta\", \"rho\",\"omega\", \"phi\", \"J-psi\", \"psi'\", \"Y\", \
            \"Z\", \"all\". Put the chosen string in quotes.")
    args = parser.parse_args()

    # Creating the logger and setting its level
    logger = utils.set_logger("Dimuon_invm")

    logger.info("Starting the analysis of the root file nanoAOD...")

    ########################################################################
    #test_begin(f"{args.file}", "Muon_pt")
    #test_begin(f"{args.file}", "Muon_eta")
    #test_begin(f"{args.file}", "Muon_phi")
    #test_begin(f"{args.file}", "Muon_mass")


    # LEPTONS ANALYSIS

    # If "leptons_analysis" is not run (it already has been done), dimu_cached
    # and diel_cached are set to None, so in the following the other functions
    # take data directly from tha data saved by the snapshot.
    # But "leptons_analysis" must be run at least one to collect data.
    dimu_cached = None
    diel_cached = None
    outfile_m = "dimuon.root"
    outfile_el = "dielectron.root"
    #dimu_cached, diel_cached = leptons_analysis(f"{args.file}")

    # DIMUON MASS SPECTRUM
    os.makedirs("Spectra", exist_ok=True)
    logger.debug("The new directory \"Spectra\" is created")
    #mumu_spectrum(outfile_m, dimu_cached)

    # DIMUON MASS SPECTRUM WITHOUT BUMP
    #mumu_spectrum_bump(outfile_m, dimu_cached)

    # DIMUON ETA DISTRIBUTION
    #mumu_eta(outfile_m, dimu_cached)

    # RESONANCES' FIT
    os.makedirs("Fit", exist_ok=True)
    logger.debug("The new directory \"Fit\" is created")
    #resonance_fit(outfile_m, f"{args.particle}")

    # PROPERTIES
    os.makedirs("Properties", exist_ok=True)
    logger.debug("The new directory \"Properties\" is created")
    # resonance_prop(outfile_m, dimu_cached, f"{args.particle}")

    # Z ANALYSIS
    os.makedirs("Z analysis", exist_ok=True)
    logger.debug("The new directory \"Z analysis\" is created")

    #############################################################
    # Auxiliary
    # fit_pt("dimuon.root")
    # fit_eta("dimuon.root")
    # fit_mass("dimuon.root")
    # fit_cos("dimuon.root")

    # Retrieve number of dimuons and dielectrons in the analysis

###############################################################################

    # Plot Z resonance in different eta ranges for both dimuons and
    # dielectrons data
    # data = ["dimuon", "dielectron"]
    # cached = [dimu_cached, diel_cached]
    # for i in range(0, len(utils.RAPIDITY_BIN), 1):
    #
    #     # Retrieve the value of eta range bin
    #     y_inf = utils.RAPIDITY_BIN[f"{i}"][0]
    #     y_sup = utils.RAPIDITY_BIN[f"{i}"][1]
    #
    #     bin = 40 #?
    #     for dilepton, cache in zip(data, cached):
    #         #mass_vs_eta(f"{dilepton}.root", y_inf, y_sup, bin, cache)
    #         #cos_vs_eta(f"{dilepton}.root", y_inf, y_sup, bin, cache)
    #         #mass_vs_eta(f"dimuon.root", y_inf, y_sup, bin, cache)
    #         #cos_vs_eta("dimuon.root", y_inf, y_sup, bin, cache)
    #         pass

    # AFB
    #weight("dimuon.root")
    # afb("dimuon_w.root")
    #
    # # Plot AFB for different cut in pt
    #
    # os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))
    # pt_plot = ["(10,120)", "(80,120)"]
    # for p in pt_plot:
    #     pathlist = Path("").glob(f"*_{p}dimuon_w.txt")
    #     with open(f"Afblist_{p}.txt", "w+") as outf:
    #         for path in pathlist:
    #             print(path, file=outf)
    #
    # ##Dimuon Afb plots
    # for p in pt_plot:
    #     afb_plot(f"Afblist_{p}.txt")
    # os.chdir(os.path.dirname(os. getcwd()))

###############################################################################


    # More inputs
    # inputs = []
    # if args.particle:
    #     for p in args.particle:
    #         inputs.extend()
    # else:
    #     logger.info("No particle in input, the default value is \"all\"".)
    #
    # if inputs:
    #     for particle in tdqm(inputs):
    #         try:
    #             main(particle)
    #         except:
    #             .....
    # else:
    #     main("all")
    #


    # Elapsed time
    timer.Stop()
    logger.info(f"Elapsed time from the beginning is: {timer.RealTime()}(real time) seconds")
    logger.info(f"Elapsed time from the beginning is: {timer.CpuTime()}(cpu time) seconds")
