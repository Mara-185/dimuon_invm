import ROOT
import argparse
import logging
import os
import math
import sys
from pathlib import Path

"""The script takes as argument:
    - the data file (URL) of dileptons (-f), for example:
      root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root";
    - the string of the particle's name (-p) to analyze and fit, among those in
      the dimuon spectrum, which are :
      "eta", "rho","omega", "phi", "J-psi", "psi'", "Y", "Z".

      There are different functions, among which:
      - "leptons_analysis" which selects the couple of muons and electrons which are interesting
        in order to create the dimuon mass spectrum and for further analysis;
      - "mumu_spectrum" which plot the dimuon mass spectrum.
      - "resonance_prop" creates different plots of main characteristics of the
        particle chosen in the spectrum;
      - "resonance_fit" solves every resonance and returns the plot with the fit
        and a txt with the fit results.
       - [........................]
"""


def leptons_analysis(infile):
    """ It takes in input a nano-AOD data file and returns two root files, named :
    "dimuon.root" and "dielectron.root" both with four columns, each one related
    to different quantity, for example:
    ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi"].
    It returns the data cached in memory, in order to speed up the analysis.
    """

    nthreads = 8
    ROOT.ROOT.EnableImplicitMT(nthreads)

    #Creating an RDataFrame for the analysis
    inf = ROOT.TFile.Open(infile, "READ")
    tree_name = inf.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, inf)
    logger.info("The RDataFrame is created.")

    #Filter on two muons and two electrons with opposite charge
    rdf_mu = rdf.Filter("nMuon==2","Selection of two muons").\
     Filter("Muon_charge[0]!=Muon_charge[1]","Selection of muons with opposite charge")
    logger.info("The cut on two muons with opposite charge is done.")
    rdf_e = rdf.Filter("nElectron==2","Selection of two electrons").\
     Filter("Electron_charge[0]!=Electron_charge[1]","Selection of electrons with opposite charge")
    logger.info("The cut on two electrons with opposite charge is done.")

    #By "Jitting", a function is created that returns the invariant mass, pt
    #(transverse momentum), phi (azimuthal angle) and eta(pseudorapidity) of the two muons.

    ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<double> dimu_vec(float pt0, float eta0, float phi0, float mass0, float pt1, float eta1, float phi1, float mass1)
    {
    	ROOT::Math::PtEtaPhiMVector p0, p1;			                            //????Definiamo un TLorentz Vector da Math perche' ottimizzato rispetto a quello della classe relativa (vedi bene questione annessa)
    	p0.SetCoordinates(pt0, eta0, phi0, mass0);
    	p1.SetCoordinates(pt1, eta1, phi1, mass1);
    	ROOT::VecOps::RVec<double> P{(p1+p0).Pt(), (p1+p0).Eta(), (p1+p0).Phi(), (p1+p0).M()};
    	return P;
    }


    ROOT::VecOps::RVec<float> cos_rap(float pt0, float eta0, float phi0, float mass0, float pt1, float eta1, float phi1, float mass1)
    {
        float pz0, pz1, E0, E1, P0_1, P0_2, P1_1, P1_2, cos, numer, denom, mll, ptt, pzll, y;
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

    logger.debug("The second c++ function is defined.")

    #By the new function, four new columns of data has been created from our data file.
    rdf_dimu = rdf_mu.\
        Define("Dimuon_mass", "dimu_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[3]").\
        Define("Dimuon_pt", "dimu_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_eta", "dimu_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]").\
        Define("Dimuon_phi", "dimu_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[2]").\
        Define("Dimuon_cos", "cos_rap(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_y", "cos_rap(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]")

    rdf_die = rdf_e.\
        Define("Dielectron_mass", "dimu_vec(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[3]").\
        Define("Dielectron_pt", "dimu_vec(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[0]").\
        Define("Dielectron_eta", "dimu_vec(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[1]").\
        Define("Dielectron_phi", "dimu_vec(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[2]").\
        Define("Dielectron_cos", "cos_rap(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[0]").\
        Define("Dielectron_y", "cos_rap(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1])[1]")
    logger.debug("The new columns of dimuons and dieletrons are defined.")


    #A list with the names of the new columns is created in order to make a snapshot of them
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi", "Dimuon_cos", "Dimuon_y"]
    branchlist_e = [ "Dielectron_mass", "Dielectron_pt", "Dielectron_eta", "Dielectron_phi", "Dielectron_cos", "Dielectron_y"]

    #Two snapshots are done to collect the useful physiscal quantity of the dimuons and dielectrons
    #in a single root file.
    logger.debug("Starting snapshots...")
    #opt = ROOT.RDF.RSnapshotOptions()
    #opt.fLazy = True
    rdf_die.Snapshot("dielectron", "dielectron.root", branchlist_e)
    mu_c = rdf_dimu.Snapshot("dimuon", "dimuon.root", branchlist_mu).Cache()
    #################################LAZY?
    logger.info("The snapshots are done.")

    logger.info("CutFlow muons:")
    rdf_mu.Report().Print()

    logger.info("CutFlow electrons:")
    rdf_e.Report().Print()

    inf.Close()
    return mu_c


def mumu_spectrum(infile, mu_c=None):
    """It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" named "dilepton.root" and plot the histogram of
    the dimuons invariant mass."""


    if mu_c == None:
        t_name = infile.replace(".root", "")
        rdf_m = ROOT.RDataFrame(t_name,infile)
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass", 50000, 0.3, 200), "Dimuon_mass")
    else:
        h = mu_c.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass", 50000, 0.3, 200), "Dimuon_mass")


    #Styling
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

    #Labels
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.DrawLatex(0.165, 0.720, "#eta")
    label.DrawLatex(0.190, 0.772, "#rho,#omega")
    label.DrawLatex(0.245, 0.775, "#phi")
    label.DrawLatex(0.400, 0.850, "J/#psi")
    label.DrawLatex(0.410, 0.700, "#psi'")
    label.DrawLatex(0.485, 0.700, "Y(1, 2, 3S)")
    label.DrawLatex(0.795, 0.680, "Z")

    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectrum')))
    c.SaveAs("dimuon_spectrum.pdf")
    c.SaveAs("dimuon_spectrum.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The file \"dimuon_spectrum.pdf\" has been created.")


def mumu_spectrum_bump(infile, mu_c=None):
    """It takes in input the root data file obtained by the function "leptons_analysis"
    named "dilepton.root" and plot the histogram of the dimuons invariant mass."""

    if mu_c == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass", 50000, 0.3, 200), "Dimuon_mass")
    else:
        rdf_m = mu_c.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass", 50000, 0.3, 200), "Dimuon_mass")


    #Styling
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

    #Labels
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

    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectrum')))
    c.SaveAs("dimuon_spectrum_bump.pdf")
    c.SaveAs("dimuon_spectrum_bump.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The file \"dimuon_spectrum_bump.pdf\" has been created.")


if __name__ == "__main__":

    ROOT.gROOT.SetBatch()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)

    #Enable the multi-threading analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    #Style of the canvas
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

    #Start the the timer
    timer = ROOT.TStopwatch()
    timer.Start()

    #Creating the parser
    parser = argparse.ArgumentParser(description = "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str,\
        help="The path of the nanoAOD file to analyze.")
    parser.add_argument("-c", "--cuts",required=False, type=int,\
        help="If True the cuts on pt and eta are done, otherwise the function selects only a values' window on the mass.")
    parser.add_argument("-p", "--particle",required=False, type=str, \
        help="The name of the particle to analyze. It is necessary if you want to see resonances' fit and properties.")
    args = parser.parse_args()

    #Creating the logger and setting its level
    logger = logging.getLogger("Dimuon_invm")
    logging.basicConfig(level=logging.DEBUG)
    handler = logging.FileHandler("Dimuon_invm.log", "w+")
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    logger.info("Starts the analysis of the root file nanoAOD...")


    #LEPTONS ANALYSIS
    mu_cache = None
    #mu_cache = leptons_analysis(f"{args.file}")

    #SPECTRUM
    os.makedirs("Spectrum", exist_ok=True)
    logger.debug("The new directory \"Spectrum\" is created")
    mumu_spectrum("dimuon.root", mu_cache)
    mumu_spectrum_bump("dimuon.root", mu_cache)

    #BUMP


    #print elapsed time
    timer.Stop()
    logger.info(f"Elapsed time from the beginning is: {timer.RealTime()}(real time) seconds")
