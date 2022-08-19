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
      - "resonance_fit" solves every resonance and returns the plot with the fit
        and a txt with the fit results.
      - "resonance_prop" creates different plots of main characteristics of the
        particle chosen in the spectrum;
      - [........................]

    In the analysis the following version have been used:
    - Python v3.8
    - ROOT v6.24
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
    function "leptons_analysis" named "dimuon.root" and plot the histogram of
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

    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_spectrum.pdf")
    c.SaveAs("dimuon_spectrum.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The file \"dimuon_spectrum.pdf\" has been created.")


def mumu_spectrum_bump(infile, mu_c=None):
    """It takes in input the root data file obtained by the function "leptons_analysis"
    named "dimuon.root" and plot the histogram of the dimuons invariant mass
    with a cut for events which have pt>20 GeV."""

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

    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_spectrum_bump.pdf")
    c.SaveAs("dimuon_spectrum_bump.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The file \"dimuon_spectrum_bump.pdf\" has been created.")


def mumu_eta(infile, mu_c=None):
    """It takes in input the root data file obtained by the function "leptons_analysis"
    named "dimuon.root" and plot the histogram of the dimuons pseudorapidity."""

    if mu_c == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon eta", "Dimuon eta", 150, -4, 4), "Dimuon_eta")
    else:
        rdf_m = mu_c.Filter("Dimuon_pt>20")
        h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Dimuon eta", "Dimuon eta", 150, -4, 4), "Dimuon_eta")

    #Styling
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

    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    c.SaveAs("dimuon_eta.pdf")
    c.SaveAs("dimuon_eta.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The file \"dimuon_eta.pdf\" has been created.")


def write_fitresults(results, filename):
    text = ROOT.std.ofstream(filename)
    results.printMultiline(text, 1111, True)
    text.close()


def resonance_fit(infile, particle):
    """It takes in input the root data file obtained with "mumu_analysis"
        and a string with the name of the resonance to fit.
        It retrieves a plot with the fit, a txt with fit results and a workspace
        ????????????????????????????????????????????????????????"""

    #An auxiliary TTree from the root data file is created.
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    #Styling
    ymax=0.8
    xmax=0.4
    name=f"#{particle}"
    bkg = "Chebychev"                                                              #The background is characterized by Chebychev polynomials
    sig = "gaus"
    low_edge_pt = 10
    upper_edge_pt = 100
    low_edge_eta = -1.2
    upper_edge_eta = 1.2
    nparam = 2

    error=1

    #Mass limits are defined in order to make a selection on data for each particle.
    #Also the shapes of singal and background are choosed.
    if particle=="eta":
        low_edge = 0.52
        upper_edge = 0.57
        mean = ROOT.RooRealVar("mean", "mean", 0.552, 0.53, 0.56)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.07, 0.01, 0.1)
        ymax=0.4
    elif particle=="rho" or particle=="omega":
        low_edge = 0.72
        upper_edge = 0.84
        mean = ROOT.RooRealVar("mean", "mean", 0.78, 0.75, 0.8)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1)
    elif particle=="phi":
        low_edge = 0.96
        upper_edge = 1.07
        mean = ROOT.RooRealVar("mean", "mean", 1.02, 0.98, 1.06)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.004, 0.0001, 0.1)
    elif particle=="J-psi":
        low_edge = 2.65
        upper_edge = 3.55
        mean = ROOT.RooRealVar("mean", "mean", 3.10, 2.9, 3.2)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.)
        sig = "Crystal ball"
        nparam = 3
        xmax=0.42
        name="J/#psi"
    elif particle=="psi'":
        low_edge = 3.55
        upper_edge = 3.85
        mean = ROOT.RooRealVar("mean", "mean", 3.7, 3.6, 3.8)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.04, 0.0001, 1.)
        nparam = 3
        xmax = 0.9
    elif particle=="Z":
        low_edge = 82
        upper_edge = 98
        mean = ROOT.RooRealVar("mean", "mean", 91, 89, 93)
        sigma = ROOT.RooRealVar("sigma", "sigma",2, 0.01, 4 )
        xmax=0.42
        name=f"{particle}"
    elif particle=="Y":
        low_edge = 8.5
        upper_edge = 11
        mean1 = ROOT.RooRealVar("mean1", "mean1", 9.4, 9.2, 9.7)
        mean2 = ROOT.RooRealVar("mean2", "mean2", 10, 9.8, 10.15)
        mean3 = ROOT.RooRealVar("mean3", "mean3", 10.3, 10., 10.5)
        sigma1 = ROOT.RooRealVar("sigma1", "sigma1", 0.005, 0.001, 1.)
        sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 0.005, 0.001, 1.)
        sigma3 = ROOT.RooRealVar("sigma3", "sigma3", 0.005, 0.001, 1.)
        sig = "3gaus"
        nparam = 3
        name=f"{particle}"
    elif particle=="all":
        for p in ["eta", "rho", "phi", "J-psi", "psi'", "Y", "Z"]:
            resonance_fit(infile, p)
    else:
        print("Invalid argument!")
        print("Possible arguments are:\"eta\", \"rho\", \"omega\", \"phi\",\
         \"J-psi\", \"psi'\", \"Y\", \"Z\" or \"all\"")
        inp = input("Insert the right string or leave it empty to exit the program:")
        if inp=="":
            sys.exit(1)
        else:
            print(f"The particle chosen is: \"{inp}\"")
            resonance_fit(infile, inp)
            error=0;

    if particle!="all" and error!=0:
        Dimuon_mass = ROOT.RooRealVar("Dimuon_mass", "Dimuon_mass", low_edge, upper_edge)
        Dimuon_pt = ROOT.RooRealVar("Dimuon_pt", "Dimuon_pt", low_edge_pt, upper_edge_pt)
        Dimuon_eta = ROOT.RooRealVar("Dimuon_eta", "Dimuon_eta", low_edge_eta, upper_edge_eta)


        cut = ROOT.RooFormulaVar("cut on mass, pt, eta",\
         f"(Dimuon_mass>{low_edge})&&(Dimuon_mass<{upper_edge})&&(Dimuon_pt>{low_edge_pt})\
         &&(Dimuon_pt<{upper_edge_pt})&&(Dimuon_eta>{low_edge_eta})&&(Dimuon_eta<{upper_edge_eta})", \
         ROOT.RooArgList(Dimuon_mass, Dimuon_pt, Dimuon_eta))
        rds = ROOT.RooDataSet("rds","rds",tree,ROOT.RooArgSet(Dimuon_mass, Dimuon_pt, Dimuon_eta),cut)
        logger.debug("The roodataset is created from the tree with cut on mass, pt, and eta.")
        num = rds.sumEntries()


        if sig=="gaus":
            sigf = ROOT.RooGaussian(sig, "The resonance", Dimuon_mass, mean, sigma)
        elif sig=="Crystal ball":
            alpha = ROOT.RooRealVar("alpha", "alpha", 1.5, 0.5, 5)
            n = ROOT.RooRealVar("n", "n", 5, 0., 10)
            sigf = ROOT.RooCBShape(sig, sig, Dimuon_mass, mean, sigma, alpha, n)
        elif sig=="3gaus":
            sig1 = ROOT.RooGaussian(f"{sig}1", "The resonance1", Dimuon_mass, mean1, sigma1)
            sig2 = ROOT.RooGaussian(f"{sig}2", "The resonance2", Dimuon_mass, mean2, sigma2)
            sig3 = ROOT.RooGaussian(f"{sig}3", "The resonance3", Dimuon_mass, mean3, sigma3)

        if bkg=="Chebychev":
            a0 = ROOT.RooRealVar("a0", "a0", -0.17, -1, 1)
            a1 = ROOT.RooRealVar("a1", "a1", 0.4, -1, 1)
            if nparam == 2:
                bkgf= ROOT.RooChebychev(bkg, bkg, Dimuon_mass, ROOT.RooArgList(a0, a1))
            elif nparam == 3:
                a2 = ROOT.RooRealVar("a2", "a2", -1, 1)
                bkgf= ROOT.RooChebychev(bkg, bkg, Dimuon_mass, ROOT.RooArgList(a0, a1, a2))

        #Extended Fit
        nbkg = ROOT.RooRealVar("nbkg", "backgrounds events", num/2, 0, num)

        if particle=="Y":
            nsig1 = ROOT.RooRealVar("nsig1", "signal events", num/6, 0, num)
            nsig2 = ROOT.RooRealVar("nsig2", "signal events", num/7, 0, num)
            nsig3 = ROOT.RooRealVar("nsig3", "signal events", num/7, 0, num)
            tot = ROOT.RooAddPdf("Model", "The total pdf", ROOT.RooArgList(sig1, sig2, sig3, bkgf),\
                ROOT.RooArgList(nsig1, nsig2, nsig3, nbkg))
        else:
            nsig = ROOT.RooRealVar("nsig", "signal events", num/5, 0, num)
            tot = ROOT.RooAddPdf("tot", "The total pdf", ROOT.RooArgList(sigf, bkgf),\
                ROOT.RooArgList(nsig, nbkg))

        opt_list = ROOT.RooLinkedList()
        opt_list.Add(ROOT.RooFit.Extended(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.Save())
        opt_list.Add(ROOT.RooFit.BatchMode(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.NumCPU(4))
        #opt_list.Add(ROOT.RooFit.Timer(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.PrintLevel(-1))

        results = tot.fitTo(rds,opt_list)

        #Plot and styling
        mframe = Dimuon_mass.frame()
        mframe.SetTitle(f"{name} mass")
        rds.plotOn(mframe,ROOT.RooFit.Name("Data"),ROOT.RooFit.MarkerColor(861), ROOT.RooFit.Name("Data"))
        tot.plotOn(mframe, ROOT.RooFit.Name("Model"),ROOT.RooFit.LineColor(798))
        tot.plotOn(mframe, ROOT.RooFit.Name(f"{bkg}"),ROOT.RooFit.Components(bkg), ROOT.RooFit.LineStyle(ROOT.kDashed),\
            ROOT.RooFit.LineColor(ROOT.kRed))

        if particle =="Y":
            tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}1"), ROOT.RooFit.Components(f"{sig}1"), \
             ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))
            tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}2"), ROOT.RooFit.Components(f"{sig}2"), \
             ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))
            tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}3"), ROOT.RooFit.Components(f"{sig}3"), \
             ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))
        else:
            tot.plotOn(mframe, ROOT.RooFit.Name(f"{sig}"), ROOT.RooFit.Components(sig), \
             ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(418))

        c = ROOT.TCanvas(f"{name} Resonance", f"{name} Resonance")
        chi = mframe.chiSquare()

        mframe.GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
        mframe.GetXaxis().SetTitleSize(0.045)
        mframe.GetXaxis().CenterTitle()
        mframe.GetYaxis().CenterTitle()

        #Deaw box and Legend
        if particle=="Y":
            para = ROOT.RooArgSet(mean1, sigma1, mean2, sigma2, mean3, sigma3, ROOT.RooRealVar("chi", "chi",chi))
            tot.paramOn(mframe,ROOT.RooFit.Parameters(para), ROOT.RooFit.Format("NEU",\
                ROOT.RooFit.AutoPrecision(0)),ROOT.RooFit.Layout(0.6, 0.9, 0.9))
            mframe.getAttText().SetTextSize(0.035)
            #Labels
            ROOT.gPad.SetGrid()
            mframe.Draw()
            label = ROOT.TLatex()
            label.SetNDC(True)
            label.DrawLatex(0.3, 0.8, "Y(1S)")
            label.DrawLatex(0.5, 0.55, "Y(2S)")
            label.DrawLatex(0.7, 0.5, "Y(3S)")
        else:
            para = ROOT.RooArgSet(mean, sigma, ROOT.RooRealVar("chi", "chi",chi))
            tot.paramOn(mframe,ROOT.RooFit.Parameters(para), ROOT.RooFit.Format("NEU",\
                ROOT.RooFit.AutoPrecision(0)), ROOT.RooFit.Layout(xmax-.3, xmax, ymax))
            mframe.getAttText().SetTextSize(0.035)
            ROOT.gPad.SetGrid()
            mframe.Draw()

        if particle=="Y":
            leg = ROOT.TLegend(xmax-.3, ymax+.1, xmax-.1, ymax)
            leg.AddEntry(f"{sig}1",f"{sig}", "l")
        else:
            leg = ROOT.TLegend(xmax-.3,ymax+.1,xmax-.1,ymax)
            leg.AddEntry(f"{sig}",f"{sig}", "l")
        leg.AddEntry(f"{bkg}",f"{bkg}", "l")
        leg.AddEntry("Model","Total", "l")
        leg.SetShadowColor(0)
        leg.Draw()

        #Save results in directory "Fit"
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Fit')))
        c.SaveAs(f"{particle}_fit.pdf")
        c.SaveAs(f"{particle}_fit.png")
        write_fitresults(results, f"res_{particle}.txt")

        if particle=="Y":
            print(f"The mass and the width of the particles {name} obtained from the fit are:\
            \nm(Y1S) = {mean1.getValV()}; width(Y1S) = {sigma1.getValV()}\
            \nm(Y2S) = {mean2.getValV()}; width(Y2S) = {sigma2.getValV()}\
            \nm(Y3S) = {mean3.getValV()}; width(Y3S) = {sigma3.getValV()}\nchi2 = {chi}\n")
        else:
            print(f"The mass and the width of the particle {particle} obtained from the fit is:\
            \nm = {mean.getValV()}; width = {sigma.getValV()}\nchi2 = {chi}\n")
        os.chdir(os.path.dirname(os. getcwd()))
        f.Close()


def resonance_prop(infile,mu_c=None, particle="all"):
    """The function creates different plots of the main properties of the particles,
        visible as resonances in the plot of the dimuon invariant mass.
    It takes in input:
    - the root data file obtained by the function "leptons_analysis", which contains
        pt, phi, eta and invariant mass of the dimuons, named "dimuon.root";
    - a string with the name of the resonance; the default value is "all", in this
        case plots of all resonances' properties are created.
    """

    eta_lim=3.5
    if mu_c == None:
        t_name = infile.replace(".root", "")
        rdf_i = ROOT.RDataFrame(t_name,infile)
        rdf = rdf_i.Filter("(Dimuon_pt>0)&&(Dimuon_pt<120)").Filter(f"(Dimuon_eta>-{eta_lim})&&(Dimuon_eta<{eta_lim})")
    else:
        rdf = mu_c.Filter("(Dimuon_pt>0)&&(Dimuon_pt<120)").Filter(f"(Dimuon_eta>-{eta_lim})&&(Dimuon_eta<{eta_lim})")

    name=f"#{particle}"
    error=1

    if particle=="eta":
        rdf_cut = rdf.Filter("(Dimuon_mass>=0.52)&&(Dimuon_mass<=0.57)", "eta cut")
    elif particle=="rho" or particle=="omega":
        rdf_cut = rdf.Filter("(Dimuon_mass>=0.72)&&(Dimuon_mass<=0.84)", "ro/omega cut")
    elif particle=="phi":
        rdf_cut = rdf.Filter("(Dimuon_mass>=0.96)&&(Dimuon_mass<=1.07)", "phi cut")
    elif particle=="J-psi":
        name="J/#psi"
        rdf_cut = rdf.Filter("(Dimuon_mass>=2.65)&&(Dimuon_mass<=3.55)", "J/psi cut")
    elif particle=="psi'":
        rdf_cut = rdf.Filter("(Dimuon_mass>=3.55)&&(Dimuon_mass<=3.85)", "psi(2S) cut")
    elif particle=="Y":
        #rdf_cut = rdf.Filter("(Dimuon_mass>=8.5)&&(Dimuon_mass<=11)", "Y cut")
        eta_lim=5
        name=f"{particle}"
    elif particle=="Z":
        rdf_cut = rdf.Filter("(Dimuon_mass>=82)&&(Dimuon_mass<=98)", "Z cut")
        eta_lim=7.5
        name=f"{particle}"
    elif particle=="all":
        for p in ["eta", "rho", "phi", "J-psi", "psi'", "Y", "Z"]:
            resonance_prop(infile,mu_c, p)
    else:
        print("Invalid argument!")
        print("Possible arguments are:\"eta\", \"rho\", \"omega\", \"phi\", \"J-psi\", \"psi'\", \"Y\", \"Z\"")
        inp = input("Insert the right string or leave it empty to exit the program:")
        if inp=="":
            sys.exit(1)
        else:
            print(f"The particle chosen is: \"{inp}\"")
            resonance_prop(infile,mu_c, inp)
            error=0;


    if particle!="all" and error!=0:
        c = ROOT.TCanvas()
        c.UseCurrentStyle()
        c.SetGrid()
        ROOT.gStyle.SetOptStat("e")

        if particle=="Y":
            rdf_cut1 = rdf.Filter("(Dimuon_mass>=9.1)&&(Dimuon_mass<=9.75)", "Y cut1")
            rdf_cut2 = rdf.Filter("(Dimuon_mass>=9.75)&&(Dimuon_mass<=10.2)", "Y cut")
            rdf_cut3 = rdf.Filter("(Dimuon_mass>=10.1)&&(Dimuon_mass<=10.6)", "Y cut3")

            n1 = rdf_cut1.Count().GetValue()
            n2 = rdf_cut2.Count().GetValue()
            n3 = rdf_cut3.Count().GetValue()
            nbin1=math.floor(ROOT.TMath.Sqrt(n1))
            nbin2=math.floor(ROOT.TMath.Sqrt(n2))
            nbin3=math.floor(ROOT.TMath.Sqrt(n3))

            # for h_pt, h_eta,h_phi,rdf_cut,name,nbin in zip(hpt ,heta, hphi, rdfcut,par,Nbin):
            #     h_pt = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Transverse momentum {name}",\
            #         f"Transverse momentum {name};"+"p_{t} [MeV];Events", nbin, 0, 120), "Dimuon_pt")
            #     h_eta = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Pseudorapidity {name}", \
            #         f"Pseudorapidity {name};#eta;Events", nbin, -eta_lim, eta_lim), "Dimuon_eta")
            #     h_phi = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Azimuthal angle {name}", \
            #         f"Azimuthal angle {name};#phi [rad];Events", nbin, -3.5, 3.5), "Dimuon_phi")


            h_pt1 = rdf_cut1.Histo1D(ROOT.RDF.TH1DModel(f"Transverse momentum Y(1S)",\
                f"Transverse momentum Y(1S);"+"p_{t} [MeV];Events", nbin1, 0, 120), "Dimuon_pt") #1000
            h_eta1 = rdf_cut1.Histo1D(ROOT.RDF.TH1DModel(f"Pseudorapidity Y(1S)", \
                f"Pseudorapidity Y(1S);#eta;Events", nbin1, -eta_lim, eta_lim), "Dimuon_eta") #500
            h_phi1 = rdf_cut1.Histo1D(ROOT.RDF.TH1DModel(f"Azimuthal angle Y(1S)", \
                f"Azimuthal angle Y(1S);#phi [rad];Events", nbin1, -3.5, 3.5), "Dimuon_phi") #500

            h_pt2 = rdf_cut2.Histo1D(ROOT.RDF.TH1DModel(f"Transverse momentum Y(2S)",\
                f"Transverse momentum Y(2S);"+"p_{t} [MeV];Events", nbin2, 0, 120), "Dimuon_pt") #1000
            h_eta2 = rdf_cut2.Histo1D(ROOT.RDF.TH1DModel(f"Pseudorapidity Y(2S)", \
                f"Pseudorapidity Y(2S);#eta;Events", nbin2, -eta_lim, eta_lim), "Dimuon_eta") #500
            h_phi2 = rdf_cut2.Histo1D(ROOT.RDF.TH1DModel(f"Azimuthal angle Y(2S)", \
                f"Azimuthal angle Y(2S);#phi [rad];Events", nbin2, -3.5, 3.5), "Dimuon_phi") #500

            h_pt3 = rdf_cut3.Histo1D(ROOT.RDF.TH1DModel(f"Transverse momentum Y(3S)",\
                f"Transverse momentum Y(3S);"+"p_{t} [MeV];Events", nbin3, 0, 120), "Dimuon_pt") #1000
            h_eta3 = rdf_cut3.Histo1D(ROOT.RDF.TH1DModel(f"Pseudorapidity Y(3S)", \
                f"Pseudorapidity Y(3S);#eta;Events", nbin3, -eta_lim, eta_lim), "Dimuon_eta") #500
            h_phi3 = rdf_cut3.Histo1D(ROOT.RDF.TH1DModel(f"Azimuthal angle Y(3S)", \
                f"Azimuthal angle Y(3S);#phi [rad];Events", nbin3, -3.5, 3.5), "Dimuon_phi") #500

            yeta=[h_eta1, h_eta2, h_eta3]
            yphi=[h_phi1, h_phi2, h_phi3]
            ypt=[h_pt1, h_pt2, h_pt3]
            total = yeta+yphi+ypt
            par = ["Y(1S)","Y(2S)","Y(3S)"]

            for h in total:
                h.GetXaxis().CenterTitle()
                h.GetYaxis().CenterTitle()
                h.GetXaxis().SetTitleSize(0.04)
                h.SetFillColorAlpha(9, .8)

            for h_eta, h_phi, h_pt, particle in zip(yeta, yphi, ypt, par):
                h_eta.Draw()
                os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Properties')))
                c.Print(f"{particle}_properties.pdf[", "pdf")
                c.Print(f"{particle}_properties.pdf", f"Title:{particle} pseudorapidity")
                c.SaveAs(f"{particle} pseudorapidity.png")
                h_phi.Draw()
                c.Print(f"{particle}_properties.pdf", f"Title:{particle} azimuthal angle")
                c.SaveAs(f"{particle} azimuthal angles.png")
                h_pt.Draw()
                c.Print(f"{particle}_properties.pdf",f"Title:{particle} transverse momentum")
                c.Print(f"{particle}_properties.pdf]", f"Title:{particle} transverse momentum")
                c.SaveAs(f"{particle} transverse momentum.png")
                os.chdir(os.path.dirname(os. getcwd()))

        else:
            n = rdf_cut.Count().GetValue()
            nbin=math.floor(ROOT.TMath.Sqrt(n))
            h_pt = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Transverse momentum {name}",\
                f"Transverse momentum {name};"+"p_{t} [MeV];Events", nbin, 0, 120), "Dimuon_pt") #1000
            h_eta = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Pseudorapidity {name}", \
                f"Pseudorapidity {name};#eta;Events", nbin, -eta_lim, eta_lim), "Dimuon_eta") #500
            h_phi = rdf_cut.Histo1D(ROOT.RDF.TH1DModel(f"Azimuthal angle {name}", \
                f"Azimuthal angle {name};#phi [rad];Events", nbin, -3.5, 3.5), "Dimuon_phi") #500

            for h in [h_eta, h_pt, h_phi]:
                h.GetXaxis().CenterTitle()
                h.GetYaxis().CenterTitle()
                h.GetXaxis().SetTitleSize(0.04)
                h.SetFillColorAlpha(9, .8)

            h_eta.Draw()
            os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Properties')))
            c.Print(f"{particle}_properties.pdf[", "pdf")
            c.Print(f"{particle}_properties.pdf", f"Title:{particle} pseudorapidity")
            c.SaveAs(f"{particle} pseudorapidity.png")
            h_phi.Draw()
            c.Print(f"{particle}_properties.pdf", f"Title:{particle} azimuthal angle")
            c.SaveAs(f"{particle} azimuthal angles.png")
            h_pt.Draw()
            c.Print(f"{particle}_properties.pdf",f"Title:{particle} transverse momentum")
            c.Print(f"{particle}_properties.pdf]", f"Title:{particle} transverse momentum")
            c.SaveAs(f"{particle} transverse momentum.png")
            os.chdir(os.path.dirname(os. getcwd()))


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
    os.makedirs("Spectra", exist_ok=True)
    logger.debug("The new directory \"Spectra\" is created")
    #mumu_spectrum("dimuon.root", mu_cache)

    #BUMP
    #mumu_spectrum_bump("dimuon.root", mu_cache)

    #ETA DISTRIBUTION
    #mumu_eta("dimuon.root", mu_cache)

    #FIT
    os.makedirs("Fit", exist_ok=True)
    logger.debug("The new directory \"Fit\" is created")
    #resonance_fit("dimuon.root", f"{args.particle}")

    #PROPERTIES
    os.makedirs("Properties", exist_ok=True)
    logger.debug("The new directory \"Properties\" is created")
    resonance_prop("dimuon.root",mu_cache,f"{args.particle}")

    #Z ANALYSIS
    os.makedirs("Z analysis", exist_ok=True)
    logger.debug("The new directory \"Z analysis\" is created")


    #print elapsed time
    timer.Stop()
    logger.info(f"Elapsed time from the beginning is: {timer.RealTime()}(real time) seconds")
