"""
The script takes as argument:

    - the data file (URL) of dileptons (-f), for example: "root:
        //eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/
        Run2012B_DoubleMuParked.root";
    - the string of the particle's name (-p) to analyze and fit, among those
        in the dimuon spectrum, which are :
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

In the analysis the following version have been used:

    - Python v3.8
    - ROOT v6.24 ("source ~/root/bin/thisroot.sh" command needed before starting
        the analysis to set the environment of ROOT)

"""

#"root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root"

import argparse
import logging
import os
import math
import sys
import time
import ROOT

# Add my modules to the path
root_utils = os.path.abspath('../Utils')
sys.path.insert(0, root_utils)
import utils


def leptons_analysis(url, outfile):
    """
    It takes in input a nano-AOD data file and creates a root files, named:
    "dimuon.root" with four useful columns for further analysis:
    ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi"].
    It returns the RDataFrame cached in memory with the columns listed above.

    :param url: url of the root file to upload from the web.
    :type url: url of root file, required.
    :param outfile: name of root file in output
    :type outfile: string, required
    :return: data cached in memory
    :rtype: RDataFrame cached.

    """

    # Create an RDataFrame from the root file.
    root_file = ROOT.TFile.Open(url, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug("The RDataFrame have been created.")

    # Filter two muons and two electrons with opposite charge
    rdf_mu =rdf.Filter("nMuon==2","Selection of two muons").\
                Filter("Muon_charge[0]!=Muon_charge[1]",\
                        "Selection of muons with opposite charge").\
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

    # Print cutflows
    logger.info("CutFlow muons:")
    rdf_mu.Report().Print()
    print(f"Nslots: {rdf_mu.GetNSlots()}")

    # Save Node Graph
    ROOT.RDF.SaveGraph(rdf_mu, "rdf_mu.dot")

    # List with the names of the new columns to make a snapshot of them
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi", \
        "Dimuon_cos", "Dimuon_y"]

    # A snapshot is done to collect the useful physiscal quantity of the dimuons
    # in a single root file. Data are also cached in memory.
    logger.debug("Starting snapshots...")
    tree = outfile_m.replace(".root", "")
    dimu_cache = rdf_mu.Snapshot(tree, outfile, branchlist_mu).\
        Cache()
    logger.info("The snapshot is done and data are cached.")

    # Close the file and return data cached.
    root_file.Close()
    return dimu_cache


def mumu_spectrum(infile, mu_cached=None, bump=False):
    """
    It takes in input the root data file or the data cached, obtained by the
    function "leptons_analysis" and plot the histogram of the dimuons invariant
    mass. It is also possible to eliminate the bump between the Y and Z
    resonances by setting the argument "bump=True"; it makes a cut on pt>20 GeV.

    :param infile: name of data root file
    :type infile: string, required
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, NOT REQUIRED, default=False
    :param bump: False to plot without the bump
    :type bump: boolean, NOT REQUIRED

    """

    # Histogram of dimuons mass
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        if bump is False:
            rdf = rdf.Filter("Dimuon_pt>20")
        h = rdf.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")
    else:
        if bump is False:
            mu_cached = mu_cached.Filter("Dimuon_pt>20")
        h = mu_cached.Filter("Dimuon_pt>20").Histo1D(
            ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass", 50000, 0.3, 200),
            "Dimuon_mass")

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

    # Save results in png in the folder "Spectra"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectra')))
    if bump is False:
        c.SaveAs("dimuon_spectrum.png")
    else:
        c.SaveAs("dimuon_spectrum_wo_bump.png")


    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The plot \"dimuon_spectrum.png\" has been created.")


def mumu_eta(infile, mu_cached=None):
    """
    It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" and plot the histogram of the dimuons
    pseudorapidity.

    """

    # Histogram of dimuons pseudorapidity
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        #rdf_m = rdf.Filter("Dimuon_pt>20")
        h = rdf.Filter("Dimuon_pt>20").Histo1D(ROOT.RDF.TH1DModel(
            "Dimuon eta", "Dimuon eta", 150, -4, 4), "Dimuon_eta")
        del rdf

    else:
        #rdf_m = mu_cached.Filter("Dimuon_pt>20")
        h = rdf.Filter("Dimuon_pt>20").Histo1D(ROOT.RDF.TH1DModel(
            "Dimuon eta", "Dimuon eta", 150, -4, 4), "Dimuon_eta")

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
    """
    It takes in input the root data file obtained with "leptons_analysis"
    and a string with the name of the resonance to fit. Possible arguments
    are: \"eta\", \"rho\",\"omega\", \"phi\", \"J-psi\", \"psi'\", \"Y\",
    \"Z\", \"all\"."
    It creates a plot with the fit and a txt file with fit results, for the
    chosen resonance(s).

    :param infile: Data file to analyze
    :type infile: root file, required
    :param particle: name of the particle to fit
    :type particle: string, required
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, NOT REQUIRED, default=False

    """

    # Start the timer
    start_fit = time.time()

    # An auxiliary TTree from the root data file is created.
    if mu_cached is None:
        f = ROOT.TFile.Open(infile)
        tree_name = infile.replace(".root", "")
        tree = f.Get(f"{tree_name}")
        logger.debug("The tree is created.")
    else:
        logger.debug("Devi aggiornare la versione di root a 6.26 e usare"
            "RooDataSetHelper")

    low_edge_pt = 20
    upper_edge_pt = 100
    eta_edge = 1.2
    # low_edge_eta = -1.2
    # upper_edge_eta = 1.2

    error_flag=1

    # All background is characterized by Chebychev polynomials
    bkg = "Chebychev"
    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    if particle=="all":
        for p in arguments:
            resonance_fit(infile, p)
    elif particle not in arguments:
        logger.error("Invalid argument! \nPossible arguments are:\"eta\", "
            "\"rho\", \"omega\",\"phi\", \"J-psi\", \"psi'\", \"Y\", \"Z\" or "
            "\"all\"")
        inp = input("Insert the right string or leave it empty to exit the "
            "program:")
        if inp=="":
            sys.exit(1)
        else:
            logger.info(f"The input chosen is: \"{inp}\"")
            resonance_fit(infile, inp)
            error_flag=0

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

        # Retrieve values of mass range for each particles
        lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][0]
        upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][1]

        Dimuon_mass = ROOT.RooRealVar("Dimuon_mass", "Dimuon_mass",
            lower_mass_edge,upper_mass_edge)
        Dimuon_pt = ROOT.RooRealVar("Dimuon_pt", "Dimuon_pt", low_edge_pt,
            upper_edge_pt)
        Dimuon_eta = ROOT.RooRealVar("Dimuon_eta", "Dimuon_eta", -eta_edge,
            +eta_edge)

        # Define "signal", "background" and cuts for each particle in order
        # to fit them.
        cut = ROOT.RooFormulaVar("cut on mass, pt, eta",\
            f"Dimuon_mass>{lower_mass_edge} && Dimuon_mass<{upper_mass_edge} \
            && Dimuon_pt>{low_edge_pt} && Dimuon_pt<{upper_edge_pt} \
            && Dimuon_eta>{-eta_edge} && Dimuon_eta<{eta_edge} ", \
            ROOT.RooArgList(Dimuon_mass, Dimuon_pt, Dimuon_eta))
        rds = ROOT.RooDataSet("rds","rds",tree,
            ROOT.RooArgSet(Dimuon_mass, Dimuon_pt, Dimuon_eta),cut)
        logger.debug("The roodataset is created from the tree with cut on mass,"
            " pt, and eta.")
        num = rds.sumEntries()

        if sig=="gaus":
            sigf=ROOT.RooGaussian(sig, "The resonance", Dimuon_mass, mean, sigma)
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
            # NEW
            # eff = ROOT.RooRealVar("eff", "The efficiency", 0.20, 0.00001, 1.)
            # lum = ROOT.RooRealVar("lum", "The luminosity", 2, 0.00001, 50, "pb-1")
            # cross = ROOT.RooRealVar("cross", "The cross section", 3., 0., 40, "pb")
            # nsig = ROOT.RooFormulaVar("N", "@0*@1*@2", ROOT.RooArgList(eff, lum, cross))
            # eff.setConstant(1)
            # lum.setConstant(1)

            # OLD
            nsig = ROOT.RooRealVar("nsig", "signal events", num/5, 0, num)
            tot = ROOT.RooAddPdf("Model", "The total pdf",
                ROOT.RooArgList(sigf, bkgf), ROOT.RooArgList(nsig, nbkg))

        # Fit properties
        opt_list = ROOT.RooLinkedList()
        opt_list.Add(ROOT.RooFit.Extended(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.Save())
        opt_list.Add(ROOT.RooFit.BatchMode(ROOT.kTRUE))
        opt_list.Add(ROOT.RooFit.NumCPU(4))
        opt_list.Add(ROOT.RooFit.PrintLevel(-1))

        # Fit
        results = tot.fitTo(rds,opt_list)

        # Retrieve number of free parameters
        ndf = results.floatParsFinal().getSize()

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
        chi2 = f"chi/ndf = {round(chi*ndf,2)}/{ndf} = {round(chi,2)}"

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
            f"the fit are:\n"
            f"m(Y1S) = {mean1.getValV()}; width(Y1S) = {sigma1.getValV()}\n"
            f"m(Y2S) = {mean2.getValV()}; width(Y2S) = {sigma2.getValV()}\n"
            f"m(Y3S) = {mean3.getValV()}; width(Y3S) = {sigma3.getValV()}\n"
            f"chi2\\ndf = {chi}\n")
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
                f"from the fit are: \n"
                f"m = {mean.getValV()}; width = {sigma.getValV()}\n"
                f"chi2\\ndf = {chi}\n")

        leg.AddEntry(f"{bkg}",f"{bkg}", "l")
        leg.AddEntry("Model","Total", "l")
        leg.SetShadowColor(0)
        leg.Draw()

        #Save results in the directory "Fit"
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Fit')))
        c.SaveAs(f"{particle}_fit.png")
        utils.write_fitresults(results, f"res_{particle}.txt")

        # Return in main directory
        os.chdir(os.path.dirname(os. getcwd()))
        f.Close()
        logger.info(f"Elapsed time to fit {particle}: {time.time()-start_fit}")


def resonance_prop(infile, mu_cached=None, particle="all"):
    """
    The function creates different plots of the main properties of the
    resonances such as transverse momentum, pseudorapidity and azimuthal angle.
    It takes in input the root data file or the data cached obtained by the
    function "leptons_analysis" and a string with the name of the resonance;
    the default value is "all", in this case plots of all resonances' properties
    are created.

    :param infile: name of data root file
    :type infile: string, required
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, NOT REQUIRED, default=False
    :param particle: name of the particle to fit
    :type bump: string, default="all"

    """

    # Cuts on pt
    pt_cut_inf = 0
    pt_cut_sup = 120

    error_flag = 1

    # Rdataframe with useful cuts is created.
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_cut=rdf.Filter(f" Dimuon_pt>{pt_cut_inf} && Dimuon_pt<{pt_cut_sup}")
        del rdf
    else:
        rdf_cut=mu_cached.Filter(f" Dimuon_pt>{pt_cut_inf} && Dimuon_pt<{pt_cut_sup}")

    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    # Check on input arguments
    if particle=="all":
        for p in arguments:
            resonance_prop(infile,mu_cached, p)
    elif particle not in arguments:
        logger.error("Invalid argument! \nPossible arguments are:\"eta\", "
            "\"rho\", \"omega\",\"phi\", \"J-psi\", \"psi'\", \"Y\", \"Z\" or "
            "\"all\"")
        inp = input("Insert the right string or leave it empty to exit the "
            "program:")
        if inp=="":
            sys.exit(1)
        else:
            logger.info(f"The particle chosen is: \"{inp}\"")
            resonance_prop(infile,mu_cached, inp)
            error_flag=0

    if particle!="all" and error_flag:

        logger.info(f"Plot properties of particle {particle}...")

        # Retrieve values of mass range for each particle
        lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][0]
        upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][1]
        rdf_m = rdf_cut.Filter(f" Dimuon_mass>={lower_mass_edge} \
            && Dimuon_mass<={upper_mass_edge}", f"{particle} cut")
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
                h.Draw()
                c.Print(f"{particle} properties.pdf", f"Title:{particle} {k}")
                c.SaveAs(f"{particle} {k}.png")
            c.Print(f"{particle} properties.pdf]", "pdf")

            # Return in main directory
            os.chdir(os.path.dirname(os. getcwd()))




if __name__ == "__main__":

    # Set ROOT environment
    ROOT.gROOT.SetBatch()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)

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

    # Enable parallel analysis
    nthreads = 8
    ROOT.ROOT.EnableImplicitMT(nthreads)

    # Start the the timer
    start = time.time()
    start_cpu = time.process_time()
    print(start, start_cpu)

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file", type=str,
        help="The path of the nanoAOD file to analyze.")
    parser.add_argument("-p", "--particle",required=False, type=str, nargs="+",
        help="The name of the particle to analyze. It is necessary if you want \
            to see resonances' fit and properties. Possible arguments are:\
            \"eta\", \"rho\",\"omega\", \"phi\", \"J-psi\", \"psi'\", \"Y\", \
            \"Z\", \"all\". Put the chosen string in quotes.")
    parser.add_argument("-no--a", "--analysis", required=False, action="store_false",
        default=True, help="If True, retrieve the data file from web and "
        "does the selection. Otherwise it uses the \"dimuon.root\" file, "
        "created in a previous run of the script. Of corse in this case, it's "
        "needed to create the snapshot of the data before." )
    parser.add_argument("-o", "--outfile", required=False, type=str,
        help="if \"-no--aa\" is called to avoid the analysis, it's necessary"
        " to pass the name of the data file by this argument.")
    args = parser.parse_args()

    # Creating the logger and setting its level
    logger = utils.set_logger("Dimuon_invm", logging.DEBUG)

    logger.info("Starting the analysis of the root file nanoAOD...")

    # Load the shared library "tools.cpp" which contains some functions to
    # calculate the useful quantities for the analysis.
    ROOT.gSystem.Load('../Utils/tools_cpp.so')

    # LEPTONS ANALYSIS

    # If "leptons_analysis" is not run (it already has been done), dimu_cached
    # is set to None, so in the following the other functions take data directly
    #  from the root file. But "leptons_analysis" must be run at least one to
    # collect data.

    if args.analysis is True:
        try:
            s=args.file.rfind("/")
        except AttributeError as ex:
            logger.error("If you want to run the analaysis you have to insert "
                f"the path of the data files! \n {ex}")
            sys.exit(1)
        else:
            outfile_m = args.file[s+1:]
            dimu_cached = leptons_analysis(f"{args.file}", outfile_m)
    else:
        if not os.path.isfile(args.outfile):
            raise IOError(f"The file {outfile_m} doesn't exist. Maybe you have "
                           "to run the analysis first.")
        dimu_cached = None
        outfile_m = args.outfile

    # DIMUON MASS SPECTRUM
    os.makedirs("Spectra", exist_ok=True)
    logger.debug("The new directory \"Spectra\" is created")
    mumu_spectrum(outfile_m, dimu_cached)


    # # DIMUON ETA DISTRIBUTION
    # ACCETTANZA DELL'ESPERIMENTO DALA DISTRIBUZIONE IN ETA???
    # DIMINUISCE CON IL PT????
    mumu_eta(outfile_m, dimu_cached)

    # RESONANCE'S FIT
    os.makedirs("Fit", exist_ok=True)
    logger.debug("The new directory \"Fit\" is created")
    for p in args.particle:
        resonance_fit(outfile_m, p)

    # # PROPERTIES
    os.makedirs("Properties", exist_ok=True)
    logger.debug("The new directory \"Properties\" is created")
    for p in args.particle:
        resonance_prop(outfile_m, dimu_cached, p)

    # Elapsed time
    logger.info(f" Elapsed time from the beginning is:"
        f" {time.time()-start}(real time) seconds")
    logger.info(f" Elapsed time from the beginning is: "
        f"{time.process_time()-start_cpu}(cpu time) seconds")
