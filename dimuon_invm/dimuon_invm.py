"""
Analysis of the dimuon invariant mass spectrum.
The script takes as arguments:

    - (-f) the data file (URL) of dileptons (they have to be NanoAOD type), for
        example:"root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/
        Run2012B_DoubleMuParked.root";
    - (-p) the string of the particle's name to analyze and fit, among those
        in the dimuon spectrum, which are :
        "eta", "rho","omega", "phi", "J-psi", "psi'", "Y", "Z".
        Better put the string in quotes, because for the "psi'" the " ' "
        character gives some troubles. Default value is \"all\", so it
        returns fit and properties of all resonances.

Furthermore, there are other optional arguments in order to skip the analysis:

    - (-no--a) if we have already done the analysis to create a snapshot of the
        useful data (through \"leptons_analysis\"), we can then skip this step
        passing this string as argument. But it's also necessary to pass the
        the name of the data file, through the argument (-o).
    - (-o) the name of the output root file created from a previous analysis.

In the analysis the following version have been used:

    - Python v3.8
    - ROOT v6.24 ("source ~/root/bin/thisroot.sh" command needed before starting
        the analysis to set the environment of ROOT)

"""


import argparse
import logging
import os
import math
import sys
import time
import ROOT

# Add my modules to the path
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT_DIR)
sys.path.insert(0, os.path.abspath('../Utils'))
import utils

# pylint: disable=E1101
# (9.17/10)


def leptons_analysis(url, outfile):
    """
    It takes in input a nano-AOD data file from which it selects interesting
    muon pairs for further analysis. It also creates a root file, named
    "dimuon_*run*.root" with four new useful quantities:
    ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi"].
    It returns the RDataFrame cached in memory with the columns listed above.

    :param url: url of the root file to download from the web.
    :type url: url of root file, required.
    :param outfile: name of root file in output
    :type outfile: string, required
    :return: data cached in memory
    :rtype: RDataFrame cached.

    """

    # Create an RDataFrame from the root file.
    root_file = ROOT.TFile.Open(url, "READ")
    tree_in = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_in, root_file)
    logger.debug("The RDataFrame has been created.\n")

    # Filter two muons with opposite charge
    rdf_mu = rdf.Filter("nMuon>1 && Muon_charge[0]!=Muon_charge[1]",
                "Selection of two muons with opposite charge").\
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
                 Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[2]")

    # Print cutflows
    logger.info("CutFlow muons:")
    rdf_mu.Report().Print()

    # List with the names of the new columns
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi"]

    # A snapshot is done to collect data in a single root file and they are also
    # cached in memory. A node Graph is saved.
    ROOT.RDF.SaveGraph(rdf_mu, "rdf_mu.dot")
    tree_out = outfile.replace(".root", "")
    logger.debug("Starting snapshots...")
    dimu_cache = rdf_mu.Snapshot(tree_out, outfile, branchlist_mu).Cache()
    logger.info("The snapshot is done and data are cached.")

    # Close the file and return data cached.
    root_file.Close()
    return dimu_cache


def mumu_spectrum(infile, mu_cached=None):
    """
    It takes in input the root data file or data cached, obtained by the function
     "leptons_analysis" and plot the histogram of the dimuons invariant mass,
     named \"dimuon_spectrum.png\".

    :param infile: name of data file
    :type infile: string
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, not required, default=None

    """

    # Histogram of dimuons mass
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        h_mu = rdf.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")
    else:
        h_mu = mu_cached.Histo1D(ROOT.RDF.TH1DModel("Dimuon mass", "Dimuon mass",
            50000, 0.3, 200), "Dimuon_mass")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    c_mu = ROOT.TCanvas("dimuon spectrum", "#mu^{+}#mu^{-} invariant mass")
    c_mu.SetLogx()
    c_mu.SetLogy()
    h_mu.GetXaxis().SetTitle("m_{#mu^{+}#mu^{-}} [GeV]")
    h_mu.GetXaxis().SetTitleSize(0.04)
    h_mu.GetXaxis().CenterTitle()
    h_mu.GetYaxis().SetTitle("Events")
    h_mu.GetYaxis().SetTitleSize(0.04)
    h_mu.GetYaxis().CenterTitle()
    h_mu.SetLineColor(601)
    h_mu.Draw()

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

    # Save results in png in the folder "Spectrum"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectrum')))
    c_mu.SaveAs("dimuon_spectrum.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The plot \"dimuon_spectrum.png\" has been created.")


def mumu_eta(infile, mu_cached=None):
    """
    It takes in input the root data file or data cached obtained by the function
    "leptons_analysis" and plot the histogram of the dimuons pseudorapidity,
    named \"dimuon_eta.png\".

    :param infile: name of data file
    :type infile: string
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, not required, default=None

    """

    # Histogram of dimuons pseudorapidity
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        h_eta = rdf.Filter("Dimuon_pt<120").Histo1D(ROOT.RDF.TH1DModel(
            "Dimuon eta", "Dimuon eta", 160, -7, 7), "Dimuon_eta")
    else:
        h_eta = mu_cached.Filter("Dimuon_pt<120").Histo1D(ROOT.RDF.TH1DModel(
            "Dimuon eta", "Dimuon eta", 160, -7, 7), "Dimuon_eta")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    c_eta = ROOT.TCanvas("dimuon spectrum", "#mu^{+}#mu^{-} pt")
    c_eta.SetGrid()
    h_eta.GetXaxis().SetTitle("#eta_{#mu^{+}#mu^{-}} [GeV]")
    h_eta.GetXaxis().SetTitleSize(0.04)
    h_eta.GetXaxis().CenterTitle()
    h_eta.GetYaxis().SetTitle("Events")
    h_eta.GetYaxis().SetTitleSize(0.04)
    h_eta.GetYaxis().SetTitleOffset(1.3)
    h_eta.GetYaxis().CenterTitle()
    h_eta.SetFillColorAlpha(601, .8)
    h_eta.Draw()

    # Save results in png in the folder "Spectrum"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Spectrum')))
    c_eta.SaveAs("dimuon_eta.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("The plot \"dimuon_eta.png\" has been created.")


def resonance_fit(infile, particle="all"):
    """
    It takes in input the root data file obtained by "leptons_analysis" and
    a string with the name of the resonance to fit (default value is \"all\",
    so in this case all resonances are fitted). Possible arguments are:
    \"eta\",\"rho\",\"omega\",\"phi\",\"J-psi\",\"psi'\",\"Y\",\"Z\",\"all\"."
    It creates a plot with the fitted data and a txt file with fit results of
    the chosen resonance(s).
    If the string passed as particle to fit is an invalid argument, the function
    raise a SyntaxError.

    :param infile: name of data file to analyze
    :type infile: string, required
    :param particle: name of the particle to fit
    :type particle: string, default="all"

    """

    # Start the timer
    start_fit = time.time()

    # An auxiliary TTree from the root data file is created.
    f = ROOT.TFile.Open(infile)
    tree_name = infile.replace(".root", "")
    tree = f.Get(f"{tree_name}")
    logger.debug("The tree is created.")

    low_edge_pt = 20
    upper_edge_pt = 120
    eta_edge = 1.2

    # All background is characterized by Chebychev polynomials
    bkg = "Chebychev"
    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    # Check on input arguments
    if particle=="all":
        for p in arguments:
            resonance_fit(infile, p)
    elif particle not in arguments:
        raise SyntaxError
    elif particle!="all":
        logger.info(f"Fit of particle {particle}...")

        # Mass limits are defined in order to make a selection on data for each
        # particle. The shapes of signal and background are choosed, too.
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

        c_fit = ROOT.TCanvas(f"{name} Resonance", f"{name} Resonance")

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
        c_fit.SaveAs(f"{particle}_fit.png")
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
    If the string passed as particle is an invalid argument, the function raise
    a SyntaxError.

    :param infile: name of data root file
    :type infile: string, required
    :param mu_cached: data cached in memory
    :type mu_cached: RDataFrame, not required, default=None
    :param particle: name of the particle to fit
    :type particle: string, default="all"

    """

    # Cuts on pt
    pt_cut_inf = 0
    pt_cut_sup = 120

    # Rdataframe with useful cuts is created.
    if mu_cached is None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_cut=rdf.Filter(f" Dimuon_pt>{pt_cut_inf} && Dimuon_pt<{pt_cut_sup}")
    else:
        rdf_cut=mu_cached.Filter(f" Dimuon_pt>{pt_cut_inf} && Dimuon_pt<{pt_cut_sup}")

    arguments = ["eta", "rho", "omega", "phi", "J-psi", "psi'", "Y", "Z"]

    # Check on input arguments
    if particle=="all":
        for p in arguments:
            resonance_prop(infile, mu_cached,p)
    elif particle not in arguments:
        raise SyntaxError
    elif particle!="all":
        logger.info(f"Plot properties of particle {particle}...")

        # Retrieve values of mass range for each particle
        lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][0]
        upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"{particle}"][1]
        rdf_m = rdf_cut.Filter(f" Dimuon_mass>={lower_mass_edge} \
            && Dimuon_mass<={upper_mass_edge}", f"{particle} cut")

        # Styling
        name = utils.FIT_INIT_PARAM[f"{particle}"][6]
        c_prop = ROOT.TCanvas()
        c_prop.UseCurrentStyle()
        c_prop.SetGrid()
        ROOT.gStyle.SetOptStat("e")


        # Properties for the three Y resonances
        if particle=="Y":
            for i in range(1, 4, 1):
                # Retrieve values of mass range and do the cuts
                lower_mass_edge = utils.PARTICLES_MASS_RANGE[f"Y{i}"][0]
                upper_mass_edge = utils.PARTICLES_MASS_RANGE[f"Y{i}"][1]
                rdf_m_y = rdf_m.Filter(f"(Dimuon_mass>={lower_mass_edge})\
                    &&(Dimuon_mass<={upper_mass_edge})", f"Y cut{i}S")

                n = rdf_m_y.Count().GetValue()
                nbin_y = math.floor(ROOT.TMath.Sqrt(n)/2)

                # Create a dictonary with all histograms booked
                h_y = {}
                h_y["Transverse momentum"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Transverse momentum Y({i}S)",f"Transverse momentum Y({i}S);"
                    "p_{t} [MeV];Events", nbin_y, 0, 120), "Dimuon_pt")
                h_y["Pseudorapidity"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Pseudorapidity Y({i}S)", f"Pseudorapidity Y({i}S);#eta;"
                    "Events", nbin_y, -5.5, 5.5), "Dimuon_eta")
                h_y["Azimuthal angle"]=rdf_m_y.Histo1D(ROOT.RDF.TH1DModel(
                    f"Azimuthal angle Y({i}S)",f"Azimuthal angle Y({i}S);#phi "
                    "[rad];Events", nbin_y, -3.5, 3.5), "Dimuon_phi")

            # Do and save the histogram in the directory "Properties"
                for j, (k, h) in enumerate(h_y.items()):
                    h.GetXaxis().CenterTitle()
                    h.GetYaxis().CenterTitle()
                    h.GetXaxis().SetTitleSize(0.04)
                    h.SetFillColorAlpha(601, .9)
                    if j==0:
                        os.chdir(os.path.abspath(os.path.join(os.sep,
                            f'{os.getcwd()}', 'Properties')))
                        c_prop.Print(f"Y({i}S) properties.pdf[", "pdf")
                    h.Draw()
                    c_prop.Print(f"Y({i}S) properties.pdf", f"Title:Y({i}S) {k}")
                    c_prop.SaveAs(f"Y({i}S) {k}.png")
                c_prop.Print(f"Y({i}S) properties.pdf]", "pdf")

                # Return in main directory
                os.chdir(os.path.dirname(os. getcwd()))
        else:
            # Plots for the other particles
            n = rdf_m.Count().GetValue()
            nbin = math.floor(ROOT.TMath.Sqrt(n)/2)

            # Range of rapidity
            eta_max = rdf_m.Max("Dimuon_eta").GetValue()
            eta_lim = eta_max+(0.02*eta_max)

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
                h.SetFillColorAlpha(601, .9)
                if j==0:
                    os.chdir(os.path.abspath(os.path.join(os.sep,
                        f'{os.getcwd()}','Properties')))
                    c_prop.Print(f"{particle} properties.pdf[", "pdf")
                h.Draw()
                c_prop.Print(f"{particle} properties.pdf", f"Title:{particle} {k}")
                c_prop.SaveAs(f"{particle} {k}.png")
            c_prop.Print(f"{particle} properties.pdf]", "pdf")

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
    N_THREADS = 4                                                                # (60-70% CPU)
    ROOT.ROOT.EnableImplicitMT(N_THREADS)

    # Start the the timer
    start = time.time()

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file", type=str,
        help="The path of the nanoAOD file to analyze.")
    parser.add_argument("-p", "--particle",required=False, type=str, nargs="*",
        default=["all"],help="The name of the particle to analyze. It is necessary"
            " if you want to see resonances' fit and properties. Possible "
            "arguments are:\"eta\",\"rho\",\"omega\",\"phi\",\"J-psi\", \"psi'\","
            " \"Y\",\"Z\", \"all\". Put the chosen string in quotes.")
    parser.add_argument("-no--a", "--analysis", required=False, action="store_false",
        default=True, help="If True, retrieve the data file from web and "
        "does the selection. Otherwise it uses the \"dimuon_*run*.root\" file, "
        "creted in a previous run of the script. Of course in this case, it's "
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
    # If "leptons_analysis" is not run, dimu_cached is set to None, so in the
    # following the other functions take data directly from the root file. But
    # the analysis must be run at least one to collect data.

    if args.analysis is True:
        try:
            s=args.file.rfind("/")
        except AttributeError as ex:
            logger.error("If you want to run the analaysis you have to insert "
                f"the path of the data files! \n {ex}")
            sys.exit(1)
        else:
            outfile_m = f"dimuon_{args.file[s+1:]}"
            dimu_cached = leptons_analysis(f"{args.file}", outfile_m)
    else:
        if not os.path.isfile(args.outfile):
            raise IOError(f"The file {args.outfile} doesn't exist. Maybe you have"
                           " to run the analysis first.")
        dimu_cached = None
        outfile_m = args.outfile

    # DIMUON MASS SPECTRUM
    os.makedirs("Spectrum", exist_ok=True)
    logger.debug("The new directory \"Spectrum\" is created")
    mumu_spectrum(outfile_m, dimu_cached)

    # DIMUON ETA DISTRIBUTION
    mumu_eta(outfile_m, dimu_cached)

    # RESONANCE'S FIT & PROPERTIES
    os.makedirs("Fit", exist_ok=True)
    os.makedirs("Properties", exist_ok=True)
    logger.debug("The new directories \"Fit\" and \"Properties\" are created")

    for part in args.particle:
        try:
            resonance_fit(outfile_m, part)
            resonance_prop(outfile_m, dimu_cached, part)
        except SyntaxError:
            logger.error(f"Invalid argument: \"{part}\"!\nPossible arguments are:"
                "\"eta\", \"rho\",\"omega\",\"phi\", \"J-psi\", \"psi'\", \"Y\","
                " \"Z\" or \"all\"")

    # Elapsed time
    logger.info(f" Elapsed time from the beginning is:"
        f" {time.time()-start}(real time) seconds")
