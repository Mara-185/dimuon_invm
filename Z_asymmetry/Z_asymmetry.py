"""
In this module a few functions are implemented in order to study the
\"Forward-Backward asymmetry" for the Z boson decay.
The analysis is done in six different rapidity ranges of equal size and
twelve mass bins:

    - Mass bins : 60 < Mass < 120:

    [60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120];

    - Rapidity bins of equal size in \|rapidity\| < 2.4:

    [0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]

The script takes as arguments:

    - (-f) the txt files (they must finish with \"*_index.txt\") with the url
        of the root data files to analyze. They have to be NanoAOD type.
    - (-t) the type of dataset to analyze. [E.g. "MC" or "data"]. This string
        has to be added for each data file insert by \"-f\".
    - (-c) this argument can be passed to make the comparison plots of mass and
        cos(theta*) distributions and also a plot comparison of the Afb values.
        It's necessary to have analyzed Monte Carlo and collected data previously.

(E.g.: python3 Z_asymmetry.py -f MC_index.txt Run2012B_SingleMu_merged_index.txt
\\ Run2012C_SingleMu_merged_index.txt -t MC data data -c)

Furthermore, there are other optional arguments in order to skip the analysis:

    - (-no--a) if we have already done the analysis to create a snapshot of the
        useful data (through \"z_main\"), we can then skip this step passing
        this string as argument. In this case the name of the file data to
        analyze is read from the file "Analyzed_files.txt" (when the script is
        run w/o \"-no--a\", the names of the snapshots created are wrote down
        in it), and it takes the root files from the directory named "Snapshots".

(E.g.: python3 Z_asymmetry.py -no--a -c)

In the analysis the following versions have been used:

    - Python v3.8
    - ROOT v6.24 ("source ~/root/bin/thisroot.sh" command needed before starting
        the analysis to set the environment of ROOT)

"""

import logging
import os
import argparse
import time
import sys
import ROOT
import numpy as np

# Add my modules to the path
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT_DIR)
sys.path.insert(0, os.path.abspath('../Utils'))
import utils

# pylint: disable=E1101
# (9.39)

# Array with mass and eta bin values
MASS_BIN = np.array([60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120],
    dtype=float)
RAPIDITY_BIN=np.array([0, .4, .8, 1.2, 1.6, 2.0, 2.4], dtype=float)

# Functions

def retrieve_dataset(findex, type_d):
    """
    The function takes as argument the file index with all name of root files
    of the chosen dataset to analyze and its type (e.g. "MC" or "data").
    (The file index has to be downloaded or created manually.)
    Three lists and a standard vector are returned in order to save the time
    needed to analyze each file, its number of events, the name of the
    snapshot created in the "z_main" function and its type.

    :param findex: index of root files which make up the entire dataset
    :type findex: string of txt file, required
    :param type_d: type of file (e.g. "MC" or "data")
    :type type_d: string, required
    :return files: list of the analyzed files
    :rtype files: standard vector of string
    :return timef: list of time elpased analyzing each file
    :rtype timef: list
    :return Nf: list of events analyzed for each file
    :rtype Nf: list
    :return typef: list of data types
    :rtype typef: list

    """

    with open(findex, "r", encoding="utf-8") as inf:
        # Strings of all files in the dataset are loaded.
        data = [_.rstrip('\n') for _ in inf]
        run = findex.replace("_index.txt", "")
        logger.info(f" Total number of data files in {run}: {len(data)}\n")

        timef = []
        Nf = []
        typef=[]
        files = ROOT.std.vector("string")()

        for i,d in enumerate(data):
            logger.info(f" Analyzing file number {i+1} :\n{d}")

            start_analysis = time.time()
            string_file, N_events = z_main(d, i, run)
            stop_analysis = time.time() - start_analysis
            logger.info(f" Elapsed time : {stop_analysis}\n")

            typef.append(type_d)
            files.emplace_back(string_file)
            Nf.append(N_events)
            timef.append(stop_analysis)

    return files, timef, Nf, typef


def z_main(url, iteration, run):
    """
    The function create an RDataFrame to make the selection of the good events
    for the analysis, to compute some useful quantities and store them in new
    columns. It returns the string of the root file created by a snapshot only
    with the columns needed.

    :param url: url of the root file to upload from the web.
    :type url: url of root file, required.
    :param iteration: number of file analyzed
    :type iteration: int, required
    :param run: name of the data run
    :type run: string, required
    :return snap_name: name of the root file created
    :rtype snap_name: string
    :return N_ev: number of all events in the file
    :rtype N_ev: float

    """

    # Create an RDataFrame to make selection on data from the root file.
    root_file = ROOT.TFile.Open(url, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug(" The RDataFrame have been created.")

    # Count of all data in the file.
    N_ev = rdf.Count().GetValue()

    rdf_cut=rdf.Filter("nMuon>1","Selection on events with at least two muons").\
        Define("mask", "abs(Muon_eta)<2.4 && Muon_isGlobal==1 && Muon_gChi2<10 "
        "&& Muon_pfRelIso03_all<0.1 && abs(Muon_dxy)<0.2 && HLT_IsoMu24_eta2p1==1").\
        Define("Mu_pt","Muon_pt[mask]").Define("Mu_eta","Muon_eta[mask]").\
        Define("Mu_mass","Muon_mass[mask]").Define("Mu_phi","Muon_phi[mask]").\
        Define("Mu_charge", "Muon_charge[mask]").\
        Filter("Mu_pt[0]>25 && Mu_pt[1]>15","Trigger on leading and next leading"
            " muon pt").\
        Filter("Mu_charge[0]!=Mu_charge[1]","Selection on opposite charge").\
        Define("Dimuon_mass","dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[3]").\
        Define("Dimuon_pt","dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[0]").\
        Define("Dimuon_eta","dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]").\
        Define("Dimuon_phi","dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[2]").\
        Define("Dimuon_cos","cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_charge[0],Mu_pt[1], Mu_eta[1], Mu_phi[1], \
            Mu_mass[1])[0]").\
        Define("Dimuon_y","cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], \
            Mu_mass[0], Mu_charge[0], Mu_pt[1], Mu_eta[1], Mu_phi[1], \
            Mu_mass[1])[1]").\
        Filter("Dimuon_mass>60 && Dimuon_mass<120","Cut on Z resonance").\
        Define("w_d", "weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[0]").\
        Define("w_n", "weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[1]")

    logger.info("\nReport of all cuts:\n")
    rdf_cut.Report().Print()

    branchlist = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi",\
        "Dimuon_cos","Dimuon_y","w_d", "w_n"]
    snap_name = f"dimuon_w{iteration}[{run}].root"

    # A snapshot is done to collect the useful physiscal quantity of the
    # dileptons in a root file and a node graph is saved.
    ROOT.RDF.SaveGraph(rdf_cut, "rdf_zmain.dot")
    rdf_cut.Snapshot("dimuon_w", snap_name, branchlist)
    logger.info(f" A snapshot of the good events named {snap_name} is done. ")

    root_file.Close()

    return snap_name, N_ev


def cos_eta(rdf_all_c, data_type, pt_lim=(0,120), infile=None):
    """
    Plot the cos of Z resonance (mass range from 60 GeV to 120 GeV), in the six
    different range of rapidity. It also make a cut on pt.
    Plot is saved as "Dimuon_cos[y(rap_inf,rap_sup),pt(pt_inf,pt_sup)]_*dtype*.png"
    in a directory named "Plot".
    It also creates a root file containg all histograms, named \"hcos_*dtype*.root\",
    useful for further comparison between Monte Carlo and data.
    Data can be passed as RDataFrame or as root file (in the last case the name
    of the TTree stored has to be the same of the file).

    :param rdf_all_c: RDataFrame to analyze.
    :type rdf_all_c: RDataFrame, required
    :param data_type: type of dataset analyzed
    :type data_type: string, required
    :param pt_lim: range limits of transverse momentum
    :type pt_lim: tuple, default=(0,120)
    :param infile: data file
    :type infile: root file, not required

    """

    # RDataFrame with appropriate cuts is created.
    if infile is None:
        rdf_c = rdf_all_c.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                f"Cut on pt {pt_lim}").Define("Dimuon_y_abs", "abs(Dimuon_y)")
    else:
        t_name = infile.replace(".root", "")
        rdf_t = ROOT.RDataFrame(t_name,infile)
        rdf_c = rdf_t.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                f"Cut on pt {pt_lim}").Define("Dimuon_y_abs", "abs(Dimuon_y)")

    rdf_c.Report().Print()
    logger.debug(" The RDataFrame is created.")

    # Booking 2D histogram of cos(theta*) in 12 rapidity bin (from -2.4 to 2.4)
    h_all_c=rdf_c.Histo2D(ROOT.RDF.TH2DModel("Cos(theta*) in defined range rapidity",
        "Cos(theta*) vs. rapidity",6, RAPIDITY_BIN, 40, -1, 1),"Dimuon_y_abs",
        "Dimuon_cos")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    ROOT.gStyle.SetTitleBorderSize(2)
    ROOT.gStyle.SetTitleFillColor(19)
    ROOT.gStyle.SetTitleSize(0.0, "t")
    ROOT.gStyle.SetTitleOffset(0.0, "t")
    c_cos = ROOT.TCanvas("cos", "cos")

    # Change directory to save results in "Plot"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Plot')))
    file_cos = ROOT.TFile(f"hcos_{data_type}.root", "recreate")
    # Drawing histograms for each rapidity bin
    for i in range(0, 6, 1):
        ROOT.gStyle.SetOptStat("e")
        h_cos = h_all_c.ProjectionY(f"hc{i}", i+1, i+1)
        h_cos.SetTitle(f"Cos(#theta*) in y[{RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}]"
            ";Cos(#theta*);Events")
        h_cos.SetFillColorAlpha(600, .7)
        h_cos.GetYaxis().SetMaxDigits(3)
        h_cos.Draw()
        h_cos.Write()
        c_cos.SaveAs(f"Dimuon_cos[y({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})"
            f",pt({pt_lim[0]},{pt_lim[1]})]_{data_type}.png")
        logger.info(f" Drawing the cos(theta*) distribution for "
            f"{RAPIDITY_BIN[i]}<|Dimuon rapidity|<{RAPIDITY_BIN[i+1]}")

    file_cos.Close()

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def mass_eta(rdf_all_m, data_type, pt_lim = (0,120), infile=None):
    """
    Plot the Z resonance (mass range from 60 GeV to 120 GeV), in the six
    different range of rapidity. It also make a cut on pt.
    Plot is saved as "Dimuon_mass[y(rap_inf,rap_sup),pt(pt_inf,pt_sup)]_*dtype.png"
    in a directory named "Plot".
    It also creates a root file containg all histograms, named \"hmass_*dtype.root\",
    useful for further comparison between Monte Carlo and data.
    Data can be passed as RDataFrame or as root file (in the last case the name
    of the TTree stored has to be the same of the file).

    :param rdf_all_m: RDataFrame to analyze.
    :type rdf_all_m: RDataFrame, required
    :param data_type: type of dataset analyzied
    :type data_type: string, required
    :param pt_lim: range limits of transverse momentum
    :type pt_lim: tuple, default=(0,120)
    :param infile: data file
    :type infile: root file, not required

    """

    # RDataFrame with appropriate cuts is created.
    if infile is None:
        rdf_m=rdf_all_m.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                f"Cut on pt {pt_lim}").Define("Dimuon_y_abs", "abs(Dimuon_y)")
    else:
        t_name = infile.replace(".root", "")
        rdf_t = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf_t.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                f"Cut on pt {pt_lim}").Define("Dimuon_y_abs", "abs(Dimuon_y)")

    rdf_m.Report().Print()
    logger.debug(" The RDataFrame is created.")

    # Booking 2D histogramof the mass in 12 rapidity bin (from -2.4 to 2.4)
    h_all_m = rdf_m.Histo2D(ROOT.RDF.TH2DModel("Mass in defined range rapidity",
        "Mass vs. rapidity",6, RAPIDITY_BIN, 60, 60, 120), "Dimuon_y_abs",
        "Dimuon_mass")

    # Styling
    ROOT.gStyle.SetOptStat("e")
    ROOT.gStyle.SetTitleBorderSize(2)
    ROOT.gStyle.SetTitleFillColor(19)
    ROOT.gStyle.SetTitleSize(0.0, "t")
    ROOT.gStyle.SetTitleOffset(0.0, "t")
    c_mass = ROOT.TCanvas("mass", "mass")

    # Change directory to save results in "Plot"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Plot')))
    file_mass = ROOT.TFile(f"hmass_{data_type}.root", "recreate")
    # Drawing histograms for each rapidity bin
    for i in range(0, 6, 1):
        h_mass = h_all_m.ProjectionY(f"hm{i}", i+1, i+1)
        h_mass.SetTitle(f"Mass in y[{RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}];"
            "Mass [GeV];Events")
        h_mass.SetFillColorAlpha(601, .7)
        ROOT.gPad.SetLogy()
        h_mass.Draw()
        h_mass.Write()
        c_mass.SaveAs(f"Dimuon_mass[y({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})"
            f",pt({pt_lim[0]},{pt_lim[1]})]_{data_type}.png")
        logger.info(f" Drawing the mass distribution for {RAPIDITY_BIN[i]}<"
                    f"|Dimuon rapidity|<{RAPIDITY_BIN[i+1]}")

    file_mass.Close()

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def afb(rdf_all, data_type, pt_lim=(0,120)):
    """
    The function calculates the value of Afb (Forward-backward Asymmetry)
    using the \"angular event weighting\" in the six different ranges of
    rapidity for each bin mass.
    It creates six different plots (one for each rapidity range)
    (named \"Afb[y(rap_inf,rap_sup),pt(pt_inf,pt_sup)]_*dype.png\") and a plot
    with all of them (named \"afb_y_*dtype.png\").

    :param rdf_all: RDataframe to compute Afb
    :type rdf_all: RDataFrame, required
    :param data_type: type of dataset analyzed
    :type data_type: string, required
    :param pt_lim: range limits of pt
    :type pt_lim: tuple

    """

    logger.info(" Starting the calculation of Afb...")

    # Create RDataFrame and cut on pt
    rdf_afb=rdf_all.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
        f"Cut on pt {pt_lim}").Define("Dimuon_y_abs", "abs(Dimuon_y)")

    # Create 4 different histograms:
    #   - h_Nf, forward events with weights w_n
    #   - h_Df, forward events with weights w_d
    #   - h_Nb, backward events with weights w_n
    #   - h_Db, backward events with weights w_d

    rdf_f = rdf_afb.Filter("Dimuon_cos>0", "Forward events")
    rdf_f.Report().Print()

    h_Nf = rdf_f.Histo2D(ROOT.RDF.TH2DModel("N_f", "N_f", 6, RAPIDITY_BIN,
        12, MASS_BIN), "Dimuon_y_abs", "Dimuon_mass", "w_n")

    h_Df = rdf_f.Histo2D(ROOT.RDF.TH2DModel("D_f", "D_f", 6, RAPIDITY_BIN,
        12, MASS_BIN), "Dimuon_y_abs", "Dimuon_mass", "w_d")

    rdf_b = rdf_afb.Filter("Dimuon_cos<0", "Backward events")
    rdf_b.Report().Print()

    h_Nb = rdf_b.Histo2D(ROOT.RDF.TH2DModel("N_b", "N_b", 6, RAPIDITY_BIN,
        12, MASS_BIN), "Dimuon_y_abs", "Dimuon_mass", "w_n")

    h_Db = rdf_b.Histo2D(ROOT.RDF.TH2DModel("D_b", "D_b", 6, RAPIDITY_BIN,
        12, MASS_BIN), "Dimuon_y_abs", "Dimuon_mass", "w_d")

    # Compute numerator N = Nf-Nb
    h_Nf.Add(h_Nb.GetPtr(), -1)

    # Compute denominator D = Df+Db
    h_Df.Add(h_Db.GetPtr(),1)

    # Compute Division A4 = D/N
    h_Nf.Divide(h_Df.GetPtr())

    # Compute Afb = (3/8)*A4
    h_Nf.Scale(3/8)

    # Line to draw at Afb = 0
    line = ROOT.TLine(60, 0, 120, 0)
    line.SetLineWidth(1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)

    ROOT.gStyle.SetTitleBorderSize(2)
    ROOT.gStyle.SetTitleFillColor(19)
    ROOT.gStyle.SetTitleSize(0.0, "t")
    ROOT.gStyle.SetTitleOffset(0.0, "t")

    # Plot results for each range of rapidity
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Plot')))
    c_afb = ROOT.TCanvas("c_afb", "c_afb")
    c_afb.Print(f"afb_{data_type}.pdf[", "pdf")

    file_afb = ROOT.TFile(f"hafb_{data_type}.root", "recreate")

    ROOT.gStyle.SetOptStat("mr")
    for i in range(0, 6, 1):
        ROOT.gStyle.SetOptStat(0)
        h_yn = h_Nf.ProjectionY(f"hy_{i}", i+1, i+1)
        h_yn.SetTitle(f"Afb[y({round(RAPIDITY_BIN[i],1)},"
            f"{round(RAPIDITY_BIN[i+1], 1)})]")
        h_yn.Draw("E")
        h_yn.Write()
        line.Draw()
        c_afb.Print(f"afb_{data_type}.pdf", f"Title:Afb[y({RAPIDITY_BIN[i]}"
            f",{RAPIDITY_BIN[i+1]})")
        c_afb.SaveAs(f"Afb[y({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}),"
            f"pt({pt_lim[0]},{pt_lim[1]})]_{data_type}.png")

    c_afb.Print(f"afb_{data_type}.pdf]", f"Title:Afb[y({RAPIDITY_BIN[i]},"
        f"{RAPIDITY_BIN[i+1]})")
    file_afb.Close()

    # Summary plot
    c_div = ROOT.TCanvas("divided", "divided", 1400, 900)
    c_div.Divide(6, 1, 0, 0)

    tex = ROOT.TLatex()
    tex.SetTextAlign(13)
    tex.SetTextSize(0.04)
    tex.DrawLatex(0.44, 0.96, "A_{FB} vs Mass")

    tex2 = ROOT.TLatex()
    tex2.SetTextSize(0.020)
    tex2.DrawLatex(0.03, 0.83, "A_{FB}")

    for k in range(1, 7, 1):
        c_div.cd(k)
        ROOT.gStyle.SetTitleBorderSize(0)
        ROOT.gStyle.SetTitleFillColor(0)
        ROOT.gStyle.SetTitleSize(0.08, "t")
        ROOT.gStyle.SetTitleOffset(0.4, "t")

        h_y2 = h_Nf.ProjectionY(f"hyd_{k}", k, k)
        h_y2.SetAxisRange(-0.8,0.8, "Y")
        h_y2.GetXaxis().SetLabelSize(0.05)
        h_y2.GetXaxis().SetTitle("Mass [GeV]")
        h_y2.GetXaxis().CenterTitle()
        h_y2.GetXaxis().SetTitleSize(0.05)
        h_y2.GetXaxis().ChangeLabel(1, -1, 0)
        h_y2.GetXaxis().ChangeLabel(7, -1, 0)
        h_y2.SetTitle("#scale[1.2]{"f"{round(0+(k-1)*0.4, 1)}<|y|<"
            f"{round(0.4+(k-1)*0.4 ,1)}""}")

        p = c_div.GetPad(k)
        p.SetPad(f"{k}",f"{k}",((k-1)*0.16)+(0.02),.05,((k-1)*0.16)+.18,.9,0,0,0)

        if k==1:
            p.SetLeftMargin(0.065)
            p.SetRightMargin(0)
            h_y2.GetYaxis().SetLabelOffset(0.01)
            h_y2.GetYaxis().SetLabelSize(0.04)
        elif k==7:
            p.SetLeftMargin(0)
        else:
            p.SetRightMargin(0)
            p.SetLeftMargin(0)

        h_y2.DrawClone("E")
        ROOT.gStyle.SetOptStat(0)
        line.Draw()

    c_div.SaveAs(f"afb_y_{data_type}.png")
    os.chdir(os.path.dirname(os. getcwd()))
    logger.info("Plots of Forward-Backward asymmetry have been saved.")


def comparison_cos(mc_file, data_file):
    """
    The function takes the histograms of Cos(theta*) obtained by \"cos_eta()\"
    function, for Monte Carlo and data. So previously, it's necessary to run the
    analysis for both. It rescales the entries of Monte Carlo with the data ones
    and plot both in the same canvas.
    Plots are named \"MC_data_cos(rap_inf,rap_sup).png\".

    :param mc_file: file containing histograms from MC
    :type mc_file: root file, required
    :param data_file: file containing histograms from data
    :type data_file: root file, required

    """

    # Open files containing histograms
    try:
        f_mc = ROOT.TFile(mc_file)
        f_all = ROOT.TFile(data_file)
    except OSError:
        logger.error("Root files containing histograms don't exist. Maybe you "
            "have to run the analysis first (w/o \"-no--a\" as argument).")
        sys.exit()

    c_cc = ROOT.TCanvas("MC+data_cos", "MC+data_cos")
    c_cc.Print("MC_data_cos.pdf[", "pdf")

    for i in range(0, 6, 1):
        ROOT.gStyle.SetTitleSize(0.045, "t")

        # Retrieve hisograms of MC and DATA
        h_mc = f_mc.Get(f"hc{i}")
        n_mc = h_mc.GetEntries()
        h_all = f_all.Get(f"hc{i}")
        n_all=h_all.GetEntries()

        # Normalized MC histo to DATA histo
        h_all.Scale(n_mc/n_all)

        # Draw and style
        h_mc.Draw()
        h_all.Draw("SAME P*")
        ROOT.gStyle.SetOptStat("mr")
        h_max = h_mc.GetMaximum()
        h_mc.SetAxisRange(0, h_max+(0.25*h_max), "Y")
        h_mc.SetFillColorAlpha(ROOT.kOrange+1, .7)
        h_all.SetMarkerColorAlpha(ROOT.kBlack, .9)

        leg = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
        leg.AddEntry(h_mc, "Monte Carlo histogram", "f")
        leg.AddEntry(h_all, "Data", "p")
        leg.SetShadowColor(0)
        leg.Draw()
        c_cc.Print("MC_data_cos.pdf",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")
        c_cc.SaveAs(f"MC_data_cos({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}).png")

    c_cc.Print("MC_data_cos.pdf]",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")

    # Close files
    f_mc.Close()
    f_all.Close()
    logger.info("Plots with the comparison between data and MC cos(theta*) "
        "distributions have been saved.")

def comparison_mass(mc_file, data_file):
    """
    The function takes the histograms of Mass obtained by \"mass_eta()\"
    function, for Monte Carlo and data. So previously, it's necessary to run the
    analysis for both. It rescales the entries of Monte Carlo with the data ones
    and plot both in the same canvas.
    Plots are named \"MC_data_mass(rap_inf,rap_sup).png\"

    :param mc_file: file containing histograms from MC
    :type mc_file: root file, required
    :param data_file: file containing histograms from data
    :type data_file: root file, required

    """

    # Open files containing histograms
    try:
        f_mc = ROOT.TFile(mc_file)
        f_all = ROOT.TFile(data_file)
    except OSError:
        logger.error("Root files containing histograms don't exist. Maybe you"
            "have to run the analysis first (w/o \"-no--a\" as argument.)")
        sys.exit()

    c_cm = ROOT.TCanvas("MC+data_mass", "MC+data_afb")
    c_cm.Print("MC_data_mass.pdf[", "pdf")

    for i in range(0, 6, 1):
        # Retrieve hisograms of MC and DATA
        h_mc = f_mc.Get(f"hm{i}")
        n_mc = h_mc.GetEntries()
        h_mc.SetName("Monte Carlo")
        h_all = f_all.Get(f"hm{i}")
        n_all=h_all.GetEntries()

        # Normalized MC histo to DATA histo
        h_all.Scale(n_mc/n_all)

        # Draw and style
        h_mc.Draw()
        h_all.Draw("SAME P*")
        ROOT.gPad.SetLogy()
        ROOT.gStyle.SetOptStat("mr")
        h_mc.SetFillColorAlpha(ROOT.kAzure+2, .7)
        h_all.SetMarkerColorAlpha(ROOT.kBlack, 1)
        h_all.SetMarkerSize(.8)
        leg = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
        leg.AddEntry(h_mc, "Monte Carlo histogram", "f")
        leg.AddEntry(h_all, "Data", "p")
        leg.SetShadowColor(0)
        leg.Draw()
        c_cm.Print("MC_data_mass.pdf",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")
        c_cm.SaveAs(f"MC_data_mass({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}).png")

    # Close files
    c_cm.Print("MC_data_mass.pdf]",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")
    f_mc.Close()
    f_all.Close()
    logger.info("Plots with the comparison between data and MC mass distribution"
        " have been saved.")

def comparison_afb(mc_file, data_file):
    """
    The function takes the histograms of Afb obtained by \"afb()\" function,
    for Monte Carlo and data. So previously, it's necessary to  run the analysis
    for both. It rescales the entries of Monte Carlo with the data ones
    and plot both in the same canvas.
    Plots are named \"MC_data_afb(rap_inf,rap_sup).png\"

    :param mc_file: file containing histograms from MC
    :type mc_file: root file, required
    :param data_file: file containing histograms from data
    :type data_file: root file, required

    """

    # Open files containing histograms
    try:
        f_mc = ROOT.TFile(mc_file)
        f_all = ROOT.TFile(data_file)
    except OSError:
        logger.error("Root files containing histograms don't exist. Maybe you"
            "have to run the analysis first (w/o \"-no--a\" as argument.)")
        sys.exit()


    c_cm = ROOT.TCanvas("MC+data_afb", "MC+data_afb")
    c_cm.Print("MC_data_afb.pdf[", "pdf")

    for i in range(0, 6, 1):
        # Retrieve hisograms of MC and DATA
        h_mc = f_mc.Get(f"hy_{i}")
        h_mc.SetName("Monte Carlo")
        h_all = f_all.Get(f"hy_{i}")

        # Draw and style
        h_all.Draw()
        h_mc.Draw("SAME")
        ROOT.gStyle.SetOptStat("mr")
        h_max = h_mc.GetMaximum()
        h_min = h_all.GetMinimum()
        h_all.SetAxisRange(h_min-(0.55*abs(h_min)), h_max+(0.45*h_max), "Y")
        h_mc.SetMarkerColorAlpha(ROOT.kAzure+2, .8)
        h_mc.SetMarkerStyle(ROOT.kFullTriangleUp)
        h_all.SetMarkerColor(ROOT.kBlack)
        leg = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
        leg.AddEntry(h_mc, "Monte Carlo", "p")
        leg.AddEntry(h_all, "Data", "p")
        leg.SetShadowColor(0)
        leg.Draw()
        c_cm.Print("MC_data_afb.pdf",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")
        c_cm.SaveAs(f"MC_data_afb({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]}).png")

    # Close files
    c_cm.Print("MC_data_afb.pdf]",f"Title:({RAPIDITY_BIN[i]},{RAPIDITY_BIN[i+1]})")
    f_mc.Close()
    f_all.Close()
    logger.info("Plots with the comparison between data and MC of Forward-Backward"
        " asymmetry have been saved.")



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

    # Load the shared library "tools.cpp" which contains some functions to
    # calculate the useful quantities for the analysis.
    ROOT.gSystem.Load('../Utils/tools_cpp.so')

    # Start the the timer
    start = time.time()

    # Creating the parser
    parser = argparse.ArgumentParser(description = "Processing the root file of data.")
    parser.add_argument("-f", "--file", type=str, nargs="+",
        help="The name of the file index which contains the path of the different"
             " root files to analyze. E.g.:\"Run2012B_SingleMu_index.txt\""
             "Put the string \"_index\" before the extension \".txt\".")
    parser.add_argument("-t", "--d_type", type=str, nargs="+",
        help="The type of data to analyzed. In this analysis this parameter is "
             "intended to be among \"MC (Monte Carlo)\" and \"data (actual data)\".")
    parser.add_argument("-no--a", "--analysis", action="store_false",
        default=True, help="If True, retrieve data from internet, otherwise not."
        "If False, it takes data files already analyzed (from the directory "
        "\"Snapshot\") and listed in \Analyzed_files.txt\".Default value is True.")
    parser.add_argument("-c", "--comparison", action="store_true",
        default=False, help="If True, it runs the functions that make plots with"
        " the comparison between MC and actual data (cos(theta*), mass and Afb "
        "distributions). It's necessary to have run the functions \"cos_eta\", "
        "\"mass_eta\" and \"afb\" separatly on them.")
    args = parser.parse_args()

    # Creating the logger
    logger = utils.set_logger("Z analysis", logging.DEBUG)

    # Z ANALYSIS
    logger.info(" Starting the analysis of the Forward-Backward asymmetry"
                " of the Z resonance.")

    # A new folder to collect plots and data is created
    os.makedirs("Plot", exist_ok=True)
    os.makedirs("Snapshots", exist_ok=True)
    logger.debug(" The new directories \"Plot\" and \"Snapshots\" are created.")

    # Enable the multi-threading analysis
    N_THREADS = 8
    ROOT.ROOT.EnableImplicitMT(N_THREADS)

    # Retrieve dataset from the web and start the selection on good events.
    if args.analysis is True:
        root_files=ROOT.std.vector("string")()
        times = []
        N_times = []
        dtypes = []
        for f, ty in zip(args.file, args.d_type):
            rootf_ana, ftime, N_time, dtypef = retrieve_dataset(f, ty)
            root_files+=rootf_ana
            times+=ftime
            N_times+=N_time
            dtypes+=dtypef
        with open("Analyzed_files.txt", "w+", encoding="utf-8") as outf1:
            print("file/C:type/C", file=outf1)
            for r, tye in zip(root_files, dtypes):
                print(r, tye, file=outf1)
        logger.debug("\"Analyzed_files.txt\" has been saved.")

        # Save and plot histogram time vs N° events analyzed for each file
        with open("times_vs_N.txt", "w+", encoding="utf-8") as outf2:
            print("time:N_events", file=outf2)
            for t, N_t in zip(times, N_times):
                print(t, N_t, file=outf2)
        logger.debug("\"times_vs_N.txt\" has been saved.")

        h_times = ROOT.TGraph("times_vs_N.txt")
        h_times.SetTitle(f"Times [n_threads set = {N_THREADS}];t [s];Events")
        c_times = ROOT.TCanvas("Times", "Times")
        h_times.Draw("AP")
        c_times.SaveAs("times_vs_N.png")
        logger.debug(" Plot times vs. N has been saved.")
    elif args.analysis is False:
        if not os.path.isfile('Analyzed_files.txt'):
            raise FileNotFoundError("The file \"Analyzed_files.txt\" doesn't "
            "exist. Maybe you have to run the analysis first.")
        logger.info("No analysis on data has been done. Files from "
            "\"Analyzed_files.txt\" are considered in the following...")

    # Create RDataFrames
    t_ana=ROOT.TTree()
    try:
        nline=t_ana.ReadFile("Analyzed_files.txt")
    except OSError as err1:
        logger.error(f"{err1} \nSomething went wrong with the file "
            "\"Analyzed_files.txt\". Try to run again the analysis.")

    if nline==0:
        logger.error("\"Analyzed_files.txt\" is empty! No files to analyze.")
        sys.exit()

    #t_ana.Scan()
    rfile_temp=[]
    d_types_temp = []

    for t in t_ana:
        rfile_temp.append(t_ana.file)
        d_types_temp.append(t_ana.type)
    d_types = np.asarray(d_types_temp)

    tdata = {}
    uniques = np.unique(d_types)
    for u in uniques:
        root_files = ROOT.std.vector("string")()
        for f, t in zip(rfile_temp,d_types):
            if t==u:
                root_files.emplace_back(f)
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Snapshots')))
        all_data = ROOT.RDataFrame("dimuon_w", root_files).Cache()
        N_tot = all_data.Count().GetValue()
        logger.info(f" Total number of events that passed the cuts is {N_tot}.")
        os.chdir(os.path.dirname(os. getcwd()))
        tdata[f"{u}"]= all_data

    # Move the snapshots in a dedicated directory
    if args.analysis is True:
        for f in rfile_temp:
            os.replace(f'{os.getcwd()}/{f}', f'Snapshots/{f}')

    # Distributions of mass and cos(theta*) in six different range of rapidity.
    # "pt_lim" is a tuple to set the limits of transverse momentum.
    for dt, da in zip(tdata.keys(), tdata.values()):
        pt_lim = (0,round(da.Max("Dimuon_pt").GetValue(),2))
        # mass_eta(da,dt,pt_lim)
        # cos_eta(da,dt,pt_lim)

        # Compute AFB (Forward_Backward Asymmetry) and plot results
        logger.info(" Computing Forward-Backward asymmetry...")
        # afb(da,dt,pt_lim)

    # Comparison plot
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Plot')))
    if args.comparison:
        try:
            comparison_cos(f"hcos_{uniques[0]}.root", f"hcos_{uniques[1]}.root")
            comparison_mass(f"hmass_{uniques[0]}.root", f"hmass_{uniques[1]}.root")
            comparison_afb(f"hafb_{uniques[0]}.root", f"hafb_{uniques[1]}.root")
        except IndexError as err2:
            logger.error(f"{err2}:\nThese functions are able to compare results"
                " for two set of data, for example Monte Carlo and collected data."
                "It seems that you have analyzed only one or none.")
            sys.exit()

    os.chdir(os.path.dirname(os. getcwd()))

    # Elapsed time
    logger.info(f" Elapsed time from the beginning is:"
        f" {time.time()-start}(real time) seconds")
