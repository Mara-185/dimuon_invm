"""
In this module a few functions are implemented in order to study the
\"Forward-Backward asymmetry" for the Z boson decay.
The analysis is done in six different rapidity ranges of equal size and
twelve mass bins:

    - Mass bins : 60 < Mass < 120:

    [60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120];

    - Rapidity bins of equal size for \|rapidity\| < 2.4:

    [0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]

The script takes as argument:

    - the data file (URL) of dileptons (-f), for example: "root:
        //eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/
        Run2012B_DoubleMuParked.root";
    - the range of files to analyze (-l).

In the analysis the following version have been used:

    - Python v3.8

    - ROOT v6.24 ("source ~/root/bin/thisroot.sh" command needed before starting
        the analysis to set the environment of ROOT)

"""

import ROOT
import logging
import os
import argparse
import time
import sys
from pathlib import Path

# Add my modules to the path
root_utils = os.path.abspath('../Utils')
sys.path.insert(0, root_utils)
import utils

#logger = logging.getLogger(__name__)

def retrieve_dataset(findex, start, stop):
    """
    The function takes as argument the file index with all name of root files
    of the chosen dataset and the two extremities of the files' range to analyze.
    (The file index of the chosen dataset has to be downloaded)
    Two lists and a standard vector are returned in order to save the time
    needed to analyze each file, its number of events and the name of the
    snapshot created in the "z_main" function.

    :param findex: index of root files which make up the entire dataset
    :type findex: string of txt file, required
    :param start: first index of the range of files to analyze
    :type start: int, required
    :param stop: second index of the range of files to analyze
    :type stop: int, required
    :return files: list of the files name analyzed.
    :rtype files: standard vector of string
    :return times: list of time elpased analyzing each file
    :rtype times: list
    :return N: list of events analyzed for each file
    :rtype N: list

    """

    with open(findex) as inf:
        # Strings of all files in the dataset are loaded.
        all_dataf = [_.rstrip('\n') for _ in inf]
        run = findex.replace("_index.txt", "")
        logger.info(f" Total number of data files in {run}: {len(all_dataf)}\n")

        # Only the chosen range of files are selected.
        data = all_dataf[start:stop]
        times = []
        N = []
        files = ROOT.std.vector("string")()

        for i,d in enumerate(data):
            logger.info(f" Analyzing file number {i+1} :\n{d}")

            start_analysis = time.time()
            string_file, N_events = z_main(d, i, run)
            stop_analysis = time.time() - start_analysis
            logger.info(f" Elapsed time : {stop_analysis}")

            files.push_back(string_file)
            N.append(N_events)
            times.append(stop_analysis)

    return files, times, N



def z_main(url, iteration, run):
    """
    The function create an RDataFrame to make the selection of the good events
    for the analysis, to compute some useful quantities and store them in new
    columns. It returns the string of the root file created by a snapshot only
    of the columns needed.

    :param url: url of the root file to upload from the web.
    :type url: url of root file, required.
    :param iteration: number of file analyzed
    :type iteration: int, required
    :return snap_name: name of the root file created
    :rtype snap_name: string
    :return N_ev: number of all events in the file
    :rtype N_ev: float

    """

    logger.info("Start the analysis...")

    # Create an RDataFrame to make selection on data from the root file.
    root_file = ROOT.TFile.Open(url, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug(" The RDataFrame have been created.")

    # Count of all data in the file.
    N_ev = rdf.Count().GetValue()

    rdf_cut=rdf.Filter("nMuon>1","Selection on events with at least two muons").\
        Define("mask", "Muon_isTracker==1 && Muon_isGlobal==1 && "
        "Muon_gChi2<10 && abs(Muon_dxyBest)<0.2 && Muon_pfRelIso03_all<0.1 && "
        "abs(Muon_eta)<2.4").\
        Define("Mu_pt", "Muon_pt[mask]").\
        Define("Mu_pterr", "Muon_ptErr[mask]").\
        Define("Mu_mass", "Muon_mass[mask]").\
        Define("Mu_charge", "Muon_charge[mask]").\
        Define("Mu_phi", "Muon_phi[mask]").\
        Define("Mu_eta", "Muon_eta[mask]").\
        Filter("Mu_pt.size()>=nMuon").\
        Filter("Mu_pt[0]>25 && Mu_pt[1]>15",
               "Trigger on the leading and the next leading muons").\
        Filter("Mu_charge[0]!=Mu_charge[1]",
               "Selection on muons with opposite charge").\
        Define("Dimuon_mass", \
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[3]").\
        Define("Dimuon_pt",\
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[0]").\
        Define("Dimuon_eta", \
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]").\
        Define("Dimuon_phi",\
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[2]").\
        Define("Dimuon_cos",\
            "cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[0]").\
        Define("Dimuon_y",\
            "cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]").\
        Filter("Dimuon_mass>60 && Dimuon_mass<120","Cut on Z resonance").\
        Define("w_d", f"weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[0]").\
        Define("w_n", f"weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[1]")#.\
        #Snapshot("dimuon_w", snap_name, branchlist)

    logger.info("\nReport of all cuts:\n")
    rdf_cut.Report().Print()

    branchlist = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi",\
        "Dimuon_cos","Dimuon_y","w_d", "w_n"]
    snap_name = f"dimuon_w{iteration}[{run}].root"

    # A snapshot is done to collect the useful physiscal quantity of the
    # dileptons in a root file and a node graph is saved.
    ROOT.RDF.SaveGraph(rdf_cut, "rdf_main2.dot")
    rdf_cut.Snapshot("dimuon_w", snap_name, branchlist)
    logger.info(" The snapshot is done.")

    root_file.Close()

    return snap_name, N_ev


def mass_eta(rdf_all, rap_lim, pt_lim = (0,120), infile=None):
    """
    Plot the Z resonance (mass range from 60 GeV to 120 GeV), in the six
    different range of rapidity. It also make a cut on pt.
    Plot is saved as "Dimuon_mass(rap_inf,rap_sup)_(pt_inf,pt_sup).png" in a
    directory named "Z analysis".
    Data can be passed as RDataFrame or as root file, but the name of the TTree
    has to be the same of the file.

    :param rdf_all: RDataFrame to analyze.
    :type rdf_all: RDataFrame, required
    :param rap_lim: range limits of rapidity
    :type rap_lim: tuple, required
    :param pt_lim: range limits of transverse momentum
    :type pt_lim: tuple, required, default=(0,120)
    :param infile: data file
    :type infile: root file, not required

    """

    # RDataFrame with appropriate cuts is created.
    if infile == None:
        rdf_m=rdf_all.\
            Filter(f"abs(Dimuon_y)<{rap_lim[1]} && abs(Dimuon_y)>{rap_lim[0]}",\
                    "Cut on rapidity").\
            Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                    "Cut on pt")
        rdf_m.Report().Print()
    else:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.\
            Filter(f"abs(Dimuon_y)<{rap_lim[1]} && abs(Dimuon_y)>{rap_lim[0]}"). \
            Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}")
        del rdf

    logger.debug(" The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Mass in defined range rapidity",
        f"Mass in y [{rap_lim[0]},{rap_lim[1]}];Mass [MeV];Events", 60, 60, 120),
        "Dimuon_mass")
    logger.info(f" Drawing the mass distribution for {rap_lim[0]}< "
                f"|Dimuon rapidity| < {rap_lim[1]}")

    # Styling
    c = ROOT.TCanvas("Z mass", "Z mass")
    c.SetGrid()
    ROOT.gPad.SetLogy()
    ROOT.gStyle.SetOptStat("e")
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"Dimuon_mass({rap_lim[0]},{rap_lim[1]})_({pt_lim[0]},{pt_lim[1]}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def cos_eta(rdf_all, rap_lim, pt_lim=(0,120), infile=None):
    """
    Plot the cos of Z resonance (mass range from 60 GeV to 120 GeV), in the six
    different range of rapidity. It also make a cut on pt.
    Plot is saved as "Dimuon_cos(rap_inf,rap_sup)_(pt_inf,pt_sup).png" in a
    directory named "Z analysis".
    Data can be passed as RDataFrame or as root file, but the name of the TTree
    has to be the same of the file.

    :param rdf_all: RDataFrame to analyze.
    :type rdf_all: RDataFrame, required
    :param rap_lim: range limits of rapidity
    :type rap_lim: tuple, required
    :param pt_lim: range limits of transverse momentum
    :type pt_lim: tuple, required, , default=(0,120)
    :param infile: data file
    :type infile: root file, not required

    """

    # RDataFrame with appropriate cuts is created.
    if infile == None:
        rdf_m=rdf_all.\
            Filter(f"abs(Dimuon_y)<{rap_lim[1]} && abs(Dimuon_y)>{rap_lim[0]}",\
                    "Cut on rapidity").\
            Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}",
                    "Cut on pt")
        rdf_m.Report().Print()
        del rdf_all
    else:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.\
            Filter(f"abs(Dimuon_y)<{rap_lim[1]} && abs(Dimuon_y)>{rap_lim[0]}"). \
            Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}")
        del rdf

    logger.debug(" The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Cos in defined range rapidity",
        f"cos(#theta *) in y [{rap_lim[0]},{rap_lim[1]}];Cos;Events",40,-1,1),
        f"Dimuon_cos")
    logger.info(f" Drawing the cos(theta*) distribution for {rap_lim[0]}< "
                f"|Dimuon rapidity| < {rap_lim[1]}")

    # Styling
    c = ROOT.TCanvas("Z cos", "Z cos")
    c.SetGrid()
    ROOT.gStyle.SetOptStat("e")
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"Dimuon_cos({rap_lim[0]},{rap_lim[1]})_({pt_lim[0]},{pt_lim[1]}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def afb(rdf0, pt_lim):
    """
    The function calculates the value of Afb (Asymmetry forward-backward)
    using the \"angular event weighting\" in the six different ranges of
    rapidity, for each mass.
    It creates different \".txt\" files with the results.

    :param rdf0: RDataframe to compute Afb
    :type rdf0: RDataFrame, required
    :param pt_lim: range limits of pt
    :type pt_lim: tuple, required

    """

    logger.info(" Starting the calculation of Afb...")
    t_name="dimuon_w"

    # Create RDataFrame
    rdf=rdf0.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}")
    rdf.Report().Print()

    logger.debug(" The RDataFrame has been created.")


    for j in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve rapidity range bin from dictionary RAPIDITY_BIN
        rap_lim = utils.RAPIDITY_BIN[f"{j}"]

        rdf_y = rdf.Filter(f"abs(Dimuon_y)<{rap_lim[1]} && abs(Dimuon_y)>{rap_lim[0]}")

        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))
        with open(f"Afb[y({rap_lim[0]},{rap_lim[1]})]_({pt_lim[0]}, {pt_lim[1]})"
            f"{t_name}.txt", "w+") as outf:
            print("Afb:M:A4", file=outf)
            os.chdir(os.path.dirname(os. getcwd()))
            for i in range(0, len(utils.MASS_BIN), 1):

                # Retrive mass range bin from dictonary MASS_BIN
                m_lim = utils.MASS_BIN[f"{i}"]

                rdf_m = rdf_y.Filter(f"Dimuon_mass>{m_lim[0]} && "
                    f"Dimuon_mass<{m_lim[1]}")

                rdf_f = rdf_m.Filter(f"Dimuon_cos>0")
                D_f = rdf_f.Sum("w_d")
                N_f = rdf_f.Sum("w_n")

                rdf_b = rdf_m.Filter(f"Dimuon_cos<0")
                D_b = rdf_b.Sum("w_d")
                N_b = rdf_b.Sum("w_n")

                #######################GESTISCI ZeroDivisionError
                try:
                    Afb = (3*(N_f.GetValue()-N_b.GetValue()))/(8*(D_f.GetValue()+D_b.GetValue()))
                    A4 = Afb*(8/3)
                except ZeroDivisionError:
                    Afb = -100
                    A4 = -100
                    logger.error(" A ZeroDivisionError occured. Values are set"
                                 " to \"-100\".")

                M = m_lim[0] +((m_lim[1]-m_lim[0])/2)

                print(Afb, M, A4,file=outf)

    #del rdf_y, rdf_m, rdf_f, rdf_b, rdf


def afb_plot(infile, pt_lim):
    """
    Plot the values of Afb obtain with the function \"afb\".
    It takes as input a txt file with the name of the txt files created by
    function \"afb\". This file has to be created before using this method.

    :param infile: txt with the list of files you want to reproduce Afb results
    :type infile: txt file, required
    :param pt_lim: range limits of pt
    :type pt_lim: tuple, required

    """

    with open(infile) as ft:
        data = [_.rstrip('\n') for _ in ft]
    c = ROOT.TCanvas("Afb")
    file_name = infile.replace("txt", "pdf")
    c.Print(f"{file_name}[")
    for inf in data:
        # Retrieve data
        name = inf.replace(f"_({pt_lim[0]}, {pt_lim[1]})dimuon_w.txt", "")
        t = ROOT.TTree("tree", "tree")
        t.ReadFile(inf)

        # Styling
        c.UseCurrentStyle()
        c.SetGrid()
        ROOT.gStyle.SetOptStat("e")
        t.Draw("Afb:M","Afb!=-100")
        graph = ROOT.TGraph()
        graph = ROOT.gPad.GetPrimitive("htemp")
        graph.GetXaxis().SetTitle("M [GeV]")
        graph.GetXaxis().CenterTitle()
        graph.GetYaxis().SetTitle("Afb [a.u.]")
        graph.GetYaxis().CenterTitle()
        graph.SetTitle(f"{name}")

        c.Print(f"{file_name}", f"Title:{name}")
        c.SaveAs(f"{name}.png")
    c.Print(f"{file_name}]", f"Title:{name}")



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
    start_cpu = time.process_time()

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str, nargs="+",
        help="The name of the file index to pass to the script. "
             "It has to contain the path of the different root files to analyze,"
             "which are part of the dataset.E.g.:\"Run2011A_SingleMu_index.txt\"")
    parser.add_argument("-d", "--display", type=bool, help="Show RDataFrame")
    parser.add_argument("-l", "--limit", nargs="+",type=int,
        help="Range of files to analyze, taken from the file index."
             " Format: -l start stop. \"stop\" is excluded.")
    parser.add_argument("-no--a", "--analysis", action="store_false", default=True,
        help="If True, retrieve data from internet, otherwise not. "
        "Default value is True. The range of files selected must be already "
        "analyzed!! Otherwise the script don't consider them.")
    args = parser.parse_args()

    # Creating the logger
    logger = utils.set_logger("Z analysis", logging.DEBUG)

    # Z ANALYSIS
    logger.info(" Starting the analysis of the Forward-Backward asymmetry"
                " of the Z resonance.")

    # A new folder to collect plots and data is created
    os.makedirs("Z analysis", exist_ok=True)
    os.makedirs("Snapshots", exist_ok=True)
    logger.debug(" The new directories \"Z analysis\" and \"Snapshots\" are created.")

    # Enable the multi-threading analysis
    if not args.display:
        nthreads = 10
        ROOT.ROOT.EnableImplicitMT(nthreads)

    # Check on the validity of the inserted input (-l).
    #???????????????????????????????????????????????
    lim=[]
    N_tot_files = 0
    try:
        for n in range(0, len(args.limit), 2):
            limits = (args.limit[n], args.limit[n+1])
    except ValueError:
        logger.error("The inserted values are not right. Check on the input "
                      "signature in the help's parser. First index has to be"
                      "greater than the second.")
    else:
        lim.append(limits)
        N_tot_files+=(args.limit[n+1] - args.limit[n])

    # Check the correspondance between number of index files and limits inserted.
    if len(args.file)!=len(lim):
        logger.error(" There isn't correspondance between number of index "
            "files and ranges inserted. Please, check on help's parser.")
        sys.exit(1)

    # Retrieve dataset from the web and start the selection on good events.
    if args.analysis==True:
        root_files=ROOT.std.vector("string")() #allocare gli slot giusti?
        root_files.reserve(N_tot_files) ####????
        times = []
        N_times = []
        for f, limits in zip(args.file, lim):
            root_files0, ftime, N_time = retrieve_dataset(f, limits[0], limits[1])
            root_files+=root_files0
            times+=ftime
            N_times+=N_time
            with open("Analyzed_files.txt", "w") as outf1:
                for r in root_files:
                    print(r, file=outf1)

        # Save and plot histogram time vs N° events anayzed for each file
        with open(f"times_vs_N.txt", "w+") as outf2:
            print("time:N_events", file=outf2)
            for t, N_t in zip(times, N_times):
                print(t, N_t, file=outf2)

        h_times = ROOT.TH1D("h_times", f"Times [n_threads set = {nthreads}];"
            "t [s];Events", 10, 0, 200)
        for t, N_time in zip(times, N_times):
            h_times.Fill(t,N_time)
        c_times = ROOT.TCanvas("Times", "Times")
        h_times.Draw()
        c_times.SaveAs("times_vs_N.png")
        logger.info(" Plot times saved.")

        # Create a RDataFrame with all the selected data.
        all_data = ROOT.RDataFrame("dimuon_w", root_files).Cache()
        logger.info(f" Total number of events that passed the cuts is "
                    f"{all_data.Count().GetValue()}.")
        print(all_data.GetNSlots())

    elif args.analysis==False:
        root_files = ROOT.std.vector("string")()
        with open("Analyzed_files.txt") as inf:
            for line in inf:
                line = line.replace("\n", "")
                root_files.push_back(line)

        # Create a RDataFrame with all the selected data.
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Snapshots')))
        all_data = ROOT.RDataFrame("dimuon_w", root_files).Cache()
        N_tot = all_data.Count().GetValue()
        logger.info(f" Total number of events that passed the cuts is {N_tot}.")
        os.chdir(os.path.dirname(os. getcwd()))

    else:
        logger.error(" Invalid input for \"-a\" (analysis) option. "
                    "Possibilities are: True or False.")
        sys.exit(1)

    # Move the snapshots in a dedicated directory
    if args.analysis==True:
        start_move = time.time()
        for f in root_files:
            os.replace(f'{os.getcwd()}/{f}', f'Snapshots/{f}')
        logger.info(f" Time to move the files: {time.time()-start_move}")

    # Plot the distributions of mass and cos(theta*) in six different range of
    # rapidity. "pt_lim" is a tuple to set the limits on transverse momentum.
    pt_lim = (10,100)

    for i in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve the values of eta range bins (y_inf, y_sup) from "utils.py"
        y_lim = utils.RAPIDITY_BIN[f"{i}"]
        mass_eta(all_data,y_lim,pt_lim)
        cos_eta(all_data,y_lim,pt_lim)

    # Compute AFB (Forward_Backward Asymmetry)
    afb(all_data, pt_lim)

    # Plot AFB results
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))

    # Create a txt with a list of the files from which takes values.
    pathlist = Path("").glob(f"*_{pt_lim}dimuon_w.txt")
    with open(f"Afblist_{pt_lim}.txt", "w+") as outf:
        for path in pathlist:
            print(path, file=outf)

    logger.info(" Plotting Afb results...")
    afb_plot(f"Afblist_{pt_lim}.txt", pt_lim)

    os.chdir(os.path.dirname(os. getcwd()))

    # Elapsed time
    logger.info(f" Elapsed time from the beginning is:"
        f" {time.time()-start}(real time) seconds")
    logger.info(f" Elapsed time from the beginning is:"
        f" {time.process_time()-start_cpu}(cpu time) seconds")
