"""In this module a few functions are implemented, in order to study the
    \"Forward-Backward asymmetry" for the Z boson decay.
    The analysis is done on six different rapidity ranges of equal size and
    twelve mass bins:
    - Mass bins from article : 60 < M < 120:
        [60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120];
    - Rapidity bins of equal size for |y| < 2.4:
        [0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
"""

import ROOT
import logging
import os
import utils
import argparse
import time
from pathlib import Path


#logger = logging.getLogger(__name__)

def retrieve_dataset(findex, start, stop):
    """
    The function takes as argument the file index with all name of root files
    of the chosen dataset and the two extremities of the file's range to analyze.
    The file index of the chosen dataset has to be downloaded.

    :param findex : index of root files which make up the entire dataset
    :type findex: txt file, required.
    :param start: first index of the range of files to analyze
    :type start: int, required
    :param stop: second index of the range of files to analyze
    :type stop: int, required
    """

    with open(findex) as inf:
        all_data = [_.rstrip('\n') for _ in inf]
        run = findex.replace("_index.txt", "")
        logger.info(f"Total number of data files in run {run}: {len(all_data)}\n")
        data = all_data[start:stop]
        times = []

        root_files = ROOT.std.vector("string")()
        for i,d in enumerate(data):
            logger.info(f"Analyzing file number {i+1} :\n{d}")
            start1 = time.time()
            root_files.push_back(z_main(d, i))
            stop1 = time.time()
            times.append(stop1-start1)
            logger.info(f"Time : {stop1-start1}")

    return root_files, times


def z_analysis(infile, iteration=0):
    """
    The function create an RDataFrame to do the right cuts on data and it creates
    a snapshot of the events that passed the cuts.
    ????????????????????????????????????????????????????????

        - infile: a vector of strings which contains the files' names to analyze.
            This is the output of the function \"retrieve_dataset\".
        - iteration: the index of the file which is been analyzing ???

    :param infile: vector of strings with the name of root files
    :type infile: vector, required
    :param iteration: number of file analyzed
    :type iteration: int, required
    """

    logging.info("Start the analysis...")

    # Create an RDataFrame of useful data from the root file.
    root_file = ROOT.TFile.Open(infile, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug("The RDataFrame have been created.")

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
               "Selection on muons with opposite charge")

######################CHECK ON MU.PT SIZE and PTERR #####################

    ###########################################################################
    ####################    SHARED LIBRARY  ###################################
    # Devo per forza fare il load altrimenti la dichiaro ma non la definisco; se non metto "gInterpreter"
    # invece, è come se non la dichiarassi.
    # Ha senso fare la shared library più del jitting? sì, perchè g++ permette di decidere
    # l'ottimizzazione, il jitting usa la più bassa per essere più veloce.
    # Ma io che ho fatto in ROOT come la specifico?
    if iteration==0:
        ROOT.gInterpreter.ProcessLine('#include "tools.h"')
        ROOT.gSystem.Load('./tools_cpp.so')
        logging.info("The shared library has been uploaded.")

    # A list of the columns to cache.
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi",
        "Dimuon_cos", "Dimuon_y"]

    rdf_dimu = rdf_cut.\
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
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]").Cache(branchlist_mu)
    del rdf

    ROOT.RDF.SaveGraph(rdf_dimu, "rdf_mu.dot")
    rdf_dimu.Report().Print()
    root_file.Close()

    return rdf_dimu


def weight(data_cached, iteration):
    """ The function computes the numerator and denominator weights for each
    dileptons, which are necessary to obtain the forward-backward asymmetry.

    :param infile: Root file obtained by the function \"z_analysis\".
    :type infile: root file, required
    :param iteration: number of file analyzed
    :type iteration: int, required

    """
    # Filter on mass in the range of Z resonance.
    rdf=data_cached.Filter("Dimuon_mass>60 && Dimuon_mass<120","Cut on Z resonance")

    rdf.Report().Print()

    rdf_w = rdf.\
        Define("w_d", f"weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[0]").\
        Define("w_n", f"weights(Dimuon_pt,Dimuon_mass, Dimuon_cos)[1]")
    del data_cached, rdf

    # A list with the names of the new columns is created in order to make a
    # snapshot of them
    branchlist = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi",\
        "Dimuon_cos","Dimuon_y","w_d", "w_n"]

    # A snapshot is done to collect the useful physiscal quantity of the dileptons
    # in a single root file.
    logger.debug("Starting snapshot...")
    rdf_w.Snapshot("dimuon_w", f"dimuon_w{iteration}.root", branchlist)
    logger.info("The snapshot is done.")

    return str(f"dimuon_w{iteration}.root")


def z_main(url, n):
    """The function makes the selection of the good events for the analysis and
    also compute some useful values. It returns a string with the name of the
    root file created at the end of the analysis.

    :param url: the url of the root file uploaded from the web.
    :type url: url of root file, required.
    :param n: number of the file analyzed
    :type n: int, required
    """
    # Selection on events
    rdf = z_analysis(url, n)

    # Cumpute the weights
    data_string = weight(rdf,n)

    return data_string


def mass_vs_eta(rdf0, rap_lim, pt_lim, bin, infile=None):
    """
    Plot the Z resonance (mass range from 60 GeV to 120 GeV), for the six
    different range of rapidity.
    Data can be passed as RDataFrame, as data cached or as root file. In every
    case, the name of the TTree has to be "dimuon_w".

        - rdf0: an RDataFrame wich contains all data, from the different file.
            It has to be created outside the function.
        ...

    :param rdf0: RDataFrame to analyze.
    :type rdf0: RDataFrame, required
    :param rap_lim: range limits of rapidity
    :type rap_lim: tuple, required
    :param pt_lim: range limits of pt
    :type pt_lim: tuple, required
    """

    # RDataFrame with appropriate cuts is created.
    if infile == None:
        rdf_m=rdf0.\
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

    logger.debug("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Mass in defined range rapidity",
        "Mass in defined range rapidity;Mass [MeV];Events", bin, 60, 120),
        "Dimuon_mass")
    logger.info(f"Drawing the mass distribution for {rap_lim[0]}< "
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

def cos_vs_eta(rdf0, rap_lim, pt_lim, bin, infile=None):
    """
    Plot the cos of Z resonance (mass range from 60 GeV to 120 GeV), for the six
    different range of rapidity.
    Data can be passed as RDataFrame, as data cached or as root file. In every
    case, the name of the TTree has to be "dimuon_w".

        - rdf0: an RDataFrame wich contains all data, from the different file.
            It has to be created outside the function.
        ...

    :param rdf0: RDataFrame to analyze.
    :type rdf0: RDataFrame, required
    :param rap_lim: range limits of rapidity
    :type rap_lim: tuple, required
    :param pt_lim: range limits of pt
    :type pt_lim: tuple, required
    :param data_cached: range limits of rapidity
    :type rap_lim: tuple, required
    """

    # RDataFrame with appropriate cuts is created.
    if infile == None:
        rdf_m=rdf0.\
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

    logger.debug("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("Cos in defined range rapidity",
        "cos(#theta *) in defined range rapidity;Cos;Events",bin,-1,1),f"Dimuon_cos")
    logger.info(f"Drawing the cos(theta*) distribution for {rap_lim[0]}< "
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



def afb(rdf0, pt_lim, infile=None):
    """ The function calculates the value of Afb (Asymmetry forward-backward)
     using the \"angular event weighting\" in the six different ranges of
     rapidity, for each mass.
     It creates different \".txt\" files with the results.
    """

    logger.info("Starting the calculation of Afb...")
    t_name="dimuon_w"

    # Create RDataFrame
    if infile==None:
        rdf=rdf0.Filter(f"Dimuon_pt>{pt_lim[0]} && Dimuon_pt<{pt_lim[1]}")
        rdf.Report().Print()
    else:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name, infile).Filter(f"Dimuon_pt>{pt_lim[0]} && "
            f"Dimuon_pt<{pt_lim[1]}")

    logger.debug("The RDataFrame has been created.")

    #
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
                    logger.error("A ZeroDivisionError occured. Values are set"
                                 " to \"-100\".")

                M = m_lim[0] +((m_lim[1]-m_lim[0])/2)

                print(Afb, M, A4,file=outf)

    del rdf_y, rdf_m, rdf_f, rdf_b


def afb_plot(infile, pt_lim):
    """ Plot the values of Afb obtain with the function \"afb\".
    It takes as input a txt file with the name of the txt files created from
    \"afb\". This file is created in the script.
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

    # Start the the timer
    timer = ROOT.TStopwatch()
    timer.Start(ROOT.kFALSE)

    start =time.time()

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str,
        help="The name of the file index to pass to the script. "
             "It has to contain the path of the different root files to analyze,"
             "which are part of the dataset.E.g.:\"Run2011A_SingleMu_index.txt\"")
    parser.add_argument("-d", "--display", type=bool, help="Show RDataFrame")
    parser.add_argument("-l", "--limit", nargs="+",type=int,
        help="Range of files to analyze, taken from the file index."
             " Format: -l start stop. \"stop\" is excluded.")
    parser.add_argument("-a", "--analysis", type=str, default=str(True),
        nargs="?", help="If True, retrieve data from internet, otherwise not. "
        "Default value is True.")
    args = parser.parse_args()

    # Creating the logger
    logger = utils.set_logger("Z analysis")

    # Z ANALYSIS
    logger.info("Starting the analysis of the Forward-Backward asymmetry"
                " of the Z resonance.")

    # A new folder to collect plots and data is created
    os.makedirs("Z analysis", exist_ok=True)
    os.makedirs("Snapshots", exist_ok=True)
    logger.debug("The new directories \"Z analysis\" and \"Snapshots\" are created.")

    # Enable the multi-threading analysis
    if not args.display:
        nthreads = 10
        ROOT.ROOT.EnableImplicitMT(nthreads)

############### CONTROLLA INPUT #######
    # Check on the validity of the inserted input.
    try:
        limits = tuple(args.limit)
    except ValueError:
        logging.error("The inserted values are not right. Check on signature "
                      "in the help's parser.")

    # Retrieve dataset from the web and start the selection on good events.
    if args.analysis==str(True):
        root_files, times = retrieve_dataset(args.file, limits[0], limits[1])
        with open(f"times.txt", "w+") as outf:
            for t in times:
                print(t, file=outf)

        # Create a RDataFrame with all the selected data.
        all_data = ROOT.RDataFrame("dimuon_w", root_files).Cache()
        N = all_data.Count().GetValue()
        logger.info(f"Total number of events that passed the cuts is {N}.")
    elif args.analysis==str(False):
        root_files = ROOT.std.vector("string")()
        root_files2 = ["dimuon_w0.root", "dimuon_w1.root", "dimuon_w2.root", "dimuon_w3.root", "dimuon_w4.root",
         "dimuon_w5.root", "dimuon_w6.root", "dimuon_w7.root", "dimuon_w8.root", "dimuon_w9.root", "dimuon_w10.root",
          "dimuon_w11.root", "dimuon_w12.root", "dimuon_w13.root", "dimuon_w14.root", "dimuon_w15.root", "dimuon_w16.root",
          "dimuon_w17.root", "dimuon_w18.root", "dimuon_w19.root", "dimuon_w20.root", "dimuon_w21.root", "dimuon_w22.root",
          "dimuon_w23.root", "dimuon_w24.root", "dimuon_w25.root", "dimuon_w26.root", "dimuon_w27.root", "dimuon_w28.root",
          "dimuon_w29.root", "dimuon_w30.root", "dimuon_w31.root", "dimuon_w32.root", "dimuon_w33.root", "dimuon_w34.root",
          "dimuon_w35.root", "dimuon_w36.root", "dimuon_w37.root", "dimuon_w38.root", "dimuon_w39.root", "dimuon_w40.root",
          "dimuon_w41.root", "dimuon_w42.root", "dimuon_w43.root", "dimuon_w44.root", "dimuon_w45.root", "dimuon_w46.root",
          "dimuon_w47.root", "dimuon_w48.root", "dimuon_w49.root", "dimuon_w50.root", "dimuon_w51.root", "dimuon_w52.root",
          "dimuon_w53.root", "dimuon_w54.root", "dimuon_w55.root", "dimuon_w56.root", "dimuon_w57.root", "dimuon_w58.root",
          "dimuon_w59.root", "dimuon_w60.root", "dimuon_w61.root", "dimuon_w62.root", "dimuon_w63.root", "dimuon_w64.root",
          "dimuon_w65.root", "dimuon_w66.root", "dimuon_w67.root", "dimuon_w68.root", "dimuon_w69.root", "dimuon_w70.root",
          "dimuon_w71.root", "dimuon_w72.root", "dimuon_w73.root", "dimuon_w74.root", "dimuon_w75.root", "dimuon_w76.root",
          "dimuon_w77.root", "dimuon_w78.root", "dimuon_w79.root", "dimuon_w80.root", "dimuon_w81.root", "dimuon_w82.root",
          "dimuon_w83.root", "dimuon_w84.root", "dimuon_w85.root", "dimuon_w86.root", "dimuon_w87.root", "dimuon_w88.root",
          "dimuon_w89.root", "dimuon_w90.root", "dimuon_w91.root", "dimuon_w92.root", "dimuon_w93.root", "dimuon_w94.root",
          "dimuon_w95.root", "dimuon_w96.root", "dimuon_w97.root", "dimuon_w98.root"]

        file_start = limits[0]
        file_stop = limits[1]
        for t in root_files2[file_start:file_stop]:
            root_files.push_back(t)

        # Create a RDataFrame with all the selected data.
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Snapshots')))
        all_data = ROOT.RDataFrame("dimuon_w", root_files).Cache()
        N = all_data.Count().GetValue()
        logger.info(f"Total number of events that passed the cuts is {N}.")
        os.chdir(os.path.dirname(os. getcwd()))

    else:
        logger.error("Invalid input for \"-a\" option. "
                    "Possibilities are: True or False.")
        exit()

    # Move the snapshots in a dedicated directory
    if args.analysis==str(True):
        start_move = time.time()
        for f in root_files:
            os.replace(f'{os.getcwd()}/{f}', f'Snapshots/{f}')

        logger.info(f"Time to move the files: {time.time()-start_move}")

    # Plot the distributions of mass and cos(theta*) in different range of
    # rapidity.
    pt_lim = (10,100)
    #pt_lim = (20,120)
    #pt_lim = (80,100)

    for i in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve the values of eta range bins (y_inf, y_sup)
        y_lim = utils.RAPIDITY_BIN[f"{i}"]

        ######### BINNAGGIO con NUMERO DI DATI?????????????????????
        bin = 40 #binning
        mass_vs_eta(all_data,y_lim,pt_lim,bin)
        cos_vs_eta(all_data,y_lim,pt_lim,bin)

    # Compute AFB (Forward_Backward Asymmetry)
    afb(all_data, pt_lim)

    # Plot AFB for different cut in pt
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))


    # Create a txt with a list of the files from which takes values.
    pt_plot = [f"{pt_lim}"]
    for p in pt_plot:
        pathlist = Path("").glob(f"*_{p}dimuon_w.txt")
        with open(f"Afblist_{p}.txt", "w+") as outf:
            for path in pathlist:
                print(path, file=outf)

    # Dimuon Afb plots
    logger.info("Plotting Afb results...")
    for p in pt_plot:
        afb_plot(f"Afblist_{p}.txt", pt_lim)
    os.chdir(os.path.dirname(os. getcwd()))

    # poi vedi come cambiare i programmi  per riprodurre la roba su più file

    # Elapsed time
    timer.Stop()
    print(f"Elapsed: {time.time()-start} seconds")
    logger.info(f"Elapsed time from the beginning is: {timer.RealTime()}(real time) seconds")
    logger.info(f"Elapsed time from the beginning is: {timer.CpuTime()}(cpu time) seconds")
