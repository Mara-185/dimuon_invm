import ROOT
import os
import utils

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
