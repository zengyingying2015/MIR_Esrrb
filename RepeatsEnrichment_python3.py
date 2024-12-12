#!/home/jllab/miniconda3/bin/ python
#--coding:utf-8 --
"""
RepeatsEnrichment_bg.py
Parse the sub-family, binomial test p-values added, combined p-values added.
Take kgg as input file, first column is TErepeatsname, second column is sample/Celltype/Factor. without header.
chr1|3000001|3000156	ESC


TE reference file format:
chr	start	end	strand	repname	repclass	repfamily	repeatsname
chr1	3000001	3000156	-	L1_Mur2	LINE	L1	chr1|3000001|3000156
chr1	3000237	3000733	-	L1_Mur2	LINE	L1	chr1|3000237|3000733
chr1	3000733	3000766	+	(TTTG)n	Simple_repeat	Simple_repeat	chr1|3000733|3000766

"""


#systematic library
import glob,os
from collections import Counter

#3rd library
#plot setting
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
#import pylab
sns.set_style("whitegrid")
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
#computating and stat setting.
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues as cps
import Orange
from orangecontrib.bio.utils.stats import Hypergeometric as hy
from orangecontrib.bio.utils.stats import Binomial as bin
from orangecontrib.bio.utils.stats import FDR as fdr
#from joblib import Parallel, delayed



def getRepFamily(rs, reps):
    """
    Get family count for each input file.
    """
    fams = list(reps[rs].values)
    cs = dict(Counter(fams))
    return cs


def repST(fe, fin, pre, reps):
    fout = pre + "_st.txt"
    if os.path.exists(fout):
        print( "%s has been generated,return" % fout)
        return
    fg = getRepFamily(fe, reps)
    bg = getRepFamily(fin, reps)
    N = sum(bg.values())
    n = sum(fg.values())
    data = {}
    fgkeys=[key for key in fg.keys() if not pd.isnull(key)]
    for key in fgkeys:
        m = bg[key]
        k = fg[key]
        h = hy()
        b = bin()
        hp = h.p_value(k, N, m, n)
        bp = b.p_value(k, N, m, n)
        es = float(k) / float(m) / float(n) * float(N)
        cp = cps([hp, bp], method="stouffer")[1]
        data[key] = {
            "bg": m,
            "fg": k,
            "hy_p": hp,
            "bin_p": bp,
            "es": es,
            "combinedP": cp
        }
    data = pd.DataFrame(data).T
    ps = data["combinedP"].values
    qs = fdr(ps)
    data["FDR"] = qs
    data = data.sort_values("FDR")
    data.to_csv(fout, sep="\t", index_label="repFamily")



def preMat(fs=glob.glob("c*st.txt"), qcut=1e-10,cut=1):
    reps = set()
    mats = {}
    for f in fs:
        pre = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        mats[pre] = mat
        qs = mat["FDR"]
        qs = qs[qs < qcut]
        rs = set(qs.index)
        reps.update(rs)
    reps = list(reps)
    mat_qs = {}
    mat_es = {}
    for pre, mat in mats.items():
        mat_qs[pre] = mat.loc[reps, "FDR"]
        mat_es[pre] = mat.loc[reps, "es"]
    mat_qs = pd.DataFrame(mat_qs)
    mat_es = pd.DataFrame(mat_es)
    mat_qs[mat_qs < 1e-100] = 1e-100
    mat_qs[mat_qs > qcut] = 1
    a = mat_qs.sum(axis=1)
    a.sort_values(ascending=True)
    mat_qs = 0.0 - np.log10(mat_qs)
    mat_qs = mat_qs.loc[a.index,:]
    #filter only one situation
    #mat_qs = filterRs(mat_qs,cut=cut)
    #sort all samples
    #mat_qs = sortCs(mat_qs)
    mat_es = mat_es.loc[mat_qs.index,mat_qs.columns]
    mat_qs.to_csv("Module_FDR.csv", sep="\t", index_label="rep")
    mat_es.to_csv("Module_EnrichmentScore.csv", sep="\t", index_label="rep")


def main():

    repf = "/home/jllab/Analysis/TEProject/TEs_mm9.txt"
    reps = pd.read_table(repf, index_col=7)
    reps = reps["repfamily"]
    selkgg = "E2-enhancers_TEs.kgg"
    selkgs = pd.read_csv(selkgg,sep="\t")
    bg = list(reps.index)
    for k in set(selkgs.iloc[:,1]):
        rs = list(selkgs[selkgs.iloc[:,1]==k].iloc[:,0])
        repST(rs, bg, "c_%s"%k, reps)
    #####qcut value can be changed accordingly.
    preMat(fs=glob.glob("c*st.txt"), qcut=1e-5,cut=1)



if __name__ == "__main__":
    main()
