import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import ranksums, mannwhitneyu, pearsonr, zscore

import phenograph
import umap
# import dmaps
from scipy import sparse

from scopeseq.utils import subsample_adata

"""
umap, phenograph
"""


def run_phenograph(matrix, metric='correlation', k=20):
    pgs, graph, q = phenograph.cluster(matrix, k=k, primary_metric=metric)

    n = len(set(pgs))

    if n < 11:
        pgcolors = mpl.cm.tab10.colors
    elif n < 21:
        pgcolors = mpl.cm.tab20.colors
    else:
        pgcolors = [name for name, hex in mpl.colors.cnames.items()]
        pgcolors.reverse()
    pgcolors = pgcolors[:n]

    return pgs, n, pgcolors


def run_umap(matrix, metric='correlation', k=20):
    umap_model = umap.UMAP(metric=metric, n_neighbors=k).fit(matrix)
    return umap_model.embedding_


"""
aneoploidy
"""


def filt_chr_sets(gids, matrix, chr_sets, thresh, exclude):
    new_chr_sets = [[] for c in chr_sets]
    for i in range(len(gids)):
        gid = gids[i]
        if sparse.issparse(matrix):
            p = (matrix.getrow(i).count_nonzero() > thresh)
        else:
            p = (np.count_nonzero(matrix[i]) > thresh)
        if p and gid not in exclude:
            for j in range(len(chr_sets)):
                if gid in chr_sets[j]:
                    new_chr_sets[j].append(gid)
                    break
    new_chr_sets = [set(c) for c in new_chr_sets]
    return new_chr_sets


def cptt(matrix):
    norm = matrix.sum(axis=0)
    return np.log2(matrix/norm*10000.0+1.0)


def get_ploidy(chr_sets, gids, matrix):
    cptt_matrix = np.array(cptt(matrix))
    cmatrix = []
    for c in chr_sets:
        cvec = np.array([0.0 for cell in range(len(cptt_matrix[0]))])
        j = 0
        for i in range(len(gids)):
            gid = gids[i]
            if gid in c:
                pt = cptt_matrix[i, :]
                cvec += pt
                j += 1
        cmatrix.append(cvec/float(j))
    return np.array(cmatrix)


def run_ploidy(gmt_file, gids, matrix, chr_use=None, amplify=None):
    if chr_use is None:
        chr_use = [6, 9]
    if amplify is None:
        amplify = [True, False]

    with open(gmt_file) as f:
        chr_sets = [set(line.split()[2::]) for line in f]

    exclude = {'ENSG00000223865.10', 'ENSG00000230795.3', 'ENSG00000206503.11', 'ENSG00000196126.10',
                'ENSG00000228078.1', 'ENSG00000204622.11', 'ENSG00000204632.11', 'ENSG00000223534.1',
                'ENSG00000224557.7', 'ENSG00000235290.1', 'ENSG00000235301.1', 'ENSG00000204642.13',
                'ENSG00000243753.5', 'ENSG00000231461.1', 'ENSG00000231389.7', 'ENSG00000261548.1',
                'ENSG00000204592.8', 'ENSG00000198502.5', 'ENSG00000204287.13', 'ENSG00000237398.1',
                'ENSG00000204525.15', 'ENSG00000232629.8', 'ENSG00000241106.6', 'ENSG00000234745.9',
                'ENSG00000224372.1', 'ENSG00000225851.1', 'ENSG00000231130.1', 'ENSG00000242574.8',
                'ENSG00000214922.9', 'ENSG00000196735.11', 'ENSG00000206341.7', 'ENSG00000229391.7',
                'ENSG00000179344.16', 'ENSG00000237541.3', 'ENSG00000226030.1', 'ENSG00000204252.12',
                'ENSG00000181126.13', 'ENSG00000196301.3', 'ENSG00000204257.14'}

    filt_chr = filt_chr_sets(gids, matrix, chr_sets, 100, exclude)

    ploidy = get_ploidy(filt_chr, gids, matrix)

    zploidy_score = None
    for i in range(len(chr_use)):
        if amplify[i]:
            score = zscore(ploidy.T[:, chr_use[i]])
        else:
            score = zscore(ploidy.T[:, chr_use[i]])*(-1)
        if zploidy_score is None:
            zploidy_score = score
        else:
            zploidy_score = zploidy_score + score

    return ploidy, zploidy_score


"""
diffusion map

def diffusion_map():
    dist = np.empty([data.shape[0], data.shape[0]])
    for i in range(data.shape[0]):
        for j in range(data.shape[0]):
            dist[i, j] = pearsonr(data[i, :], data[j, :])[0]

    # Compute top three eigenvectors.
    # Here we assume a good value for the kernel bandwidth is known.
    dmap = dmaps.DiffusionMap(dist)
    dmap.set_kernel_bandwidth(0.5)
    dmap.compute(3)

    # Plot result. Scale by top eigenvector.
    v = dmap.get_eigenvectors()

    diffusion_emb = pd.DataFrame([v[:, 1]/v[:, 0], v[:, 2]/v[:, 0]]).T

    return diffusion_emb
"""


def neftel_subtype(adata_norm, genesets_infile, group=None, b=100):
    adata_norm_subpgs = subsample_adata(adata_norm, group=group, n=0)
    df_nmatrix_subpgs = pd.DataFrame(adata_norm_subpgs.X.todense().T)
    genes_subpgs = adata_norm_subpgs.var['Gene']
    ind = (df_nmatrix_subpgs > 0).any(axis=1).values
    gexp_subpgs = genes_subpgs[ind]
    df_nmatrix_subpgs = df_nmatrix_subpgs[ind]
    df_nmatrix_subpgs.index = gexp_subpgs

    # aggregate gene expression to 30 bin
    agg_exp_subpgs = np.sum(df_nmatrix_subpgs, axis=1).tolist()
    df_gexp_subpgs = pd.DataFrame({'Gene': gexp_subpgs, 'Agg_Exp': agg_exp_subpgs})
    df_gexp_subpgs = df_gexp_subpgs.set_index('Gene')
    exp_bin_subpgs = pd.qcut(df_gexp_subpgs['Agg_Exp'], b, labels=False,
                             duplicates='drop')  # equal element/sample in each bin
    df_gexp_subpgs['bin'] = exp_bin_subpgs

    # read genesets
    df_gs = pd.read_csv(genesets_infile, sep='\t', header=0)
    df_gs = df_gs.dropna(how='all')
    # genes in genesets not expressed
    gs_exclude = []
    gs_include = []
    ctrl_mask = []
    for i in range(df_gs.shape[1]):
        gs_include_temp = []
        for j in range(df_gs.shape[0]):
            if df_gs.loc[j][df_gs.columns[i]] not in gexp_subpgs.values:
                gs_exclude.append(df_gs.loc[j][df_gs.columns[i]])
            else:
                gs_include_temp.append(df_gs.loc[j][df_gs.columns[i]])
                ctrl_mask.append(df_gs.loc[j][df_gs.columns[i]])
        gs_include.append(gs_include_temp)
    print(set(gs_exclude), len(gs_include), len(set(ctrl_mask)))

    genesets_dict = dict(zip(df_gs.columns.tolist(), gs_include))

    # exclude genes in geneset from df_gexp, use as control genes
    mask = df_gexp_subpgs.index.isin(ctrl_mask)
    df_ex_subpgs = df_gexp_subpgs[~mask]

    # define control genes for each gene set
    ctrlgene_allsets = []
    for key in df_gs.columns.tolist():
        genes_gs = genesets_dict[key]
        df_set = df_gexp_subpgs.loc[genes_gs]
        bin_list = df_set['bin'].tolist()
        # random choice 100 genes from df_ex with the same bin of each gene in df_set
        ctrl_genes = []
        for i in range(len(bin_list)):
            gbin = bin_list[i]
            df_binctrl = df_ex_subpgs[df_ex_subpgs['bin'] == gbin]
            ctrl_r_temp = list(np.random.choice(np.array(df_binctrl.index), 100, replace=False))
            # ctrl_genes.append(ctrl_r_temp)
            ctrl_genes = ctrl_genes + ctrl_r_temp
        print(len(bin_list), len(ctrl_genes))
        ctrlgene_allsets.append(ctrl_genes)
    len(ctrlgene_allsets)

    ctrlsets_dict = dict(zip(df_gs.columns.tolist(), ctrlgene_allsets))

    # normalize and filter out non-zero expression genes
    df_nmatrix = pd.DataFrame(adata_norm.X.todense().T)
    genes = adata_norm.var['Gene']
    ind = (df_nmatrix > 0).any(axis=1).values
    gexp = genes[ind]
    df_nmatrix = df_nmatrix[ind]
    df_nmatrix.index = gexp

    # calculate SC for each cell from nmatrix_exp
    SC_all = pd.DataFrame()
    for key in df_gs.columns.tolist():
        genes_gs = genesets_dict[key]
        genes_ctrl = ctrlsets_dict[key]
        df_nmatrix_gs = df_nmatrix.loc[genes_gs]
        Er_gs = df_nmatrix_gs.mean(axis=0).tolist()  # The list of average[Er(Gj,i)] order by cells
        df_nmatrix_ctrls = df_nmatrix.loc[genes_ctrl]
        Er_ctrl = df_nmatrix_ctrls.mean(axis=0).tolist()  # The list of average[Er(Gjctrl,i)] order by cells
        SC_j = [x - y for x, y in zip(Er_gs, Er_ctrl)]
        df_SC_j = pd.DataFrame({key: SC_j})
        SC_all = pd.concat([SC_all, df_SC_j], axis=1)

    # Two-dimensional representaion
    # OPC/NPC vs AC/MES
    x = []
    y = []
    for i in range(SC_all.shape[0]):
        D = max(SC_all.iloc[i]['OPC'], SC_all.iloc[i]['NPC1'], SC_all.iloc[i]['NPC2']) - max(SC_all.iloc[i]['AC'],
                                                                                             SC_all.iloc[i]['MES1'],
                                                                                             SC_all.iloc[i]['MES2'])

        if D > 0:
            y.append(D)
            SC_on = max(SC_all.iloc[i]['NPC1'], SC_all.iloc[i]['NPC2']) - SC_all.iloc[i]['OPC']
            x.append(SC_on)

        else:
            y.append(D)
            SC_am = max(SC_all.iloc[i]['MES1'], SC_all.iloc[i]['MES2']) - SC_all.iloc[i]['AC']
            x.append(SC_am)

    hetero_2D = pd.DataFrame({'x': x, 'y': y})

    return hetero_2D, SC_all

