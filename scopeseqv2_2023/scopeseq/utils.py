import numpy as np
import pickle
import anndata
import loompy
from scipy.stats import zscore


def load_matrix(filename, subsample):
    gids = []
    genes = []
    matrix = []
    with open(filename) as f:
        for line in f:
            llist = line.split()
            gids.append(llist[0])
            genes.append(llist[1])
            matrix.append([int(val) for val in llist[2::]])
    gids = np.array(gids)
    ind = np.argsort(gids)
    gids = gids[ind]
    genes = np.array(genes)[ind]
    matrix = np.array(matrix)[ind]
    if subsample > 0:
        sind, smatrix = subsample_matrix(matrix, subsample)
        return gids, genes, smatrix, sind
    else:
        return gids, genes, matrix


def subsample_matrix(matrix, n):
    # subsample n cells from matrix
    ind = np.random.choice(matrix.shape[1], n, replace=False)
    return ind, matrix[:, ind]


def write_matrix(filename, gids, genes, matrix):
    with open(filename, 'w') as g:
        for gid, gene, vector in zip(gids, genes, matrix):
            st = gid+'\t'+gene+'\t'+'\t'.join([str(val) for val in vector])+'\n'
            g.write(st)
    return 0


def pickle_load(path):
    f = open(path, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj


def pickle_dump(obj, path):
    f = open(path, 'wb')
    pickle.dump(obj, f)
    f.close()
    return 0


def merge_anndata(dirs, output_dir):
    """
    merge multiple count matrix in anndata form
    :param dirs: a list of the directory path where anndata is stored
    :param output_dir: the directory path of the output folder
    :return: combined, merged anndata
            will output count_matrix.loom as the merged count matrix in the output folder;
            will output annotation as the merged cell annotation in the output folder
    """
    combined = None
    for d in dirs:
        # read in anndata count_matrix
        path = d+'count_matrix.loom'
        count_matrix = anndata.read_loom(path, obs_names='', var_names='')
        # merge anndata count_matrix
        if combined is None:
            combined = count_matrix
        else:
            combined = anndata.concat([combined, count_matrix], join='outer', merge='first')
    combined.obs = combined.obs.reset_index().drop(['index'], axis=1)
    if 'level_0' in combined.obs.columns:
        combined.obs = combined.obs.drop('level_0', axis=1)
    # output merged count matrix
    combined.write_loom(output_dir + 'count_matrix.loom')
    combined.obs.to_csv(output_dir + 'annotation.txt', sep='\t', index=False)
    return combined


def subsample_cells_anndata(ann, output_dir, levels=None):
    """
    subsample
    so every sample (defined by batch and sample in the anndata obs annotation) will have the same number of cells
    anndata.obs must have columns 'batch' 'sample' and 'CellID'
    :param ann: count matrix in anndata form
    :param output_dir: output directory
    :return: ann_sub, subsampled anndata
            will output count_matrix_sub.loom as subsampled anndata to the output folder
            will output annotation_sub.txt as subsamples anndata cell anntation to the output folder
    """
    if not levels:
        levels = ['batch', 'sample']

    ann.obs['subsample'] = None
    # get the number of cells for each sample
    summary = ann.obs.groupby(levels).agg('count')['CellID'].reset_index()
    # get the number of cells for the smallest sample, ncells
    ncells = summary['CellID'].min()
    for i in range(summary.shape[0]):
        data = ann.obs['CellID'].notna()
        for level in levels:
            data = data & (ann.obs[level] == summary.loc[i, level])
            index = ann.obs.index[data]
        # random select ncells of each sample, store subsample info in the 'subsample' column
        n = summary.loc[i, 'CellID']
        rmask = np.zeros(n, dtype=int)
        rmask[:ncells] = 1
        np.random.shuffle(rmask)
        rmask = rmask.astype(int)
        ann.obs.loc[index, 'subsample'] = rmask
    # subsample cells
    ann_sub = ann[ann.obs['subsample'] == 1]

    # output subsampled anndata
    ann_sub.write_loom(output_dir + 'count_matrix_sub.loom')
    ann_sub.obs.to_csv(output_dir + 'annotation_sub.txt', sep='\t', index=False)

    return ann_sub


def filter_cells_anndata(ann, marker_file):
    marker = pd.read_csv(marker_file,sep='\t',header=None)
    ann_T = ann.copy().T
    ann_filter = ann_T[np.isin(ann_T.obs['Accession'],marker[0])].copy().T
    ann_filter = ann_filter[ann_filter.X.sum(axis=1)>0]
    ann_single = ann_filter[ann_filter.obs['sample']!='Doublet']
    ann_single = ann[np.isin(ann.obs.index.values,ann_single.obs.index.astype('int').values)]
    return ann_single


def load_matrix_toloom(matrix_file, loom_file, pg_file=None, umap_file=None, subsample=0):
    gid, genes, matrix = load_matrix(filename=matrix_file, subsample=subsample)
    var = pd.DataFrame({'Accession': gids, 'Gene': genes})
    var.index = var.index.astype('str')

    if pg_file:
        pgs = pd.read_csv(pg_file,sep='\t',header=None)
        obs = pd.DataFrame({'pgs': pgs[0].to_list()})
    else:
        obs=None

    adata = anndata.AnnData(X=matrix.T, obs=obs, var=var)

    if umap_file:
        umap_emb = pd.read_csv(umap_file,sep='\t',header=None)
        adata.obsm['umap_emb'] = umap_emb.to_numpy()

    adata.write_loom(loom_file, write_obsm_varm=True)


def loom_addcluster(loom_file, loom_out_file, pg_file=None, umap_file=None):
    adata = anndata.read_loom(loom_file, obs_names='', var_names='')

    if pg_file:
        pgs = pd.read_csv(pg_file,sep='\t',header=None)
        adata.obs['pgs']= pgs[0].to_list()

    if umap_file:
        umap_emb = pd.read_csv(umap_file,sep='\t',header=None)
        adata.obsm['umap_emb'] = umap_emb.to_numpy()

    adata.write_loom(loom_out_file, write_obsm_varm=True)


def cptt(matrix):
    norm = matrix.sum(axis=0)
    return np.log2(matrix/norm*10000.0+1.0)


def loom_norm(loom_file, norm_loom_file, sparse=True):
    adata = anndata.read_loom(loom_file, obs_names='', var_names='')

    if sparse:
        matrix_norm = cptt(adata.X.todense().T)
    else:
        matrix_norm = cptt(adata.X.T)
    adata_norm = anndata.AnnData(X=matrix_norm.T, obs=adata.obs, var=adata.var)

    if adata.obsm:
        adata_norm.obsm = adata.obsm
    adata_norm.write_loom(norm_loom_file, write_obsm_varm=True)


def load_genesets_fromtxt(genesets_file):
    df_gs=pd.read_csv(genesets_file,sep='\t',header=0)
    df_gs=df_gs.dropna(how='all')
    genesets_dict = {}
    for key in df_gs.columns:
        genesets_dict[key] = df_gs[key].dropna().tolist()

    return genesets_dict


def map(x, ref_before, ref_after):
    mapping = pd.DataFrame(ref_after, index=ref_before)
    mapped = mapping.loc[x].values
    return mapped


def marker_expression(adata_norm, marker, group, sparse=True):
    if sparse:
        matrix_norm = adata_norm.X.todense().T
    else:
        matrix_norm = adata_norm.X.T

    matrix_marker = pd.DataFrame()
    for i in marker:
        if i in adata_norm.var['Gene'].values:
            index = np.where(adata_norm.var['Gene']==i)[0][0]
            matrix_marker[i] = matrix_norm[index,:].tolist()[0]
    matrix_marker_tmp = matrix_marker.copy()
    matrix_marker_tmp[group] = adata_norm.obs[group].values
    matrix_marker_mean = matrix_marker_tmp.groupby([group]).mean()
    matrix_marker_mean_zscore = matrix_marker_mean.apply(zscore,axis=0)

    return matrix_marker, matrix_marker_mean, matrix_marker_mean_zscore


def geneset_expression(adata_norm, geneset, group, sparse=True):
    if sparse:
        matrix_norm = np.array(adata_norm.X.todense().T)
    else:
        matrix_norm = adata_norm.X.T

    matrix_geneset = pd.DataFrame()
    for key in geneset.keys():
        matrix_marker = pd.DataFrame()
        for i in geneset[key]:
            if i in adata_norm.var['Gene'].values:
                index = np.where(adata_norm.var['Gene']==i)[0][0]
                matrix_marker[i] = matrix_norm[index,:].tolist()
        matrix_geneset[key] = matrix_marker.mean(axis=1)
    matrix_geneset_tmp = matrix_geneset.copy()
    matrix_geneset_tmp[group] = adata_norm.obs[group].values
    matrix_geneset_mean = matrix_geneset_tmp.groupby([group]).mean()
    matrix_geneset_mean_zscore = matrix_geneset_mean.apply(zscore,axis=0)

    return matrix_geneset, matrix_geneset_mean, matrix_geneset_mean_zscore


def plt_marker_heatmap(pdf_outfile, matrix_marker_mean, matrix_marker_mean_zscore):
    with PdfPages(pdf_outfile) as pdf:
        sns.heatmap(matrix_marker_mean.T, cmap='Reds',cbar_kws={'label':'Average normalized expression'})
        plt.xlabel('')
        plt.ylabel('')
        pdf.savefig(bbox_inches='tight')
        plt.close()

        sns.heatmap(matrix_marker_mean_zscore.T, cmap='coolwarm',cbar_kws={'label':'zscore(Average normalized expression)'},
                vmin=-2.5, vmax=2.5)
        plt.xlabel('')
        plt.ylabel('')
        pdf.savefig(bbox_inches='tight')
        plt.close()


def plt_marker_umap(pdf_outfile, adata_norm, matrix_marker):
    with PdfPages(pdf_outfile) as pdf:
        for i in matrix_marker.columns:
            plt_map(pdf, x=adata_norm.obsm['umap_emb'], name=i,
                    color=matrix_marker[i].values, cmap='Reds', t='float', s=5)


def subsample_adata(adata, group=None, n=0):
    if group:
        group_set = set(adata.obs[group])
        if n>0:
            n_sub=min(n,adata.obs[group].value_counts().min())
        else:
            n_sub=adata.obs[group].value_counts().min()
        index_sub=[]
        for g in group_set:
            index_group = adata[adata.obs[group]==g].obs.index.values
            index_group_sub = np.random.choice(index_group,n_sub, replace=False).tolist()
            index_sub = index_sub+index_group_sub
        adata_sub = adata[index_sub]
    else:
        if n>0:
            n_sub=min(n,adata.shape[0])
        else:
            n_sub=adata.shape[0]
        index_group = adata.obs.index.values
        index_sub = np.random.choice(index_group,n_sub, replace=False).tolist()
        adata_sub = adata[index_sub]
    print(n_sub)
    return adata_sub

