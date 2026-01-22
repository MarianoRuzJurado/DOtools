import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.simplefilter("ignore", FutureWarning)
import scanpy as sc
import pandas as pd
import scib
import anndata as ad
import platform


#adata = sc.read_h5ad("/Users/mariano/Desktop/Work_PhD/ProjectOrtho/R_analysis/Dotools/sce_data.h5ad")
def eval_integration(
        adata: ad.AnnData,
        label_key: str,
        batch_key: str,
        embed: str,
        use_rep: str,
        covariate: str,
        type_: str = "embed",
        n_comps: int = 30,
        n_cores: int = 10,
        scale: str = True,
        verbose: bool = False,
        return_all_silhouette: bool = False,
        ):
    os_system = platform.system()
    adata.obsm[embed] = adata.obsm[embed].to_numpy()
    adata.obsm["X_pca"] = adata.obsm["PCA"]

    # scib package metric calculation for principal component regression
    pcr_comparison = scib.metrics.pcr_comparison(adata_pre=adata,
                                                 adata_post=adata,
                                                 covariate=covariate,
                                                 n_comps=n_comps,
                                                 embed=embed,
                                                 scale=scale,
                                                 verbose=verbose)
    # and silhouette score per batch
    silhouetteBatch = scib.metrics.silhouette_batch(adata=adata,
                                                    batch_key=batch_key,
                                                    label_key=label_key,
                                                    embed=embed,
                                                    return_all=return_all_silhouette,
                                                    scale=scale,
                                                    verbose=verbose)
    # graph connectivity
    sc.pp.neighbors(adata, use_rep=embed)
    graph_connectivity = scib.metrics.graph_connectivity(adata, label_key=label_key)

    #LISI scoring if plattform is not MacOS, should work on Ubuntu
    if(os_system == "Linux"):
        LISI_score = scib.metrics.ilisi_graph(adata,
                                              batch_key=batch_key,
                                              type_=type_,
                                              use_rep=use_rep,
                                              scale=scale,
                                              n_cores=n_cores,
                                              verbose=verbose)
    else:
        raise ValueError(
            f"LISI Score is only supported on Linux!"
            f"Detected platform: {os_system}"
        )

    # Build Dataframe with results for OrthoIntegrate
    output_df = pd.DataFrame()
    output_df.loc['graph_connectivity', embed] = graph_connectivity
    output_df.loc[
        'pcr_comparison', embed] = pcr_comparison  # lower scores mean less impact of each sample to our pca analysis, so it become more robust to individual samples
    output_df.loc['silhouette_batch', embed] = silhouetteBatch

    if (os_system == "Linux"):
        output_df.loc['LISI', embed] = LISI_score

    return output_df