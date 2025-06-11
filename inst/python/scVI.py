import scvi
import anndata as ad
def run_scvi(
    adata: ad.AnnData,
    batch_key: str,
    layer_counts: str = "counts",
    layer_logcounts: str = "logcounts",
    categorical_covariates: list = None,
    continuos_covariates: list = None,
    n_hidden: int = 128,
    n_latent: int = 30,
    n_layers: int = 3,
    dispersion: str = "gene-batch",
    gene_likelihood: str = "zinb",
    get_model: bool = False,
    **kwargs
) -> None:
    """Run scVI.

    Run scVI to integrate sc/snRNA more information on
    `scvi-tools <https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html>`_.

    :param adata: annotated dt matrix.
    :param batch_key: `.obs` column with batch information.
    :param layer_counts: layer with counts. Raw counts are required.
    :param layer_logcounts: layer with log-counts. Log-counts required for calculation of HVG.
    :param categorical_covariates: `.obs` column names with categorical covariates for scVI inference.
    :param continuos_covariates: `.obs` column names with continuous covariates for scVI inference.
    :param n_hidden: number of hidden layers.
    :param n_latent: dimensions of the latent space.
    :param n_layers: number of layers.
    :param dispersion: dispersion mode for scVI.
    :param gene_likelihood: gene likelihood.
    :param get_model: return the trained model.
    :param kwargs: additional arguments for `scvi.model.SCVI`.
    :return: None or the model, the latent space is saved in the anndata under X_scVI.
    """

    print("Run scVI")
    assert layer_logcounts in adata.layers, "logcounts layer not in anndata"
    assert layer_counts in adata.layers, "counts layer not in anndata"

    # Set-up anndata and model
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer_counts,
        batch_key=batch_key,
        continuous_covariate_keys=continuos_covariates,
        categorical_covariate_keys=categorical_covariates,
    )

    model_scvi = scvi.model.SCVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dispersion=dispersion,
        gene_likelihood=gene_likelihood,
        **kwargs,
    )

    model_scvi.view_anndata_setup()
    model_scvi.train()  # Train
    adata.obsm["X_scVI"] = model_scvi.get_latent_representation()

    if get_model:
        return model_scvi
    else:
        del model_scvi
        return None
