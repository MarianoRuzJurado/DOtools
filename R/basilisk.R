DOtoolsEnv <- basilisk::BasiliskEnvironment(
    envname = "DOtools_env",
    pkgname = "DOtools",
    packages = "python=3.11",
    pip = c(
        "scvi-tools==1.3.0",
        "celltypist==1.6.3",
        "scanpro==0.3.2",
        "scipy==1.15.3",
        "scib==1.1.7"
    )
)

CellBenEnv <- basilisk::BasiliskEnvironment(
    envname = "CellBen_env",
    pkgname = "DOtools",
    packages = "python=3.7",
    pip = c(
        "cellbender==0.3.0",
        "lxml_html_clean==0.4.2"
    )
)
