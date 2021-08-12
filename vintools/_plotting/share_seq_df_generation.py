
def shareseq_palette_df():
    
    """Defines the dataframe for colors pertaining to the celltypes of the SHARE-seq dataset."""

    import pandas as pd

    d = {
        "celltypes": [
            "alpha-high-CD34+ Bulge",
            "alpha-low-CD34+ Bulge",
            "Isthmus",
            "K6+ Bulge/Companion Layer",
            "TAC-1",
            "TAC-2",
            "IRS",
            "Medulla",
            "Hair Shaft - Cuticle/Cortex",
            "ORS",
            "Basal",
            "Spinous",
            "Granular",
            "Infundibulum",
            "Endothelial",
            "Dermal Fibroblast",
            "Dermal Sheath",
            "Dermal Papilla",
            "Macrophage DC",
            "Melanocyte",
            "Sebaceous Gland",
            "Schwann Cell",
            "Mixed"
        ],
        "colors": [
            "#0F532C",
            "#1D8342",
            "#42B649",
            "#69C4A5",
            "#FDD82F",
            "#FCAE1F",
            "#F57A03",
            "#EF2D1A",
            "#A10109",
            "#660A17",
            "#96CBE8",
            "#149FD7",
            "#0765AC",
            "#427FC2",
            "#8984C0",
            "#6C4CA1",
            "#98479C",
            "#A80864",
            "#491786",
            "#ECE819",
            "#FEC75C",
            "#E84C9D",
            "#161616"
        ],
    }

    shareseq_df = pd.DataFrame(d)
    
    return shareseq_df
