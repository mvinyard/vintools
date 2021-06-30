def vin_colors(random=False, subset=None, separable=False):
    
    """
    Color presets:
    --------------
    Red, Orange, Gold, Peach, Yellow, Green, Turqoise, Moss, Navy, Blue
    """
    
    import numpy as np

    color_presets = ["#E5514C", "#E37A3F", "#EA9B3F", "#E98A57","#F0C864","#9ABD76","#60A78D","#5D8F8D","#5D738D","#417B9E"]
    
    if separable==True:
        color_presets = ["#E5514C", "#9ABD76","#E37A3F", "#5D738D","#60A78D","#EA9B3F", "#5D8F8D","#E98A57","#417B9E","#F0C864"]
    
    if random == True:
        try:
            selection = np.random.choice(range(0,len(color_presets)), subset, replace=False).astype(int)
            color_presets = np.array(color_presets)[selection]
        except:
            print("Subset must be specified to select colors at random")

    return color_presets


def shareseq_colors():
    
    from .share_seq_df_generation import shareseq_palette_df
    
    df = shareseq_palette_df()
    
    return df
    