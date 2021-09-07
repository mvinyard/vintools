# made for out of package
# in package, download is not necessary

from IPython.display import display, Audio
import os

def gsutil_get(fetch, destination):
    
    """"""
    
    function = "gsutil -m cp gs://"
    gsutil = function + fetch + " " + destination
    
    os.system(gsutil)
    
class InvisibleAudio(Audio):
    '''
    IPython.display.Audio, but without it being displayed. It still takes space in the output 
    layout tough.
    '''
    def _repr_html_(self):
        audio = super()._repr_html_()
        audio = audio.replace('<audio', f'<audio onended="this.parentNode.removeChild(this)"')
        return f'<div style="display:none">{audio}</div>'
    
def ipy_cell_complete_audio(download_destination="."):
    
    """Plays a sound when a jupyter cell completes."""
    
    sound_file = os.path.join(download_destination, "business_done.wav")
    
    if os.path.exists(sound_file)==False:
        fetch = "vintools/business_done.wav"
        gsutil_get("vintools/business_done.wav", download_destination)
    
    
    display(InvisibleAudio(sound_file, autoplay=True))