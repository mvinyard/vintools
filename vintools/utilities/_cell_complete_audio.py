
from IPython.display import display, Audio
import os

class InvisibleAudio(Audio):
    '''
    IPython.display.Audio, but without it being displayed. It still takes space in the output 
    layout tough.
    '''
    def _repr_html_(self):
        audio = super()._repr_html_()
        audio = audio.replace('<audio', f'<audio onended="this.parentNode.removeChild(this)"')
        return f'<div style="display:none">{audio}</div>'
    
def _ipy_cell_complete_audio(download_destination="."):
    
    """Plays a sound when a jupyter cell completes."""
    
    sound_file = os.path.join(download_destination, "business_done.wav")
    display(InvisibleAudio(sound_file, autoplay=True))
