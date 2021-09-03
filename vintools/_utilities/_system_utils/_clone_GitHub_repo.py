import os

def _clone_GitHub_repo(repository, destination=False):
    
    """Clone a GitHub repository. """
    
    if not destination:
        destination=os.path.join(os.getcwd(), os.path.basename(repository).strip('.git'))
    
    
    executable = "git clone {} {}".format(repository, destination)
    os.system(executable)