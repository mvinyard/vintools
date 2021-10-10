
from ..._utilities._system_utils._flexible_mkdir import _flexible_mkdir

import os, datetime

def _name_compiled_10x_web_summary_dir(title):
    if title:
        title = title + "."
    else:
        title = ""
    when = datetime.datetime.now().strftime("%d%b%Y")

    return "{}compiled_10x_web_summaries.{}".format(title, when)

def _make_compiled_10x_web_summary_dir(path_to_results, title):
    
    web_sum_dir = os.path.join(path_to_results.strip("*"), _name_compiled_10x_web_summary_dir(title))
    _flexible_mkdir(web_sum_dir)
    
    return web_sum_dir

def _cp_web_summaries(ResultsDict, gcp_bucket, compiled_web_summary_title=False):
    
    """"""
    
    web_sum_dir = _make_compiled_10x_web_summary_dir(ResultsDict['run_path'], title=compiled_web_summary_title)
    for sample in ResultsDict['samples'].keys():
        web_summary = ResultsDict['samples'][sample]['outs']['web_summary.html']
        web_summary_filename = os.path.basename(web_summary)
        web_summary_renamed = os.path.join(web_sum_dir, ".".join([sample, web_summary_filename]))
        os.system("cp {} {}".format(web_summary, web_summary_renamed))
    
    os.system("gsutil -m cp -r {} gs://{}".format(web_sum_dir, gcp_bucket))