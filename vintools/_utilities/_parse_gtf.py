
def _note_header_lines(gtf, line, header_key="header"):

    """

    Parameters:
    -----------
    gtf
        type: dict

    line
        type: str (iterable)

    header_key
        default: header
        type: str

    Returns:
    --------

    Notes:
    ------

    """

    return gtf[header_key].append(line.strip("\n"))


def _format_gene_info_string(feature_info_string, field):

    """

    Parameters:
    -----------
    feature_info_string

    field

    Returns:
    --------


    Notes:
    ------
    """

    gene_annotation = (
        feature_info_string.split(field)[
            1
        ]  # split at the label; take the info right after
        .replace('"', "")
        .replace(" ", "")
        .split(";")[0]
    )

    return gene_annotation


def _get_gene_annotation(splitline):

    """
    
    Parameters:
    -----------
    """
    feature_info_string = splitline[8]

    gene_id = _format_gene_info_string(feature_info_string, "gene_id")
    gene_name = _format_gene_info_string(feature_info_string, "gene_name")
    gene_type = _format_gene_info_string(feature_info_string, "gene_type")

    return gene_id, gene_name, gene_type


def _parse_gtf_as_dictionary(gtf_file_path, feature="gene"):

    """

    Parameters:
    -----------
    gtf_file_path
        path to .gtf file.
        type: str

    feature
        default: 'gene'
        type: str

    Returns:
    --------
    gtf_dict

    Notes:
    -----_
    (1)  The only feature implemeted so far is 'gene'
    """

    gtf_dict = {}
    with open(gtf_file_path) as file:
        for n, line in enumerate(file):
            if line.startswith("##"):
                _note_header_lines(gtf_dict, line)
            else:
                splitline = line.split("\t")
                chromosome = splitline[0]

            if chromosome not in gtf_dict.keys():
                gtf[chromosome] = {}
            if splitline[2] == feature:
                gene_id, gene_name, gene_type = _get_gene_annotation(splitline)
                gtf[chromosome][gene_name] = {}
                gtf[chromosome][gene_name]["gene_id"] = gene_id
                gtf[chromosome][gene_name]["lower_bound"] = splitline[3]
                gtf[chromosome][gene_name]["upper_bound"] = splitline[4]
                gtf[chromosome][gene_name]["strand"] = splitline[6]
                gtf[chromosome][gene_name]["gene_type"] = gene_type

    return gtf_dict
