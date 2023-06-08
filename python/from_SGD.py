from pygenome import sg
import re


def add_sgd_data(sc_description):
    pattern = r"(?P<gene_name>[a-zA-Z]{3}\d{1,2})"
    try:
        sc_gene_name = re.search(pattern=pattern, string=sc_description)
    except KeyError:
        return None

    if sc_gene_name is None:
        return None

    try:
        sgd_res = sg.stdgene[sc_gene_name]
        sgd_description = sgd_res.short_description
    except KeyError:
        sgd_description = None

    return sgd_description
