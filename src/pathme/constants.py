# -*- coding: utf-8 -*-

"""Constants for PathMe.

This module contains all the string constants used in PathMe.
"""

import os

from bio2bel.utils import get_connection
from pybel.dsl import abundance, complex_abundance, protein, rna

MODULE_NAME = 'pathme'
DEFAULT_PATHME_DIR = os.path.join(os.path.expanduser('~'), '.pathme')
PATHME_DIR = os.environ.get('PATHME_DIRECTORY', DEFAULT_PATHME_DIR)

KEGG_PATHWAYS_URL = 'http://rest.kegg.jp/link/pathway/'


def get_data_dir() -> str:
    """Ensure the appropriate PathMe data directory exists for the given module, then returns the file path."""
    os.makedirs(PATHME_DIR, exist_ok=True)
    return PATHME_DIR


DATA_DIR = get_data_dir()
DEFAULT_CACHE_CONNECTION = get_connection("pathme.constants")

# Databases contained in PathMe
#: KEGG
KEGG = 'kegg'
KEGG_DIR = os.path.join(DATA_DIR, KEGG)
KEGG_BEL = os.path.join(KEGG_DIR, 'bel')
KEGG_FILES = os.path.join(KEGG_DIR, 'xml')
KEGG_CACHE = os.path.join(KEGG_DIR, 'cache')

#: Reactome
REACTOME = 'reactome'
REACTOME_DIR = os.path.join(DATA_DIR, REACTOME)
REACTOME_BEL = os.path.join(REACTOME_DIR, 'bel')
REACTOME_FILES = os.path.join(REACTOME_DIR, 'rdf')

#: WikiPathways
WIKIPATHWAYS = 'wikipathways'
WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)
WIKIPATHWAYS_BEL = os.path.join(WIKIPATHWAYS_DIR, 'bel')
WIKIPATHWAYS_FILES = os.path.join(WIKIPATHWAYS_DIR, 'rdf', 'wp')

SPIA_DIR = os.path.join(DATA_DIR, 'spia')
CX_DIR = os.path.join(DATA_DIR, 'cx')
PPI_DIR = os.path.join(DATA_DIR, 'ppi')
UNIVERSE_DIR = os.path.join(DATA_DIR, 'universe')


def ensure_pathme_folders():
    """Ensure data folders are created."""
    os.makedirs(KEGG_DIR, exist_ok=True)
    os.makedirs(REACTOME_DIR, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_DIR, exist_ok=True)
    os.makedirs(KEGG_CACHE, exist_ok=True)
    os.makedirs(SPIA_DIR, exist_ok=True)
    os.makedirs(CX_DIR, exist_ok=True)
    os.makedirs(PPI_DIR, exist_ok=True)
    os.makedirs(UNIVERSE_DIR, exist_ok=True)

    os.makedirs(KEGG_BEL, exist_ok=True)
    os.makedirs(REACTOME_BEL, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_BEL, exist_ok=True)

    os.makedirs(KEGG_FILES, exist_ok=True)
    os.makedirs(REACTOME_FILES, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_FILES, exist_ok=True)


ensure_pathme_folders()

UNKNOWN = 'unknown'

# Other namespaces

#: InterPro
INTERPRO = 'interpro'
#: Pfam
PFAM = 'pfam'
#: BRENDA
BRENDA = 'brenda'
#: ChEMBL
CHEMBL = 'chembl'
#: miRBase
MIRBASE = 'mirbase'

#: HGNC
HGNC = 'HGNC'
#: HGNC Symbol
HGNC_SYMBOL = 'HGNC symbol'
#: ENSEMBL
ENSEMBL = 'ENSEMBL'
#: EXPASY
EXPASY = 'EXPASY'
#: ENTREZ
ENTREZ = 'ENTREZ'
#: UniProt
UNIPROT = 'UniProt'
#: ChEBI
CHEBI = 'ChEBI'
#: ChEBI chemical name
CHEBI_NAME = 'ChEBI name'
#: PubChem
PUBCHEM = 'PubChem'
#: Wikipedia
WIKIPEDIA = 'WIKIPEDIA'

#: KEGG entity identifier
KEGG_ID = 'kegg_id'

#: KEGG entity name
KEGG_NAME = 'kegg_name'

#: KEGG entity type
KEGG_TYPE = 'kegg_type'

#: Protein modification types present in KEGG
KEGG_MODIFICATIONS = {
    'phosphorylation': 'Ph',
    'glycosylation': 'Glyco',
    'ubiquitination': 'Ub',
    'methylation': 'Me',
}

ACTIVITY_ALLOWED_MODIFIERS = (
    abundance,
    protein,
    complex_abundance,
    rna,
)

#: Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.
KEGG_CITATION = '10592173'
#: Fabregat, A., et al. (2018). The reactome pathway knowledgebase. Nucleic acids research, 46(D1), D649-D655.
REACTOME_CITATION = '29145629'

#: REST API to KEGG
KEGG_KGML_URL = 'http://rest.kegg.jp/get/{}/kgml'

#: Reactome RDF
RDF_REACTOME = 'ftp://ftp.ebi.ac.uk/pub/databases/RDF/reactome/r67/reactome-biopax.tar.bz2'

#: WikiPathways RDF
RDF_WIKIPATHWAYS = 'http://data.wikipathways.org/20200310/rdf/wikipathways-20200310-rdf-wp.zip'

#: Mapping to compare conversion of entities from KEGG XML (i.e. KGML) to BEL
KEGG_STATS_COLUMN_NAMES = {
    'nodes': 'BEL Nodes',
    'entities': 'XML Entities',
    'edges': 'BEL Edges',
    'interactions': 'XML Interactions',
    'Protein': 'BEL Proteins',
    'gene': 'XML Genes',
    'Composite': 'BEL Composites',
    'RNA': 'BEL RNA Entities',
    'Complex': 'BEL Complexes',
    'group': 'XML Complexes',
    'Abundance': 'BEL Compounds',
    'compound entity': 'XML Compounds',
    'BiologicalProcess': 'BEL Biological Processes',
    'map': 'XML Biological Processes',
    'Reaction': 'BEL Reactions',
    'ortholog': 'XML Orthologs',
    'increases': 'BEL Increase Relations',
    'activation': 'XML Activation Relations',
    'expression': 'XML Expression Relations',
    'reversible': 'XML Reversible Reactions',
    'irreversible': 'XML Irreversible Reactions',
    'phosphorylation': 'XML Phosphorylation Relations',
    'glycosylation': 'XML Glycosylation Relations',
    'ubiquitination': 'XML Ubiquitination Relations',
    'methylation': 'XML Methylation Relations',
    'decreases': 'BEL Decrease Relations',
    'inhibition': 'XML Inhibition Relations',
    'repression': 'XML Repression Relations',
    'dephosphorylation': 'XML Dephosphorylation Relations',
    'association': 'BEL Association Relations',
    'binding/association': 'XML Binding/Association Relations',
    'indirect effect': 'XML Association Relations',
    'hasComponent': 'BEL Component Edges',
    'hasVariant': 'BEL Variant Edges',
    'hasReactant': 'BEL Reactants Edges',
    'hasProduct': 'BEL Products Edges',
    'compound': 'XML Compound Relations',
    'dissociation': 'XML Dissociation Relations',
    'hidden compound': 'XML Hidden Compound Relations',
    'missing interaction': 'XML Missing Interaction Relations',
    'state change': 'XML State Change Relations',
    'brite': 'XML Brite Hierarchy',
}

#: Mapping of BEL terms to their definitions
BEL_STATS_COLUMN_NAMES = {
    'nodes': 'BEL Nodes',
    'edges': 'BEL Edges',
    'Protein': 'BEL Proteins',
    'Gene': 'BEL Genes',
    'Composite': 'BEL Composites',
    'RNA': 'BEL RNA Entities',
    'Complex': 'BEL Complexes',
    'Abundance': 'BEL Compounds',
    'BiologicalProcess': 'BEL Biological Processes',
    'Reaction': 'BEL Reactions',
    'increases': 'BEL Increase Relations',
    'decreases': 'BEL Decrease Relations',
    'regulates': 'BEL Regulate Relations',
    'association': 'BEL Association Relations',
    'hasComponent': 'BEL Component Edges',
    'hasVariant': 'BEL Variant Edges',
    'hasReactant': 'BEL Reactants Edges',
    'hasProduct': 'BEL Products Edges',
    'translatedTo': 'BEL Translation Edges',
}
