# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling WikiPathways RDF and data."""

import logging
import os
import re
import zipfile
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

import networkx as nx

from bio2bel_hgnc import Manager as HgncManager
from bio2bel_wikipathways import Manager as WikiPathwaysManager
from pybel import BELGraph
from ..constants import (
    BRENDA, CHEMBL, DATA_DIR, ENSEMBL, ENTREZ, EXPASY, HGNC, INTERPRO, KEGG, MIRBASE, PFAM, REACTOME, UNIPROT,
    WIKIPATHWAYS, WIKIPEDIA,
)
from ..export_utils import get_paths_in_folder

WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

logger = logging.getLogger(__name__)


def evaluate_wikipathways_metadata(metadata: Union[str, Set[str]]) -> str:
    """Evaluate metadata in wikipathways and return the string representation."""
    if isinstance(metadata, set):
        return ','.join(metadata)

    return metadata


def _get_update_alias_symbol(
    hgnc_manager: HgncManager,
    original_identifier: str,
    original_namespace: str,
) -> Tuple[str, str, str]:
    """Try to get current alias symbol.

    :param hgnc_manager: hgnc manager
    :param original_identifier:
    :param original_namespace:
    """
    query_result = hgnc_manager.get_hgnc_from_alias_symbol(original_identifier)

    if not query_result:
        logger.debug('No found HGNC Symbol for id %s in (%s)', original_identifier, original_namespace)
        return original_namespace, original_identifier, original_identifier

    return HGNC, query_result.symbol, query_result.identifier


def _validate_query(
    hgnc_manager: HgncManager,
    query_result,
    original_identifier: str,
    original_namespace: str,
) -> Tuple[str, str, str]:
    """Process and validate HGNC query.

    :param hgnc_manager: hgnc manager
    :param query_result:
    :param original_identifier:
    :param original_namespace:
    """
    # If invalid entry from HGNC, try to find updated symbol
    if not query_result and original_namespace == HGNC:
        return _get_update_alias_symbol(hgnc_manager, original_identifier, HGNC)

    # Invalid entry, proceed with invalid identifier
    if not query_result:
        logger.debug('No found HGNC Symbol for id %s in (%s)', original_identifier, original_namespace)
        return original_namespace, original_identifier, original_identifier

    # Multiple entries are returned, for UniProt identifiers
    if isinstance(query_result, list):
        if len(query_result) > 1:
            logger.debug('UniProt identifier with multiple HGNC:s %s', query_result)
        query_result = query_result[0]

    # Correct entry, use HGNC identifier
    return HGNC, query_result.symbol, query_result.identifier


def get_valid_gene_identifier(node_ids_dict, hgnc_manager: HgncManager, pathway_id) -> Tuple[str, str, str]:
    """Return protein/gene identifier for a given RDF node.

    :param dict node_ids_dict: node dictionary
    :param hgnc_manager: hgnc manager
    :param pathway_id: ID of pathway
    :return: namespace, name, identifier
    """
    # Try to get hgnc symbol
    if 'bdb_hgncsymbol' in node_ids_dict or 'hgnc' in node_ids_dict['uri_id']:

        if 'hgnc' in node_ids_dict['uri_id']:
            hgnc_entry = hgnc_manager.get_gene_by_hgnc_id(node_ids_dict['identifier'])
            if not hgnc_entry:
                hgnc_symbol = node_ids_dict['name']
            else:
                hgnc_symbol = hgnc_entry.symbol
        else:
            hgnc_symbol = check_multiple(node_ids_dict['bdb_hgncsymbol'], 'bdb_hgncsymbol', pathway_id)
            hgnc_entry = hgnc_manager.get_gene_by_hgnc_symbol(hgnc_symbol)

        return _validate_query(hgnc_manager, hgnc_entry, hgnc_symbol, HGNC)

    # Try to get ENTREZ id
    elif 'bdb_ncbigene' in node_ids_dict or 'ncbiprotein' in node_ids_dict['uri_id']:
        if 'bdb_ncbigene' in node_ids_dict:
            entrez_id = check_multiple(node_ids_dict['bdb_ncbigene'], 'bdb_ncbigene', pathway_id)
        elif 'ncbiprotein' in node_ids_dict['uri_id']:
            entrez_id = check_multiple(node_ids_dict['identifier'], 'ncbiprotein', pathway_id)
        else:
            raise ValueError(f'Missing entrez gene identifier [pathway={pathway_id}]')

        hgnc_entry = hgnc_manager.get_gene_by_entrez_id(entrez_id)

        return _validate_query(hgnc_manager, hgnc_entry, entrez_id, ENTREZ)

    # Try to get UniProt id
    elif 'bdb_uniprot' in node_ids_dict:
        uniprot_id = check_multiple(node_ids_dict['bdb_uniprot'], 'bdb_uniprot', pathway_id)
        hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(uniprot_id)

        return _validate_query(hgnc_manager, hgnc_entry, uniprot_id, UNIPROT)

    # Try to get ENSEMBL id
    elif 'bdb_ensembl' in node_ids_dict or 'ena.embl' in node_ids_dict['uri_id']:
        if 'bdb_ensembl' in node_ids_dict:
            ensembl_id = check_multiple(node_ids_dict['bdb_ensembl'], 'bdb_ensembl', pathway_id)

        elif 'ena.embl' in node_ids_dict['uri_id']:
            ensembl_id = check_multiple(node_ids_dict['identifier'], 'bdb_ensembl', pathway_id)

        else:
            raise ValueError(f'Missing ensemble identifier [pathway={pathway_id}]')

        hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(ensembl_id)

        return _validate_query(hgnc_manager, hgnc_entry, ensembl_id, ENSEMBL)

    elif 'ec-code' in node_ids_dict['uri_id']:
        ec_number = check_multiple(node_ids_dict['name'], 'ec-code', pathway_id)
        return EXPASY, ec_number, ec_number

    # Only wikipathways identifier is given
    elif 'bdb_wikidata' in node_ids_dict:
        name = check_multiple(node_ids_dict['name'], 'wikidata', pathway_id)

        # Find out whether the name is a valid HGNC symbol
        hgnc_entry = hgnc_manager.get_gene_by_hgnc_symbol(name)

        # Correct entry, use HGNC identifier
        if hgnc_entry:
            return HGNC, hgnc_entry.symbol, hgnc_entry.identifier

        logger.debug('Adding WikiPathways node %s (%s)', name, WIKIPATHWAYS)
        return WIKIPATHWAYS, name, name

    elif WIKIPEDIA.lower() in node_ids_dict['uri_id']:
        wiki_name = check_multiple(node_ids_dict['identifier'], 'wikipedia_id', pathway_id)
        wiki_id = check_multiple(node_ids_dict['name'], 'wikipedia_name', pathway_id)
        logger.debug('Adding Wikipedia node %s (%s)', wiki_name, WIKIPATHWAYS)
        return WIKIPEDIA, wiki_name, wiki_id

    elif KEGG.lower() in node_ids_dict['uri_id']:
        kegg_id = check_multiple(node_ids_dict['identifier'], 'kegg_id', pathway_id)
        kegg_name = check_multiple(node_ids_dict['name'], 'kegg_name', pathway_id)
        logger.debug('Adding KEGG node %s ', kegg_id)
        return KEGG, kegg_name, kegg_id

    elif INTERPRO.lower() in node_ids_dict['uri_id']:
        interpro_id = check_multiple(node_ids_dict['identifier'], 'interpro_id', pathway_id)
        interpro_name = check_multiple(node_ids_dict['name'], 'interpro_name', pathway_id)
        logger.debug('Adding InterPro node %s ', interpro_id)
        return INTERPRO, interpro_name, interpro_id

    elif PFAM.lower() in node_ids_dict['uri_id']:
        pfam_id = check_multiple(node_ids_dict['identifier'], 'pfam_id', pathway_id)
        pfam_name = check_multiple(node_ids_dict['name'], 'pfam_name', pathway_id)
        logger.debug('Adding Pfam node %s ', pfam_id)
        return PFAM, pfam_name, pfam_id

    elif 'mirbase.mature' in node_ids_dict['uri_id']:
        mirbase_id = check_multiple(node_ids_dict['identifier'], 'mirbase_id', pathway_id)
        mirbase_name = check_multiple(node_ids_dict['name'], 'mirbase_name', pathway_id)
        logger.debug('Adding miRBase node %s ', mirbase_id)
        return MIRBASE, mirbase_name, mirbase_id

    elif 'chembl.compound' in node_ids_dict['uri_id']:
        chembl_id = check_multiple(node_ids_dict['identifier'], 'chembl_id', pathway_id)
        chembl_name = check_multiple(node_ids_dict['name'], 'chembl_name', pathway_id)
        logger.debug('Adding ChEMBL node %s ', chembl_id)
        return CHEMBL, chembl_name, chembl_id

    elif 'brenda' in node_ids_dict['uri_id']:
        brenda_id = check_multiple(node_ids_dict['identifier'], 'brenda', pathway_id)
        brenda_name = check_multiple(node_ids_dict['name'], 'brenda', pathway_id)
        logger.debug('Adding BRENDA node %s ', brenda_id)
        return BRENDA, brenda_name, brenda_id

    elif 'insdc' in node_ids_dict['uri_id']:
        indsc_id = check_multiple(node_ids_dict['identifier'], 'insdc', pathway_id)
        indsc_name = check_multiple(node_ids_dict['name'], 'insdc', pathway_id)
        return HGNC, indsc_name, indsc_id

    # Nodes from reactome pointing to a gene
    elif 'reactome' in node_ids_dict['uri_id']:
        return REACTOME, node_ids_dict['name'], node_ids_dict['identifier']

    raise Exception('Unknown identifier for node %s', node_ids_dict)


MULTIPLE_RE = re.compile('^[A-Z0-9]+$')


def check_multiple(element, element_name, pathway_id):
    """Check whether element is iterable.

    :param element: variable to check
    :param element_name: name to print
    :param pathway_id: ID of pathway
    :return:
    """
    if isinstance(element, (set, list)):
        _pid = pathway_id.split('/')[-1]
        logger.debug(f'Multiple values for "{element_name}": {element} [{_pid}]')
        # TODO: print the WikiPathways bps that return a set because they are probably wrong.
        if len(element) == 1:
            return list(element)[0]

        if len(element) > 1:
            for subelement in element:
                if MULTIPLE_RE.match(subelement):
                    return subelement

            return list(element)[0]

        logger.debug('Empty list/set %s', element)

    return element


def merge_two_dicts(dict1, dict2):
    """Merge two dictionaries.

    :param dict dict1:
    :param dict dict2:
    :returns: merged_dict
    :rtype: dict
    """
    merged_dict = dict1.copy()  # start with x's keys and values
    merged_dict.update(dict2)  # modifies z with y's keys and values & returns None
    return merged_dict


def convert_to_nx(
    nodes: Dict[str, Dict],
    interactions: List[Tuple[str, str, Dict]],
    pathway_info: Dict,
) -> nx.MultiDiGraph:
    """Generate a NetworkX Graph from a network data structure (dict with nodes and edges).

    :param nodes: Node id as keys and Node attributes as values
    :param interactions: list of interactions
    :param pathway_info: pathway info dictionary
    """
    graph = nx.MultiDiGraph(graph_att=pathway_info)

    for node, attributes in nodes.items():
        graph.add_node(node, attr_dict=attributes)

    for subj, obj, interaction in interactions:
        graph.add_edge(subj, obj, attr=interaction)

    return graph


def debug_pathway_info(bel_graph: BELGraph, pathway_path: str, **kwargs):
    """Debug information about the pathway graph representation.

    :param bel_graph: bel graph
    :param pathway_path: path of the pathway
    """
    logger.debug('Pathway id: %s', os.path.basename(pathway_path))

    pathway_name = bel_graph.name
    logger.debug('Pathway Name: %s', pathway_name)

    bel_nodes = bel_graph.number_of_nodes()
    bel_edges = bel_graph.number_of_edges()

    logger.debug('Nodes imported to BEL: %s', bel_nodes)
    logger.debug('Edges imported to BEL: %s', bel_edges)

    if 'statistics' in kwargs:
        statistics = kwargs.get('statistics')
        logger.debug('RDF Nodes statistics: %s', statistics['RDF nodes'])
        logger.debug('RDF Edges statistics: %s', statistics['RDF interactions'])


def debug_global_statistics(global_statistics):
    """Debug pathway statistics.

    :param dict global_statistics: pathway statistics
    """
    for statistics_type, rdf_types in global_statistics.items():
        logger.debug('Total statistics for %s', statistics_type)

        for rdf_type, value in rdf_types.items():
            logger.debug('%s: %s', rdf_type, value)


"""Download utilities"""


def get_file_name_from_url(url: str) -> str:
    """Get the last part of an URL."""
    return url.rsplit('/', 1)[1]


def unzip_file(file_path: str, export_folder: str):
    """Unzip file into a destination folder.

    :param file_path: name of the file
    :param export_folder: name of the file
    """
    zip_ref = zipfile.ZipFile(file_path, 'r')
    zip_ref.extractall(export_folder)
    zip_ref.close()


def filter_wikipathways_files(file_names: Iterable[str]) -> List[str]:
    """Filter files that have not 'ttl' extension or not start with 'WP'."""
    return [
        file_name
        for file_name in file_names
        if file_name.startswith('WP') and file_name.endswith('.ttl')
    ]


def iterate_wikipathways_paths(
    directory: str,
    connection: Optional[str] = None,
    only_canonical: bool = True,
) -> List[str]:
    """Get WikiPathways RDF files in folder.

    :param directory: folder path
    :param connection: database connection
    :param only_canonical: only identifiers present in WP bio2bel db
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(
            f'{directory} does not exist. Please ensure you have downloaded WikiPathways using '
            f'the "pathme wikipathways download" command or you have passed the right argument.',
        )

    paths = get_paths_in_folder(directory)

    # Filter files in folder that have no turtle extension or do not start with 'WP'
    paths = filter_wikipathways_files(paths)

    # Skip files not present in wikipathways bio2bel db -> stuffs from reactome and so on...
    if only_canonical:
        wikipathways_manager = WikiPathwaysManager(connection)
        if not wikipathways_manager.is_populated():
            wikipathways_manager.populate()

        wikipathways_identifiers = {
            pathway.resource_id
            for pathway in wikipathways_manager.get_all_pathways()
        }

        paths = [
            path
            for path in paths
            if path.split('.')[0] in wikipathways_identifiers
        ]

    return paths
