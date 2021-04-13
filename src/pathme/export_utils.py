# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import Iterable, List, Optional, TextIO, Tuple, Union
from urllib.parse import urljoin

import click
import networkx as nx
import pandas as pd
import requests
from diffupath.utils import get_dir_list, get_or_create_dir
from tqdm import tqdm

from bio2bel import ensure_path
from bio2bel_kegg.constants import KEGG_ORGANISM_URL, MODULE_NAME
from bio2bel_reactome import Manager as ReactomeManager
from bio2bel_reactome.models import Pathway

from networkx.utils import open_file

import pybel
from pybel import BELGraph, from_pickle, to_pickle, union
from pybel.constants import ANNOTATIONS, NAME, RELATION
from pybel.struct import add_annotation_value, count_functions, remove_isolated_list_abundances
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes
from pybel_tools.analysis.spia import bel_to_spia_matrices, spia_matrices_to_excel

from .constants import KEGG, KEGG_BEL, KEGG_FILES, KEGG_KGML_URL, KEGG_PATHWAYS_URL, \
    PATHME_DIR, REACTOME, REACTOME_BEL, REACTOME_FILES, UNIVERSE_DIR, WIKIPATHWAYS, \
    WIKIPATHWAYS_BEL, WIKIPATHWAYS_FILES
from .normalize_names import normalize_graph_names
from .pybel_utils import flatten_complex_nodes

logger = logging.getLogger(__name__)


def add_annotation_key(graph: BELGraph):
    """Add annotation key in data (in place operation)."""
    for u, v, k in graph.edges(keys=True):
        if ANNOTATIONS not in graph[u][v][k]:
            graph[u][v][k][ANNOTATIONS] = {}


def get_all_pickles(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
) -> Tuple[List[str], List[str], List[str]]:
    """Return a list with all pickle paths."""
    kegg_pickles = get_paths_in_folder(kegg_path or KEGG_BEL)
    if not kegg_pickles:
        logger.warning('No KEGG files found. Please create the BEL KEGG files')

    reactome_pickles = get_paths_in_folder(reactome_path or REACTOME_BEL)
    if not reactome_pickles:
        logger.warning('No Reactome files found. Please create the BEL Reactome files')

    wp_pickles = get_paths_in_folder(wikipathways_path or WIKIPATHWAYS_BEL)
    if not wp_pickles:
        logger.warning('No WikiPathways files found. Please create the BEL WikiPathways files')

    return kegg_pickles, reactome_pickles, wp_pickles


def get_universe_graph(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
    flatten: bool = True,
    normalize_names: bool = True,
) -> BELGraph:
    """Return universe graph."""
    universe_graphs = iterate_universe_graphs(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names,
    )
    # Just keep the graph and not the source
    universe_graphs = (graph for _, _, graph in universe_graphs)
    logger.info('Merging all into a hairball...')
    return union(universe_graphs)


@open_file(1, mode='w')
def export_ppi_tsv(graph: BELGraph, path: Union[str, TextIO]):
    """Export PPI like tsv-file."""
    for u, v, edge_data in graph.edges(data=True):
        # Only export if both node names are present
        if NAME not in u or NAME not in v:
            continue
        print(  # noqa: T001
            u[NAME], edge_data[RELATION], v[NAME],
            sep='\t',
            file=path,
        )


def export_helper(
    *,
    output: str,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
    fmt: str = 'spia',
) -> None:
    """Export helper of PathMe.

    It exporter to SPIA-like format or csv file with names.

    :param output: output directory
    :param kegg_path: directory to KEGG pickles
    :param reactome_path: directory to Reactome pickles
    :param wikipathways_path: directory to WikiPathways pickles
    :param fmt: Export format
    """
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
    )

    paths = kegg_pickles + reactome_pickles + wp_pickles

    logger.info(f'A total of {len(paths)} will be exported')

    paths = tqdm(paths, desc='Exporting PPI files')

    # Call Reactome manager and check that is populated
    reactome_manager = ReactomeManager()
    if not reactome_manager.is_populated():
        logger.warning('Reactome Manager is not populated')

    # Load each pickle and export it as excel file
    for path in paths:
        if not path.endswith('.pickle'):
            continue

        if path in kegg_pickles:
            pathway_graph = from_pickle(os.path.join(kegg_path, path))
            normalize_graph_names(pathway_graph, KEGG)

        elif path in reactome_pickles:
            # Load BELGraph
            pathway_graph = from_pickle(os.path.join(reactome_path, path))

            # Check if pathway has children to build the merge graph
            pathway_id = path[:-len('.pickle')]

            # Look up in Bio2BEL Reactome
            pathway = reactome_manager.get_pathway_by_id(pathway_id)

            # Log if it is not present
            if not pathway:
                logger.warning(f'{pathway_id} not found in database')
                continue

            # Check if there are children and merge them on the fly
            for child in yield_all_children(pathway):

                child_file_path = os.path.join(reactome_path, f"{child.resource_id}.pickle")
                if not os.path.exists(child_file_path):
                    logger.warning(f'{child.resource_id} pickle does not exist')
                    continue

                # Load the pickle and union it
                child_graph = pybel.from_pickle(child_file_path)
                pathway_graph += child_graph

            # Normalize graph names
            normalize_graph_names(pathway_graph, REACTOME)

        elif path in wp_pickles:
            pathway_graph = from_pickle(os.path.join(wikipathways_path, path))
            normalize_graph_names(pathway_graph, WIKIPATHWAYS)

        else:
            logger.warning(f'Unknown pickle file: {path}')
            continue

        # Explode complex nodes
        flatten_complex_nodes(pathway_graph)

        # Collapse nodes
        collapse_all_variants(pathway_graph)
        collapse_to_genes(pathway_graph)

        if fmt == 'spia':
            # Default SPIA exporter
            spia_matrices = bel_to_spia_matrices(pathway_graph)

            _name = path[:-len('.pickle')]
            output_file = os.path.join(output, f"{_name}.xlsx")

            if os.path.isfile(output_file):
                continue

            # Export excel file representing the connectivity matrix of the BEL Graph
            spia_matrices_to_excel(spia_matrices, output_file)

        elif fmt == 'ppi':
            _name = path[:-len('.pickle')]
            output_file = os.path.join(output, f"{_name}.tsv")
            export_ppi_tsv(pathway_graph, output_file)
        else:
            raise ValueError(f'Unknown export format: {fmt}')


def iterate_indra_statements(**kwargs) -> Iterable['import indra as indra]:
    """Iterate over INDRA statements for the universe."""
    for _, _, graph in iterate_universe_graphs(**kwargs):
        yield from pybel.to_indra_statements(graph)


def iterate_universe_graphs(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
    flatten: bool = True,
    normalize_names: bool = True,
) -> Iterable[Tuple[str, str, BELGraph]]:
    """Return universe graph."""
    kegg_pickle_paths, reactome_pickle_paths, wp_pickle_paths = get_all_pickles(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
    )

    n_paths = len(kegg_pickle_paths) + len(reactome_pickle_paths) + len(wp_pickle_paths)
    logger.info(f'{n_paths} graphs will be put in the universe')

    yield from _iterate_wp(wp_pickle_paths, wikipathways_path, flatten, normalize_names)
    yield from _iterate_kegg(kegg_pickle_paths, kegg_path, flatten, normalize_names)
    yield from _iterate_reactome(reactome_pickle_paths, reactome_path, flatten, normalize_names)


def _iterate_wp(wp_pickle_paths, wikipathways_path, flatten, normalize_names):
    for path in tqdm(wp_pickle_paths, desc=f'Loading WP pickles from {wikipathways_path}'):
        if not path.endswith('.pickle'):
            continue

        graph = from_pickle(os.path.join(wikipathways_path, path), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, WIKIPATHWAYS)

        _update_graph(graph, path, WIKIPATHWAYS)
        yield WIKIPATHWAYS, path, graph


def _iterate_kegg(kegg_pickle_paths, kegg_path, flatten, normalize_names):
    for path in tqdm(kegg_pickle_paths, desc=f'Loading KEGG pickles from {kegg_path}'):
        if not path.endswith('.pickle'):
            continue
        graph = from_pickle(os.path.join(kegg_path, path), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, KEGG)

        _update_graph(graph, path, KEGG)
        yield KEGG, path, graph


def _iterate_reactome(reactome_pickle_paths, reactome_path, flatten, normalize_names):
    for file in tqdm(reactome_pickle_paths, desc=f'Loading Reactome pickles from {reactome_path}'):
        if not file.endswith('.pickle'):
            continue

        graph = from_pickle(os.path.join(reactome_path, file), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, REACTOME)

        _update_graph(graph, file, REACTOME)
        yield REACTOME, file, graph


def _update_graph(graph, file, database):
    graph.annotation_list['database'] = {KEGG, REACTOME, WIKIPATHWAYS}
    add_annotation_key(graph)
    add_annotation_value(graph, 'database', database)
    graph.annotation_pattern['PathwayID'] = '.*'
    _name = file[:-len('.pickle')]
    add_annotation_value(graph, 'PathwayID', _name)


def _munge_node_attribute(node, attribute='name'):
    """Munge node attribute."""
    if node.get(attribute) is None:
        return str(node)
    else:
        return node.get(attribute)


def to_gml(graph: pybel.BELGraph, path: str = PATHME_DIR) -> None:
    """Write this graph to GML  file using :func:`networkx.write_gml`."""
    rv = nx.MultiDiGraph()

    for node in graph:
        rv.add_node(
            _munge_node_attribute(node, 'name'),
            namespace=str(node.get('namespace')),
            function=node.get('function'),
        )

    for u, v, key, edge_data in graph.edges(data=True, keys=True):
        rv.add_edge(
            _munge_node_attribute(u),
            _munge_node_attribute(v),
            interaction=str(edge_data[RELATION]),
            bel=str(edge_data),
            key=str(key),
        )

    nx.write_gml(rv, path)


def get_paths_in_folder(directory: str) -> List[str]:
    """Return the files in a given folder.

    :param directory: folder path
    :return: file names in folder
    """
    return [
        path
        for path in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, path))
    ]


def yield_all_children(pathway: Pathway) -> Iterable[Pathway]:
    """Transverse recursively the Reactome hierarchy and return all children for a given pathway."""
    if pathway.children:
        for child in pathway.children:
            yield child
            yield from yield_all_children(child)


def download_kgml_files(kegg_pathway_ids, path=KEGG_FILES):
    """Download KEGG KGML files by querying the KEGG API.

    :param list kegg_pathway_ids: list of kegg ids

    :param path: Location of KEGG files
    """
    for kegg_id in tqdm.tqdm(kegg_pathway_ids, desc='Downloading KEGG files'):
        request = requests.get(KEGG_KGML_URL.format(kegg_id))
        with open(os.path.join(path, f'{kegg_id}.xml'), 'w+') as file:
            file.write(request.text)


def get_kegg_pathway_ids(connection=None, populate=False, species='hsa'):
    """Return a list of all pathway identifiers stored in the KEGG database.

    :param Optional[str] connection: connection to the database
    :param populate: Whether to populate
    :param species: Defaults to Homo Sapiens
    :returns: list of all kegg_pathway_ids
    :rtype: list
    """
    return None


def get_organisms_df(url: Optional[str] = None) -> pd.DataFrame:
    """Convert tab separated txt files to pandas Dataframe.

    :param url: url from KEGG tab separated file
    :return: dataframe of the file
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(
        url or ensure_path(MODULE_NAME, KEGG_ORGANISM_URL, path='organisms.tsv'),
        sep='\t',
        header=None,
        names=[
            'kegg_id',
            'kegg_code',
            'name',
            # fourth column is the taxonomy hierarchy
        ],
        usecols=[0, 1, 2],
    )
    df['common_name'] = df['name'].map(lambda name: name.replace(')', '').split(' (')[1].capitalize() if len(
        name.replace(')', '').split(' (')) > 1 else '')
    df['name'] = df['name'].map(lambda name: name.replace(')', '').split(' (')[0].capitalize())
    return df


def get_df_value(df, comparation_column, comparation_value, name_column):
    """Populate pathways for A SINGLE specie."""
    return df.loc[df[comparation_column] == comparation_value][name_column].values[0]


def get_pathways_kegg_id(specie_id):
    """Get patwhway KEGG ID for A SINGLE specie."""
    df_pathway_species = get_organisms_df()

    specie_id = specie_id.replace('_', ' ')

    if specie_id.lower() in list(df_pathway_species['kegg_code']):
        kegg_id = specie_id

    elif specie_id.capitalize() in list(df_pathway_species['name']):
        kegg_id = df_pathway_species.loc[df_pathway_species['name'] == specie_id]['kegg_code'].values[0]

    elif specie_id.capitalize() in list(df_pathway_species['common_name']):
        kegg_id = df_pathway_species.loc[df_pathway_species['common_name'] == specie_id]['kegg_code'].values[0]

    else:
        raise Warning(f'Organism id {specie_id} not found in KEGG.')

    return kegg_id


def get_common_or_name_specie_id(specie_id, common=True):
    """Get common name or name (from its oposite combination) from KEGG mapping for A SINGLE specie.
       If common_name given, name returned, if name given common_name returned."""
    df_pathway_species = get_organisms_df()

    specie_id = specie_id.replace('_', ' ')

    if specie_id.capitalize() in list(df_pathway_species['name']):
        specie_name = df_pathway_species.loc[df_pathway_species['name'] == specie_id]['common_name'].values[0]

    elif specie_id.capitalize() in list(df_pathway_species['common_name']):
        specie_name = df_pathway_species.loc[df_pathway_species['common_name'] == specie_id]['name'].values[0]

    elif specie_id.lower() in list(df_pathway_species['kegg_code']):
        if common:
            specie_name = df_pathway_species.loc[df_pathway_species['kegg_code'] == specie_id]['common_name'].values[0]
        else:
            specie_name = df_pathway_species.loc[df_pathway_species['kegg_code'] == specie_id]['name'].values[0]

    else:
        raise Warning(f'Organism id {specie_id} not found in KEGG.')

    return specie_name


def get_pathway_kegg_url(specie_id):
    """Get patwhway URL for A SINGLE specie."""
    return urljoin(KEGG_PATHWAYS_URL, get_pathways_kegg_id(specie_id))


def get_all_pathways_organism(url=None) -> list:
    """Convert tab separated txt files to pandas Dataframe.

    :param url: url from KEGG tab separated file
    :return: dataframe of the file
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(
        url,
        sep='\t',
        header=None,
        names=[
            'pathway_id',
            # fourth column is the taxonomy hierarchy
        ],
        usecols=[1],
    )

    return [pathway_id.replace('path:', '') for pathway_id in df['pathway_id']]


def generate_universe(kegg_path=KEGG_FILES,
                      reactome_path=REACTOME_FILES,
                      wikipathways_path=WIKIPATHWAYS_FILES,
                      output=UNIVERSE_DIR,
                      no_flatten=False,
                      no_normalize_names=False,
                      specie='Homo_sapiens'):
    """Export harmonized PathMe universe."""
    flatten = not no_flatten
    normalize_names = not no_normalize_names

    # Specie name treatment.
    specie = specie.replace(' ', '_').capitalize()
    specie_altern_name = get_common_or_name_specie_id(specie).replace(' ', '_').capitalize()

    if not flatten:
        click.secho('Complexes and Reactions will be not be flatten to single nodes')

    if not normalize_names:
        click.secho('Names will not be normalized to lower case')

    # KEGG specie processing
    kegg_species_dir_list = get_dir_list(kegg_path, True)
    kegg_path = os.path.join(kegg_path, specie)

    if specie not in kegg_species_dir_list and specie_altern_name not in kegg_species_dir_list:
        kegg_ids = get_all_pathways_organism(get_pathway_kegg_url(specie))
        click.secho(
            'You are about to download KGML files from KEGG.\n'
            'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
            ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
            'We (PathMe developers) are not responsible for the end use of this data.\n',
        )
        os.makedirs(kegg_path)
        download_kgml_files(kegg_ids, path=kegg_path)

    # Reactome specie processing
    specie_file = f'{specie}.owl'
    specie_alt_file = f'{specie_altern_name}.owl'

    reactome_species_file_list = get_or_create_dir(reactome_path)

    if specie_file in reactome_species_file_list:
        reactome_path = os.path.join(reactome_path, specie_file)
    elif specie_alt_file in reactome_species_file_list:
        reactome_path = os.path.join(reactome_path, specie_alt_file)
    else:
        click.secho('Specie not found in the populated Reactome resources.')

    # WikiPathways specie processing
    wikipath_species_dir_list = get_dir_list(wikipathways_path, True)

    if specie in wikipath_species_dir_list:
        wikipathways_path = os.path.join(wikipathways_path, specie)
    elif specie_altern_name in wikipath_species_dir_list:
        wikipathways_path = os.path.join(wikipathways_path, specie_altern_name)
    else:
        click.secho('Specie not found in the populated Wikipathways resources.')

    click.secho("Merging graphs to universe and harmonizing...(this might take a while)")

    # Not explode will flip the boolean coming from the cli
    universe_graph = get_universe_graph(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names,
    )
    click.secho(f'Number of isolates after getting universe: {nx.number_of_isolates(universe_graph)}')

    # Remove isolated list abundances
    remove_isolated_list_abundances(universe_graph)

    if flatten:
        # TODO: Remove node list solo de Reactome
        click.secho(f'Number of isolates after flattening: {nx.number_of_isolates(universe_graph)}')

    click.secho("Merging variants and genes")
    collapse_all_variants(universe_graph)
    collapse_to_genes(universe_graph)
    click.secho(f'Number of isolates after collapsing variants and to genes: {nx.number_of_isolates(universe_graph)}')

    universe_graph.name = 'PathMe Universe'

    file_name = os.path.join(output, '_'.join([specie, 'pathme_universe.pickle']))
    click.secho(f"Export BEL graph to: {file_name}")
    click.secho(universe_graph.summary_str())
    click.secho(count_functions(universe_graph))

    to_pickle(universe_graph, file_name)
