# -*- coding: utf-8 -*-

"""This module contains the methods that run SPARQL queries to convert the Reactome Pathways to BEL."""

import logging
import os
from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple, Union

import rdflib
from rdflib import URIRef
from rdflib.namespace import DC, DCTERMS, Namespace, OWL, RDF, RDFS, SKOS, XSD
from tqdm import tqdm

from pybel import BELGraph, to_pickle
from .convert_to_bel import convert_to_bel
from ..constants import REACTOME_BEL
from ..utils import get_pathway_statitics, parse_rdf, query_result_to_dict

logger = logging.getLogger(__name__)

"""SPARQL string queries"""

#: SPARQL prefixes.
PREFIXES = {
    'owl': OWL,
    'xsd': XSD,
    'rdfs': RDFS,
    'rdf': RDF,
    'dcterms': DCTERMS,
    'dc': DC,
    'skos': SKOS,
    'foaf': Namespace('http://xmlns.com/foaf/0.1/'),
    'dbpedia': Namespace('http://dbpedia.org/property/'),
    'dbpedia2': Namespace('http://dbpedia.org/'),
    'biopax3': Namespace('http://www.biopax.org/release/biopax-level3.owl#'),
}

#: SPARQL query string to get all the  primary types of entries (Pathway, BiochemicalReaction) in a pathway network.
GET_ALL_TYPES = """
SELECT DISTINCT (STRAFTER(STR(?rdf_type), str(biopax3:)) AS ?entry_type)
WHERE {
    ?uri_id rdf:type ?rdf_type .
}
"""

#: SPARQL query string to get pathway URIs and names in the RDF file.
GET_ALL_PATHWAYS = """
SELECT DISTINCT ?uri_id ?name
WHERE {
    ?uri_id rdf:type biopax3:Pathway .
    ?uri_id biopax3:displayName ?name .
}
"""

#: SPARQL query string to get all components of a pathway (predicate biopax3:pathwayComponent).
GET_ALL_PATHWAY_COMPONENTS = """
SELECT DISTINCT ?uri_id ?name ?comment (STRAFTER(STR(?uri_type), str(biopax3:)) AS ?component_type)
WHERE {
    ?pathway biopax3:pathwayComponent ?uri_id .
    ?uri_id rdf:type ?uri_type .
    optional {?uri_id biopax3:displayName ?name .}
    optional {?uri_id biopax3:comment ?comment .}
}
"""

#: SPARQL query string to get all participants in an interaction and its controlType (ACTIVATION or INHIBITION).
GET_INTERACTION_PARTICIPANTS_AND_TYPE = """
SELECT DISTINCT
    (STRAFTER(STR(?component), '#') AS ?identifier)
    ?reactant
    ?product
    (STR(?control_type) AS ?interaction_type)
WHERE {
    ?component biopax3:left ?reactant .
    ?component biopax3:right ?product .
    optional {?control biopax3:controlled ?component .}
    optional {?control biopax3:controlType ?control_type }
}
"""

#: SPARQL query to get all the possible metadate (optional statements) of an entity (Protein, Dna, Pathway...).
GET_ENTITY_METADATA = """
SELECT DISTINCT
    (STRAFTER (STR(?uri_type), str(biopax3:)) AS ?entity_type)
    (STRAFTER(STR(?entity), '#') AS ?identifier)
    (STR(?entity) AS ?uri_id)
    (STRAFTER(STR(?entity), '#') AS ?reactome_id)
    (STR(?entity) AS ?uri_reactome_id)
    (STR(?entity_reference) AS ?uri_id)
    ?name
    ?cell_locat
    ?display_name
    ?complex_components
    ?comment
WHERE {
    ?entity rdf:type ?uri_type .
    optional {?entity biopax3:comment ?comment .}
    optional {?entity biopax3:entityReference ?entity_reference .}
    optional {?entity biopax3:name ?name .}
    optional {?entity biopax3:displayName ?display_name .}
    optional {?entity biopax3:cellularLocation ?cell_locat .}
    optional {?entity biopax3:organism ?organism .}
    optional {?entity biopax3:component ?complex_components .}
}
"""

"""Queries managers"""


def _get_all_entry_types(rdf_graph: rdflib.Graph) -> Set[str]:
    """Get all entries primary types.

    :param rdf_graph:
    :return: set with all entry primary types.
    """
    types_query = rdf_graph.query(
        GET_ALL_TYPES,
        initNs=PREFIXES,
    )

    return {
        str(entry.entry_type)
        for entry in types_query
    }


def _get_pathway_metadata(pathway_uri: str, rdf_graph: rdflib.Graph) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get metadata for a pathway entry.

    :param pathway_uri: URI reference of the queried graph
    :param rdf_graph: RDF Reactome Universe graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    return query_result_to_dict(
        rdf_graph.query(
            GET_ENTITY_METADATA,
            initNs=PREFIXES,
            initBindings={'entity': pathway_uri},
        ),
        attr_empty=['display_name', 'identifier', 'uri_id', 'uri_reactome_id', 'comment'],
        id_dict=False,
    )


def _get_entity_metadata(entity: rdflib.URIRef, rdf_graph: rdflib.Graph) -> Dict[str, Union[str, Set[str]]]:
    """Get the metadata for an entity (Protein, Dna, Complex...).

    :param entity: URI reference of the queried entity
    :param rdf_graph: RDF Reactome Universe graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    entity_metadata = query_result_to_dict(
        rdf_graph.query(
            GET_ENTITY_METADATA,
            initNs=PREFIXES,
            initBindings={'entity': entity},
        ),
        attr_empty=['entity_type'],
        id_dict=False,
    )
    # Complexes might contain multiple components entities so we iterate over
    # the complex components to fetch that information
    if entity_metadata['entity_type'] == 'Complex':

        complex_components = entity_metadata.get('complex_components')
        entity_metadata['complex_components'] = []

        if isinstance(complex_components, str):
            complex_component = _get_entity_metadata(URIRef(complex_components), rdf_graph)
            entity_metadata['complex_components'].append(complex_component)

        elif complex_components:
            for complex_component in complex_components:
                complex_component = _get_entity_metadata(URIRef(complex_component), rdf_graph)
                entity_metadata['complex_components'].append(complex_component)

    return entity_metadata


def _get_reaction_participants(component_uri: str, component, rdf_graph: rdflib.Graph) -> Tuple[
    Dict[Union[str, Set[str]], Dict[str, Union[str, Set[str]]]], Dict[Any, Dict[str, Any]],
]:
    """Get reaction participants (nodes and interactions) for a given reaction.

    :param component_uri: URI reference of the queried reaction component
    :param component: Reaction component metadata
    :param rdf_graph: RDF Reactome Universe graph object
    :return: returns the reaction participants as entities (Proteins, Complex, SmallMolecule...)
     and proteins (the reaction link)
    """
    interactions = {}
    nodes = {}

    spaqrl_reaction_participants = rdf_graph.query(
        GET_INTERACTION_PARTICIPANTS_AND_TYPE,
        initNs=PREFIXES,
        initBindings={'component': component_uri},
    )

    for interaction in spaqrl_reaction_participants:

        reactant_metadata = _get_entity_metadata(interaction.reactant, rdf_graph)
        product_metadata = _get_entity_metadata(interaction.product, rdf_graph)

        reactant_id = reactant_metadata['identifier']
        product_id = product_metadata['identifier']

        nodes[reactant_id] = reactant_metadata
        nodes[product_id] = product_metadata

        if interaction.identifier not in interactions:
            interactions[interaction.identifier] = {'metadata': component}

        if 'participants' not in interactions[interaction.identifier]:
            interactions[interaction.identifier]['participants'] = (reactant_id, product_id)

        else:
            interaction_participants = interactions[interaction.identifier]['participants']
            if isinstance(interaction_participants, tuple):
                interactions[interaction.identifier]['participants'] = {
                    'reactants': {interaction_participants[0], reactant_id},
                    'products': {interaction_participants[1], product_id},
                }
            else:
                interactions[interaction.identifier]['participants']['reactants'].add(reactant_id)
                interactions[interaction.identifier]['participants']['products'].add(product_id)

        interactions[interaction.identifier]['metadata']['interaction_type'] = str(interaction.interaction_type)

    return nodes, interactions


def _get_pathway_components(
    pathway_uri: rdflib.URIRef,
    rdf_graph: rdflib.Graph,
) -> Tuple[Dict[str, Dict[str, Union[str, Set[str]]]], List[Dict[str, Union[str, Set[str]]]]]:
    """Get components (nodes and interactions) for a given pathway.

    :param pathway_uri: URI reference of the queried pathway
    :param rdf_graph: RDF Reactome Universe graph object
    :return: returns the pathway components as entities (Proteins, Complex, SmallMolecule...) and proteins (their links)
    """
    interactions = {}
    nodes = {}

    spaqrl_pathway_components = rdf_graph.query(
        GET_ALL_PATHWAY_COMPONENTS,
        initNs=PREFIXES,
        initBindings={'pathway': pathway_uri},
    )

    pathway_components = query_result_to_dict(spaqrl_pathway_components)

    for component_uri, component in pathway_components.items():

        if component['component_type'] == 'BiochemicalReaction':
            component_nodes, component_interactions = _get_reaction_participants(component_uri, component, rdf_graph)

            nodes.update(component_nodes)
            interactions.update(component_interactions)

        elif component['component_type'] == 'Pathway':
            pathway_metadata = _get_pathway_metadata(component_uri, rdf_graph)

            nodes[pathway_metadata['uri_reactome_id']] = pathway_metadata

    return nodes, list(interactions.values())


def get_reactome_statistics(resource_file, hgnc_manager, chebi_manager):
    """Get types statistics for Reactome.

    :param str resource_file: RDF file
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    """
    logger.info('Parsing Reactome RDF file')
    rdf_graph = parse_rdf(resource_file, fmt='xml')

    spaqrl_all_pathways = rdf_graph.query(GET_ALL_PATHWAYS, initNs=PREFIXES)

    global_statistics = defaultdict(lambda: defaultdict(int))

    for pathway_uri, _pathway_title in tqdm(spaqrl_all_pathways, desc='Generating Reactome Statistics'):
        nodes, edges = _get_pathway_components(pathway_uri, rdf_graph)
        pathway_metadata = _get_pathway_metadata(pathway_uri, rdf_graph)

        nodes_types = [
            node['entity_type'] for node in nodes.values()
        ]
        edges_types = [
            edge['metadata']['interaction_type'] for edge in edges
        ]

        bel_graph = convert_to_bel(nodes, edges, pathway_metadata, hgnc_manager, chebi_manager)

        global_statistics, pathway_statistics = get_pathway_statitics(
            nodes_types, edges_types, bel_graph, global_statistics=global_statistics,
        )

    return global_statistics


def reactome_pathway_to_bel(pathway_uri, rdf_graph, hgnc_manager, chebi_manager) -> BELGraph:
    """Convert a Reactome pathway to BEL.

    :param pathway_uri: URI reference of the queried graph
    :param rdf_graph: RDF graph object
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    """
    pathway_metadata = _get_pathway_metadata(pathway_uri, rdf_graph)

    nodes, interactions = _get_pathway_components(pathway_uri, rdf_graph)

    return convert_to_bel(nodes, interactions, pathway_metadata, hgnc_manager, chebi_manager)


def reactome_to_bel(resource_file: str, hgnc_manager, chebi_manager, export_folder=REACTOME_BEL):
    """Create Reactome BEL graphs.

    :param resource_file: rdf reactome file (there is only one)
    :param bio2bel_hgnc.Manager hgnc_manager: uniprot id to hgnc symbol dictionary
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :param export_folder: export folder
    :return:
    """
    logger.info('Parsing Reactome RDF file')
    rdf_graph = parse_rdf(resource_file, fmt='xml')

    pathways_uris_to_names = rdf_graph.query(GET_ALL_PATHWAYS, initNs=PREFIXES)

    for pathway_uri, _pathway_name in tqdm(pathways_uris_to_names, desc=f'Exporting Reactome BEL to {export_folder}'):

        # Take the identifier of the pathway which is placed at the end of the URL and also strip the number
        # next to it. (probably version of pathway)
        file_name = pathway_uri.split('/')[-1].split('.')[0]

        pickle_file = os.path.join(export_folder, f'{file_name}.pickle')

        # Skip if BEL file already exists
        if os.path.exists(pickle_file):
            continue

        bel_graph = reactome_pathway_to_bel(pathway_uri, rdf_graph, hgnc_manager, chebi_manager)

        # Export BELGraph to pickle
        to_pickle(bel_graph, pickle_file)
