# -*- coding: utf-8 -*-

"""This module contains the methods that run SPARQL queries to convert the Reactome Pathways to BEL."""

import logging
from collections import defaultdict
from typing import Set, Dict, Union, Tuple, List, Any

import rdflib
import tqdm
from rdflib import URIRef
from rdflib.namespace import Namespace, RDFS, RDF, DCTERMS, DC, OWL, XSD, SKOS

from compath_reloaded.utils import query_result_to_dict, parse_rdf, get_pathway_statitics
from .convert_to_bel import convert_to_bel

log = logging.getLogger(__name__)

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

#: SPARQL query string to get all the  primary types of entries (Pathway,  BiochemicalReaction) in a pathway network.
GET_ALL_TYPES = """
SELECT DISTINCT (STRAFTER(STR(?rdf_type), str(biopax3:)) AS ?entry_type)
WHERE
    {
        ?uri_id rdf:type ?rdf_type .
}
"""

#: SPARQL query string to get all pathways in a RDF network.
GET_ALL_PATHWAYS = """
SELECT DISTINCT ?uri_id ?name
WHERE
    {
        ?uri_id rdf:type biopax3:Pathway .
        ?uri_id biopax3:displayName ?name .
}
"""

#: SPARQL query string to get all components of a pathway (predicate biopax3:pathwayComponent).
GET_ALL_PATHWAY_COMPONENTS = """
SELECT DISTINCT ?uri_id ?name ?comment (STRAFTER(STR(?uri_type), str(biopax3:)) AS ?component_type)
WHERE
    {
        ?pathway biopax3:pathwayComponent ?uri_id .
        ?uri_id rdf:type ?uri_type .
        optional {?uri_id biopax3:displayName ?name .}
        optional {?uri_id biopax3:comment ?comment .}
    }
"""

#: SPARQL query string to get all participants in an interaction and its controlType (ACTIVATION or INHIBITION).
GET_INTERACTION_PARTICIPANTS_AND_TYPE = """
SELECT DISTINCT (STRAFTER(STR(?component), '#') AS ?identifier) ?reactant ?product (STR(?control_type) AS ?interaction_type)
WHERE
    {
        ?component biopax3:left ?reactant .
        ?component biopax3:right ?product .
        optional {?control biopax3:controlled ?component .}
        optional {?control biopax3:controlType ?control_type }
    }
"""

#: SPARQL query to get all the possible metadate (optional statements) of an entity (Protein, Dna, Pathway...).
GET_ENTITY_METADATA = """
SELECT DISTINCT (STRAFTER(STR(?entity), '#') AS ?identifier) (STR(?entity) AS ?uri_reactome_id) ?uri_id ?name ?identifier ?cell_locat ?display_name ?complex_components ?comment (STRAFTER(STR(?uri_type), str(biopax3:)) AS ?entity_type) 
WHERE
    {        
        ?entity rdf:type ?uri_type .
        ?entity biopax3:comment ?comment .
        
        optional {?entity biopax3:entityReference ?uri_id .}
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


def _get_pathway_metadata(pathway_uri: rdflib.URIRef, rdf_graph: rdflib.Graph) -> Dict[str, Union[str, Set[str]]]:
    """Get metadata for a pathway entry.

    :param pathway_uri: URI reference of the queried graph
    :param rdf_graph: RDF Reactome Universe graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    return query_result_to_dict(
        rdf_graph.query(
            GET_ENTITY_METADATA,
            initNs=PREFIXES,
            initBindings={'entity': pathway_uri}
        ),
        attr_empty=['display_name', 'identifier', 'uri_reactome_id', 'comment'],
    )


def _get_entity_metadata(entity: rdflib.URIRef, rdf_graph: rdflib.Graph) -> Dict[str, Union[str, Set[str]]]:
    """Get the metadata for an entity (Protein, Dna, Complex...).

    :param entity: URI reference of the queried entity
    :param rdf_graph: RDF Reactome Universe graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    entity_metadata = query_result_to_dict(rdf_graph.query(
        GET_ENTITY_METADATA,
        initNs=PREFIXES,
        initBindings={'entity': entity}
    ),
        attr_empty=['entity_type']
    )

    # Complexes might contain multiple components entities so we iterate over the complex components to fetch that information

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


def _get_reaction_participants(component_uri: rdflib.URIRef, component, rdf_graph: rdflib.Graph) -> Tuple[
    Dict[Union[str, Set[str]], Dict[str, Union[str, Set[str]]]], Dict[Any, Dict[str, Any]]]:
    """Get reaction participants (nodes and interactions) for a given reaction.

    :param component_uri: URI reference of the queried reaction component
    :param component: Reaction component metadata
    :param rdf_graph: RDF Reactome Universe graph object
    :return: returns the reaction participants as entities (Proteins, Complex, SmallMolecule...) and proteins (the reaction link)
    """
    interactions = {}
    nodes = {}

    spaqrl_reaction_participants = rdf_graph.query(
        GET_INTERACTION_PARTICIPANTS_AND_TYPE,
        initNs=PREFIXES,
        initBindings={'component': component_uri}
    )

    for interaction in spaqrl_reaction_participants:

        if not 'interaction_type' in interaction.labels:
            interaction_type = 'unknown'
        else:
            interaction_type = interaction.interaction_type

        if interaction.identifier not in interactions.keys():
            interactions[interaction.identifier] = {'metadata': component
                                                    }
        reactant_metadata = _get_entity_metadata(interaction.reactant, rdf_graph)
        product_metadata = _get_entity_metadata(interaction.product, rdf_graph)

        reactant_id = reactant_metadata['identifier']
        product_id = product_metadata['identifier']

        nodes[reactant_id] = reactant_metadata
        nodes[product_id] = product_metadata

        if 'participants' not in interactions[interaction.identifier].keys():
            interactions[interaction.identifier]['participants'] = (reactant_id, product_id)

        else:
            interaction_participants = interactions[interaction.identifier]['participants']
            if isinstance(interaction_participants, tuple):
                interactions[interaction.identifier]['participants'] = {
                    'reactants': {interaction_participants[0], reactant_id},
                    'products': {interaction_participants[1], product_id}
                }
            else:
                interactions[interaction.identifier]['participants']['reactants'].add(reactant_id)
                interactions[interaction.identifier]['participants']['products'].add(product_id)

        interactions[interaction.identifier]['metadata']['interaction_type'] = str(interaction_type)

    return nodes, interactions


def _get_pathway_components(pathway_uri: rdflib.URIRef, rdf_graph: rdflib.Graph) -> Tuple[
    Dict[str, Dict[str, Union[str, Set[str]]]], List[Dict[str, Union[str, Set[str]]]]]:
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
        initBindings={'pathway': pathway_uri}
    )

    pathway_components = query_result_to_dict(spaqrl_pathway_components, id_dict=True)

    for component_uri, component in pathway_components.items():
        if component['component_type'] == 'BiochemicalReaction':
            component_nodes, component_interactions = _get_reaction_participants(component_uri, component, rdf_graph)

            nodes.update(component_nodes)
            interactions.update(component_interactions)

        elif component['component_type'] == 'Pathway':
            pathway_metadata = _get_pathway_metadata(component_uri, rdf_graph)
            nodes[pathway_metadata['uri_reactome_id']] = pathway_metadata

    return nodes, list(interactions.values())


def get_reactome_statistics(resource_file, hgnc_manager):
    """Get types statistics for Reactome.

    :param str rdf_graph: primary entries type identifier (ex: DataNode or Interaction)
    :param str primary_type: primary entries type identifier (ex: DataNode or Interaction)
    """
    log.info('Parsing Reactome RDF file')
    rdf_graph = parse_rdf(resource_file, fmt='xml')

    spaqrl_all_pathways = rdf_graph.query(GET_ALL_PATHWAYS, initNs=PREFIXES)

    global_statistics = defaultdict(lambda: defaultdict(int))

    for pathway_uri, pathway_title in tqdm.tqdm(spaqrl_all_pathways, desc='Generating Reactome Statistics'):
        print(pathway_title)

        nodes, edges = _get_pathway_components(pathway_uri, rdf_graph)
        pathway_metadata = _get_pathway_metadata(pathway_uri, rdf_graph)

        nodes_types = [
            node['entity_type'] for node in nodes.values()
        ]
        edges_types = [
            edge['metadata']['interaction_type'] for edge in edges
        ]

        bel_graph = convert_to_bel(nodes, edges, pathway_metadata, hgnc_manager)

        global_statistics, pathway_statistics = get_pathway_statitics(nodes_types, edges_types, bel_graph,
                                                                      global_statistics=global_statistics)

        print(pathway_statistics)

    return global_statistics


def rdf_pathway_to_bel(pathway_uri, rdf_graph, hgnc_manager):
    """Convert WikiPathways RDF file to BEL.

    :param str filepath: path to the file
    :rtype: pybel.BELGraph
    :param bio2bel_hgnc.Manager hgnc_manager: Bio2BEL HGNC Manager
    """
    pathway_metadata = _get_pathway_metadata(pathway_uri, rdf_graph)

    nodes, interactions = _get_pathway_components(pathway_uri, rdf_graph)

    return convert_to_bel(nodes, interactions, pathway_metadata, hgnc_manager)


def reactome_to_pybel(resource_file, hgnc_manager):
    """Load Reactome graphs.

    :param str resource_file: rdf reactome file (there is only one)
    :param dict bio2bel_hgnc.Manager: uniprot id to hgnc symbol dictionary
    :return:
    """
    log.info('Parsing Reactome RDF file')
    rdf_graph = parse_rdf(resource_file, fmt='xml')

    reactome_pathways = []

    spaqrl_all_pathways = rdf_graph.query(GET_ALL_PATHWAYS, initNs=PREFIXES)

    for pathway_uri, pathway_title in tqdm.tqdm(spaqrl_all_pathways, desc='Creating Reactome BELGraphs'):
        # debug_pathway_info(bel_graph, pathway_path)

        reactome_pathways.append(rdf_pathway_to_bel(pathway_uri, rdf_graph, hgnc_manager))

    return reactome_pathways