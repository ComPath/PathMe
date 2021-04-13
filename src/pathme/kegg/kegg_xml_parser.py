# -*- coding: utf-8 -*-

"""This module contains functions to parse KGML files."""

import itertools as itt
import json
import logging
import os
from collections import defaultdict
from xml.etree.ElementTree import parse

import requests

from bio2bel_kegg.constants import API_KEGG_GET
from bio2bel_kegg.parsers import parse_description
from ..constants import CHEBI, CHEBI_NAME, HGNC, HGNC_SYMBOL, KEGG_CACHE, KEGG_ID, KEGG_TYPE, PUBCHEM, UNIPROT
from ..wikipathways.utils import merge_two_dicts

logger = logging.getLogger(__name__)


def import_xml_etree(filename):
    """Return XML tree from KGML file.

    :param str filename: path to KGML file
    :returns: XML Tree
    :rtype: xml.etree.ElementTree.ElementTree
    """
    try:
        tree = parse(filename)
    except IOError as ioerr:
        logger.warning('File error: %s', ioerr)
        return None

    return tree


"""KEGG Handling functions"""


def _post_process_api_query(node_meta_data, hgnc_manager, chebi_manager):
    """Process API query.

    :param dict[str,str] node_meta_data: JSON retrieved from the API
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :return: Standard identifiers for the protein/chemical
    :rtype: dict[str,str]
    """
    node_dict = {}

    if 'DBLINKS' in node_meta_data:

        for resource, identifier in node_meta_data['DBLINKS']:

            if resource not in {HGNC, UNIPROT, CHEBI, PUBCHEM}:
                continue

            # Get protein identifiers
            if resource == HGNC:
                hgnc_entry = hgnc_manager.get_gene_by_hgnc_id(identifier)

                if not hgnc_entry:
                    continue

                node_dict[HGNC] = identifier
                node_dict[HGNC_SYMBOL] = hgnc_entry.symbol

            elif resource == UNIPROT:
                hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(identifier)

                if not hgnc_entry:
                    continue

                hgnc_entry = hgnc_entry[0]  # Use the first element queried
                hgnc_id = str(hgnc_entry.identifier)
                node_dict[HGNC] = hgnc_id

                if hgnc_entry.symbol:
                    node_dict[HGNC_SYMBOL] = hgnc_entry.symbol

                else:
                    node_dict[UNIPROT] = identifier

            # Get chemical identifiers
            else:
                node_dict[resource] = identifier
                if resource == CHEBI:
                    # Split multiple identifiers and get their names
                    for chebi_id in identifier.split(' '):
                        chebi_entry = chebi_manager.get_chemical_by_chebi_id(chebi_id)

                        # If the id is found in the database stick the name
                        if chebi_entry:
                            node_dict[CHEBI_NAME] = chebi_entry.name
                        # Else use the default name by KEGG to ensure the name makes it into the graph
                        elif "ENTRY_NAME" in node_meta_data:
                            node_dict[CHEBI_NAME] = node_meta_data["ENTRY_NAME"]

    return node_dict


def _process_kegg_api_get_entity(entity, entity_type, hgnc_manager, chebi_manager):
    """Send a given entity to the KEGG API and process the results.

    :param str entity: A KEGG identifier
    :param str entity_type: Entity type
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :return: JSON retrieved from the API
    :rtype: dict[str,str]
    """
    _entity_filepath = os.path.join(KEGG_CACHE, f'{entity}.json')

    if os.path.exists(_entity_filepath):
        with open(_entity_filepath) as f:
            return json.load(f)

    kegg_url = API_KEGG_GET.format(entity)

    node_meta_data = parse_description(requests.get(kegg_url))

    node_dict = _post_process_api_query(node_meta_data, hgnc_manager, chebi_manager)

    node_dict[KEGG_ID] = entity
    node_dict[KEGG_TYPE] = entity_type

    with open(_entity_filepath, 'w') as f:
        json.dump(node_dict, f)

    return node_dict


def get_entity_nodes(tree, hgnc_manager, chebi_manager):
    """Find entry elements (KEGG pathway nodes) in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :return: genes with corresponding metadata (entry_id: [kegg_id, HGNC, UniProt])
    :return: compounds with corresponding metadata (entry_id: [compound_name, ChEBI])
    :return: biological processes with corresponding metadata  (entry_id: [kegg_id, map_name])
    :return: orthologs with corresponding metadata (entry_id: [kegg_id, kegg_type])
    :rtype: dict[str,str]
    """
    entry_dict = defaultdict(list)
    compound_dict = defaultdict(list)
    map_dict = defaultdict(list)
    ortholog_dict = defaultdict(list)

    for entry in tree.findall("entry"):

        entry_id = entry.get("id")
        kegg_ids = entry.get("name")
        kegg_type = entry.get("type")

        if kegg_type.startswith('gene'):
            for kegg_id in kegg_ids.split(' '):
                # Query the API/Cache to fetch information about protein
                node_info = _process_kegg_api_get_entity(kegg_id, kegg_type, hgnc_manager, chebi_manager)
                entry_dict[entry_id].append(node_info)

        elif kegg_type.startswith('compound'):
            for compound_id in kegg_ids.split(' '):

                compound_info = _process_kegg_api_get_entity(compound_id, kegg_type, hgnc_manager, chebi_manager)

                if compound_info:
                    compound_dict[entry_id].append(compound_info)

        elif kegg_type.startswith('map'):

            map_info = {KEGG_ID: kegg_ids}

            for graphics in entry.iter('graphics'):
                map_name = graphics.get('name')
                map_info['map_name'] = map_name

            map_dict[entry_id].append(map_info)

        elif kegg_type.startswith('ortholog'):

            for ortholog_id in kegg_ids.split(' '):
                ortholog_info = {
                    KEGG_ID: ortholog_id,
                    KEGG_TYPE: kegg_type,
                }

                ortholog_dict[entry_id].append(ortholog_info)

        # TODO: other, enzyme
        elif kegg_type.startswith('brite'):
            pass

    return entry_dict, compound_dict, map_dict, ortholog_dict


def get_complex_components(tree, genes_dict, flattened=False):
    """Get complex components to either construct complex or flatten relationships.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param dict genes_dict: dictionary of all genes in pathway
    :param bool flattened: True to flatten all complex participants
    :return: dictionary of complex IDs and component IDs (complex_id: [component_ids])
    :return: flattened dictionary of complex IDs and component metadata (complex_ids: [metadata_dict])
    :rtype: dict[str,list]
    """
    # Get IDs of complex components to construct complexes of protein composites (i.e. similar proteins).
    # or get dictionary of flattened lists of all proteins involved in complexes.
    component_info = defaultdict(list)
    complex_ids = defaultdict(list)
    complexes = defaultdict(list)
    all_components = []
    flattened_complexes = defaultdict(list)

    for entry in (tree.findall("entry")):
        entry_id = entry.get("id")

        for component in entry.iter("component"):
            component_id = component.get("id")

            # Get complex IDs and each of their component IDs
            complex_ids[entry_id].append(component_id)

            # Get the IDs of all components participating in complexes
            if component_id not in all_components:
                all_components.append(component_id)

    # Get node info for each component
    for k, v in genes_dict.items():
        for (component_id, node_info) in itt.product(all_components, v):
            if component_id == k:
                component_info[component_id].append(node_info)

    # Flatten lists of components in complexes
    if flattened:

        for k, v in complex_ids.items():
            for comp_id in v:
                for comp_k, comp_v in component_info.items():
                    if comp_id == comp_k:
                        complexes[k].append(comp_v)

        for k, v in complexes.items():
            for component in v:
                for info in component:
                    flattened_complexes[k].append(info)

    return complex_ids, flattened_complexes


def get_xml_types(tree):
    """Find entity and interaction types in KEGG XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: count of all entity, relation and reaction types present in XML
    :rtype: dict[str,int]
    """
    entity_types_dict = defaultdict(int)
    interaction_types_dict = defaultdict(int)

    for entry in tree.findall('entry'):
        entry_type = entry.get('type')

        if entry_type.startswith('gene'):
            gene_ids = entry.get('name')
            entity_types_dict['gene'] += len(gene_ids.split(' '))

        elif entry_type.startswith('ortholog'):
            ortholog_ids = entry.get('name')
            entity_types_dict['ortholog'] += len(ortholog_ids.split(' '))

        elif entry_type.startswith('compound'):
            entity_types_dict['compound entity'] += 1

        else:
            entity_types_dict[entry_type] += 1

    for relation in tree.findall('relation'):
        for subtype in relation.iter('subtype'):
            relation_subtype = subtype.get('name')
            interaction_types_dict[relation_subtype] += 1

    for reaction in tree.findall('reaction'):
        reaction_type = reaction.get('type')
        interaction_types_dict[reaction_type] += 1

    entity_types_dict['entities'] = sum(entity_types_dict.values())
    interaction_types_dict['interactions'] = sum(interaction_types_dict.values())

    return merge_two_dicts(entity_types_dict, interaction_types_dict)


"""Get all interactions in KEGG pathways"""


def get_all_relationships(tree):
    """Find all relationships between 2 entities.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: relationships list [(relation_entry1, relation_entry2, relation_subtype)]
    :rtype: list[tuple]
    """
    relations_list = []

    for relation in tree.findall("relation"):

        subtype_list = []

        relation_entry1 = relation.get("entry1")
        relation_entry2 = relation.get("entry2")
        relation_type = relation.get('type')

        for subtype in relation.iter('subtype'):

            relation_subtype = subtype.get("name")
            relation_value = subtype.get("value")
            subtype_list.append(relation_subtype)

            # TODO: assume association ??
            if not relation_subtype:
                logger.warning("No relation type declared")

            # Add protein-protein, protein-compound and transcription factor-target gene product interactions
            if relation_type in {'PPrel', 'PCrel', 'GErel'}:

                if relation_subtype == 'compound':
                    relations_list.append((relation_entry1, relation_entry2, 'binding/association'))

                else:
                    # Check if multiple relation subtypes present
                    if len(subtype_list) == 1:
                        relations_list.append((relation_entry1, relation_entry2, relation_subtype))
                    else:
                        relations_list.append((relation_entry1, relation_entry2, subtype_list))

            # Add enzyme-enzyme relations denoted as binding/association
            elif relation_type.startswith('ECrel'):
                relations_list.append((relation_entry1, relation_entry2, 'binding/association'))
                relations_list.append((relation_entry1, relation_value, 'binding/association'))
                relations_list.append((relation_value, relation_entry2, 'binding/association'))

            # Add interactions between a protein and a protein in another biological process
            elif relation_type.startswith('maplink'):
                relations_list.append((relation_entry1, relation_entry2, 'binding/association'))

    return relations_list


def get_all_reactions(tree, compounds_dict):
    """Get substrates and products with ChEBI or PubChem IDs participating in reactions.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param dict compounds_dict: dictionary of KEGG compound information
    :return: dictionary with substrate ids (reaction_id: [substrate_ids])
    :return: dictionary with product ids (reaction_id: [product_ids])
    :rtype: dict[str,list]
    """
    substrates_dict = defaultdict(list)
    products_dict = defaultdict(list)

    for reaction in tree.findall("reaction"):
        reaction_id = reaction.get("id")

        for k in compounds_dict:
            for substrate in reaction.iter('substrate'):
                substrate_id = substrate.get("id")
                if substrate_id == k:
                    substrates_dict[reaction_id].append(substrate_id)

            for product in reaction.iter('product'):
                product_id = product.get("id")
                if product_id == k:
                    products_dict[reaction_id].append(product_id)

    return substrates_dict, products_dict


def get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict):
    """Get reaction edges.

    :param xml.etree.ElementTree.ElementTree xml_tree: xml tree
    :param dict substrates_dict: dictionary with substrate info
    :param dict products_dict: dictionary with product info
    :return: dictionary of reaction elements (reaction_id: [(substrate_id, product_id, reaction_type)])
    :rtype: dict[str,list]
    """
    reactions_dict = defaultdict(list)

    for reaction in xml_tree.findall("reaction"):

        reaction_type = reaction.get("type")
        reaction_id = reaction.get("id")

        if substrates_dict[reaction_id]:
            reaction_substrates = substrates_dict[reaction_id]

            if products_dict[reaction_id]:
                reaction_products = products_dict[reaction_id]

                # Add edge from substrates to products with compound info
                reactions_dict[reaction_id].append((reaction_substrates, reaction_products, reaction_type))

                # If reaction is reversible, flip the reaction order and add a new edge
                if reaction_type == "reversible":
                    reactions_dict[reaction_id].append((reaction_products, reaction_substrates, reaction_type))

    return reactions_dict
