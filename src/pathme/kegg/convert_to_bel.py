# -*- coding: utf-8 -*-

"""This module contains the methods to convert a KEGG RDF network into a BELGraph."""

import logging
import os
from collections import defaultdict
from itertools import product

import tqdm

from pybel import BELGraph, to_pickle
from pybel.dsl import abundance, bioprocess, complex_abundance, composite_abundance, pmod, protein, reaction
from pybel.dsl.edges import activity
from pybel.dsl.node_classes import CentralDogma
from pybel.struct import add_annotation_value
from pybel.struct.summary import count_functions, edge_summary
from .kegg_xml_parser import (
    get_all_reactions, get_all_relationships, get_complex_components, get_entity_nodes, get_reaction_pathway_edges,
    import_xml_etree,
)
from ..constants import (
    ACTIVITY_ALLOWED_MODIFIERS, CHEBI, CHEBI_NAME, HGNC, HGNC_SYMBOL, KEGG, KEGG_BEL, KEGG_CITATION, KEGG_ID,
    KEGG_MODIFICATIONS, PUBCHEM, UNIPROT,
)
from ..export_utils import add_annotation_key
from ..utils import add_bel_metadata

__all__ = [
    'get_bel_types',
    'kegg_to_bel',
    'kegg_to_pickles',
]

logger = logging.getLogger(__name__)

"""Populate empty BEL graph with KEGG pathway entities and interactions"""


def kegg_to_bel(path, hgnc_manager, chebi_manager, flatten=False):
    """Convert KGML file to a BELGraph.

    :param str path: path to KGML file
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :param bool flatten: flat nodes
    :rtype: BELGraph
    """
    xml_tree = import_xml_etree(path)  # Load xml
    root = xml_tree.getroot()

    graph = BELGraph(
        name=root.attrib['title'],
        version='1.0.0',
        description=root.attrib['link'],
        authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    add_bel_metadata(graph)

    graph.graph['pathway_id'] = root.attrib['name']

    # Parse file and get entities and interactions
    genes_dict, compounds_dict, maps_dict, orthologs_dict = get_entity_nodes(xml_tree, hgnc_manager, chebi_manager)
    relations_list = get_all_relationships(xml_tree)

    # Get compounds and reactions
    substrates_dict, products_dict = get_all_reactions(xml_tree, compounds_dict)
    reactions_dict = get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict)

    # Get complexes
    complex_ids, flattened_complexes = get_complex_components(xml_tree, genes_dict, flattened=flatten)

    # Add nodes to graph
    nodes = xml_entities_to_bel(graph, genes_dict, compounds_dict, maps_dict, flattened=flatten)
    nodes = xml_complexes_to_bel(
        graph=graph,
        node_dict=nodes,
        complex_ids=complex_ids,
        flatten_complexes=flattened_complexes if flatten else None,
    )

    # Add edges to graph
    add_edges(graph, relations_list, nodes)
    add_reaction_edges(graph, reactions_dict, nodes)

    graph.annotation_pattern['PathwayID'] = '.*'
    add_annotation_key(graph)
    add_annotation_value(graph, 'PathwayID', f'{root.attrib["org"]}{root.attrib["number"]}')

    return graph


"""Get all entities from XML tree and convert to BEL nodes"""


def xml_entities_to_bel(graph, genes_dict, compounds_dict, maps_dict, flattened=False):
    """Convert gene and compound entities in XML to BEL nodes.

    :param pybel.BELGraph graph: BEL Graph
    :param dict[str,str] genes_dict: KEGG genes (entry_id: [kegg_id, HGNC, UniProt])
    :param dict[str,str] compounds_dict: KEGG compounds (entry_id: [compound_name, ChEBI])
    :param dict[str,str] maps_dict: KEGG pathway maps (entry_id: [kegg_id, map_name])
    :param bool flattened: True to flatten to list of similar genes grouped together
    :return: KEGG entities to BEL nodes
    :rtype: dict[str,pybel.dsl.BaseEntity]
    """
    # Create a dictionary of flattened BEL nodes
    if flattened:
        node_dict = {
            node_id: flatten_gene_to_bel_node(graph, node_att)
            for node_id, node_att in genes_dict.items()
        }
        for node_id, node_att in compounds_dict.items():
            node_dict[node_id] = flatten_compound_to_bel_node(graph, node_att)

    # Create a dictionary of un-flattened BEL nodes
    else:
        node_dict = {
            node_id: gene_to_bel_node(graph, node_att)
            for node_id, node_att in genes_dict.items()
        }
        for node_id, node_att in compounds_dict.items():
            node_dict[node_id] = compound_to_bel(graph, node_att)

    for node_id, node_att in maps_dict.items():
        node_dict[node_id] = map_to_bel_node(graph, node_att)

    return node_dict


def xml_complexes_to_bel(graph, node_dict, complex_ids, flatten_complexes=None):
    """Convert complexes in XML to BEL nodes where each complex is made up of proteins and/or composites.

    :param pybel.BELGraph graph: BEL Graph
    :param dict[str,pybel.dsl.BaseEntity] node_dict: kegg_id to BEL node dictionary
    :param dict[str,list] complex_ids: complex IDs to corresponding component IDs
    :param Optional[dict[str,list]] flatten_complexes: complex IDs and flattened list of all components
    :return: kegg_ids to BEL nodes
    :rtype: dict[str,pybel.dsl.BaseEntity]
    """
    member_dict = defaultdict(list)

    if flatten_complexes is not None:
        for node_id, node_att in flatten_complexes.items():
            node_dict[node_id] = flatten_complex_to_bel_node(graph, node_att)

    # For all complexes, add BEL node component info
    else:
        for complex_id, member_ids in complex_ids.items():
            for member in member_ids:
                member_dict[complex_id].append(node_dict[member])

        for complex_id, bel_members in member_dict.items():
            node_dict[complex_id] = complexes_to_bel_node(graph, bel_members)

    return node_dict


def complexes_to_bel_node(graph, members):
    complex_node = complex_abundance(members=members)
    graph.add_node_from_data(complex_node)

    return complex_node


def gene_to_bel_node(graph, node):
    """Create a protein or protein composite BEL node and add to BEL Graph.

    :param pybel.BELGraph graph: BEL Graph
    :param list[dict[str,str]] node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    members = list()

    # Create a protein BEL node
    if len(node) == 1:
        for attribute in node:

            if HGNC in attribute:
                protein_node = protein(namespace=HGNC, name=attribute[HGNC_SYMBOL], identifier=attribute[HGNC])
                graph.add_node_from_data(protein_node)
                return protein_node

            elif UNIPROT in attribute:
                protein_node = protein(
                    namespace=UNIPROT.upper(), name=attribute[UNIPROT],
                    identifier=attribute[UNIPROT],
                )
                graph.add_node_from_data(protein_node)
                return protein_node

            else:
                protein_node = protein(namespace=KEGG.upper(), name=attribute[KEGG_ID], identifier=attribute[KEGG_ID])
                graph.add_node_from_data(protein_node)
                return protein_node

    # Create a composite abundance BEL node
    else:
        for member in node:
            bel_node = gene_to_bel_node(graph, [member])
            members.append(bel_node)

        protein_composite = composite_abundance(members=members)
        graph.add_node_from_data(protein_composite)
        return protein_composite


def flatten_gene_to_bel_node(graph, node):
    """Create a protein or list of protein BEL nodes and add to BEL Graph.

    :param pybel.BELGraph graph: BEL Graph
    :param dict[str,str] node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    # if only 1 protein node, return corresponding BEL node
    if len(node) == 1:
        node_dict = node[0]

        if HGNC in node_dict:
            protein_node = protein(namespace=HGNC, name=node_dict[HGNC_SYMBOL], identifier=node_dict[HGNC])
            graph.add_node_from_data(protein_node)
            return protein_node

        elif UNIPROT in node_dict:
            protein_node = protein(namespace=UNIPROT.upper(), name=node_dict[UNIPROT], identifier=node_dict[UNIPROT])
            graph.add_node_from_data(protein_node)
            return protein_node

        else:
            protein_node = protein(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            graph.add_node_from_data(protein_node)
            return protein_node

    proteins_list = []

    # if multiple protein nodes, return corresponding list of BEL nodes
    for node_dict in node:

        if HGNC in node_dict:
            protein_node = protein(namespace=HGNC, name=node_dict[HGNC_SYMBOL], identifier=node_dict[HGNC])
            graph.add_node_from_data(protein_node)
            proteins_list.append(protein_node)

        elif UNIPROT in node_dict:
            protein_node = protein(namespace=UNIPROT.upper(), name=node_dict[UNIPROT], identifier=node_dict[UNIPROT])
            proteins_list.append(protein_node)

        else:
            protein_node = protein(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            graph.add_node_from_data(protein_node)
            proteins_list.append(protein_node)

    return proteins_list


def compound_to_bel(graph, node):
    """Create an abundance BEL node or composite abundances BEL node and add to BEL Graph.

    :param pybel.BELGraph graph: BEL Graph
    :param dict node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    members = list()

    # Create a compound BEL node
    if len(node) == 1:
        node_dict = node[0]

        if CHEBI in node_dict:
            compound = abundance(namespace=CHEBI.upper(), name=node_dict[CHEBI_NAME], identifier=node_dict[CHEBI])
            graph.add_node_from_data(compound)
            return compound

        elif PUBCHEM in node_dict:
            compound = abundance(namespace=PUBCHEM.upper(), name=node_dict[PUBCHEM], identifier=node_dict[PUBCHEM])
            graph.add_node_from_data(compound)
            return compound

        else:
            compound = abundance(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            graph.add_node_from_data(compound)
            return compound

    # Create a composite abundance BEL node
    else:
        for member in node:
            bel_node = compound_to_bel(graph, [member])
            members.append(bel_node)

        compound_composite = composite_abundance(members=members)
        graph.add_node_from_data(compound_composite)
        return compound_composite


def flatten_compound_to_bel_node(graph, node):
    """Create an abundance or list of abundance BEL nodes and add to BEL Graph.

    :param pybel.BELGraph graph: BEL Graph
    :param dict node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    # if only 1 compound node, return corresponding BEL node
    if len(node) == 1:
        node_dict = node[0]

        if CHEBI in node_dict:

            compound = abundance(namespace=CHEBI.upper(), name=node_dict[CHEBI_NAME], identifier=node_dict[CHEBI])
            graph.add_node_from_data(compound)
            return compound

        elif PUBCHEM in node_dict:

            compound = abundance(namespace=PUBCHEM.upper(), name=node_dict[PUBCHEM], identifier=node_dict[PUBCHEM])
            graph.add_node_from_data(compound)
            return compound

        else:
            compound = abundance(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            graph.add_node_from_data(compound)
            return compound

    compounds_list = []

    # If multiple compound nodes, return flattened list of BEL nodes
    for node_dict in node:

        if CHEBI in node_dict:

            compound_node = abundance(namespace=CHEBI.upper(), name=node_dict[CHEBI_NAME], identifier=node_dict[CHEBI])
            graph.add_node_from_data(compound_node)
            compounds_list.append(compound_node)

        elif PUBCHEM in node_dict:

            compound_node = abundance(namespace=PUBCHEM.upper(), name=node_dict[PUBCHEM], identifier=node_dict[PUBCHEM])
            graph.add_node_from_data(compound_node)
            compounds_list.append(compound_node)

        else:
            compound_node = abundance(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            graph.add_node_from_data(compound_node)
            compounds_list.append(compound_node)

    return compounds_list


def map_to_bel_node(graph, node):
    """Create a biological process BEL node.

    :param pybel.BELGraph graph: BEL Graph
    :param graph: BELGraph
    :param dict node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    for attribute in node:
        name = attribute['map_name']
        identifier = attribute[KEGG_ID]

        if not name:
            logger.debug(f"KEGG API does not provide information about {node}. Using identifier.")
            name = identifier

        if name.startswith('TITLE:'):
            name = name[len('TITLE:'):]

        bio_process = bioprocess(namespace=KEGG.upper(), name=name, identifier=identifier)
        graph.add_node_from_data(bio_process)
        return bio_process


def flatten_complex_to_bel_node(graph, node):
    """Create complex abundance BEL node.

    :param pybel.BELGraph graph: BEL Graph
    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: pybel.dsl.BaseEntity
    """
    members = list()

    for node_dict in node:

        if HGNC in node_dict:
            protein_node = protein(namespace=HGNC, name=node_dict[HGNC_SYMBOL], identifier=node_dict[HGNC])
            members.append(protein_node)

        elif UNIPROT in node_dict:
            protein_node = protein(namespace=UNIPROT.upper(), name=node_dict[UNIPROT], identifier=node_dict[UNIPROT])
            members.append(protein_node)

        else:
            protein_node = protein(namespace=KEGG.upper(), name=node_dict[KEGG_ID], identifier=node_dict[KEGG_ID])
            members.append(protein_node)

    complex_members = complex_abundance(members=members)
    graph.add_node_from_data(complex_members)

    return complex_members


"""Get edges between BEL nodes"""


def add_edges(graph, edges, nodes):
    """Add edges to BEL graph.

    :param pybel.BELGraph graph: BEL Graph
    :param list[tuple] edges: list of relationships with entity IDs and interaction types
    :param dict nodes: dictionary of BEL nodes
    """
    for source, target, relation in edges:

        # Catch KeyError if entity node in list of edges is not a BEL node
        try:
            u = nodes[source]
            v = nodes[target]

        except KeyError:
            continue

        # If subject and object are lists, create edges between all products
        if isinstance(u, list) and isinstance(v, list):
            for pair in product(u, v):
                add_simple_edge(graph, pair[0], pair[1], relation)

        # If source is protein list and target is not, add edges between members in list and target
        elif isinstance(u, list) and not isinstance(v, list):
            for member in u:
                add_simple_edge(graph, member, v, relation)

        # If source is not a list and target is proteins list, add edges between them
        elif not isinstance(u, list) and isinstance(v, list):
            for member in v:
                add_simple_edge(graph, u, member, relation)

        # If entities are not lists, add edges between them
        else:
            add_simple_edge(graph, u, v, relation)


def add_reaction_edges(graph, reaction_dict, nodes):
    """Add reaction nodes and edges from reactants to products and enzymes to reactions to BEL Graph.

    :param pybel.BELGraph graph: BEL Graph
    :param dict reaction_dict: dictionary of reaction IDs and reactant and product IDs
    :param dict nodes: dictionary of BEL nodes
    """
    for k, v in reaction_dict.items():

        # Get BEL gene node(s)
        enzyme = nodes[k]

        # Get compound nodes
        for source, target, reaction_type in v:

            reactants_list = []
            products_list = []

            # Get reactant compound node
            for source_id in source:
                substrate = nodes[source_id]
                reactants_list.append(substrate)

            # Get product compound node
            for target_id in target:
                _product = nodes[target_id]
                products_list.append(_product)

                for reactant_compound in reactants_list:
                    for product_compound in products_list:

                        # If multiple compounds represent a reactant or a product, add reaction BEL nodes to graph
                        if isinstance(reactants_list, list) and isinstance(products_list, list):
                            reaction_node = reaction(reactants=reactant_compound, products=product_compound)
                            graph.add_node_from_data(reaction_node)

                        # If multiple compounds represent a reactant, add reaction BEL node to graph
                        elif isinstance(reactants_list, list) and not isinstance(products_list, list):
                            for _reactant_compound in reactants_list:
                                reaction_node = reaction(reactants=_reactant_compound, products=products_list)
                                graph.add_node_from_data(reaction_node)

                        # If multiple compounds represent a product, add reaction BEL node to graph
                        elif not isinstance(reactants_list, list) and isinstance(products_list, list):
                            for _product_compound in products_list:
                                reaction_node = reaction(reactants=reactants_list, products=_product_compound)
                                graph.add_node_from_data(reaction_node)

                        # If reactant and product is represented by a single compound, add reaction BEL node to graph
                        else:
                            reaction_node = reaction(reactants=reactants_list, products=products_list)
                            graph.add_node_from_data(reaction_node)

                # If enzyme is a list of genes, add edges between all enzymes and reactions
                if isinstance(enzyme, list):
                    for gene_type in enzyme:
                        add_simple_edge(graph, gene_type, reaction_node, reaction_type)
                else:
                    add_simple_edge(graph, enzyme, reaction_node, reaction_type)


def add_simple_edge(graph, u, v, relation_type):
    """Add corresponding edge type to BEL graph.

    :param pybel.BELGraph graph: BEL Graph
    :param u: source node
    :param v: target node
    :param list relation_type: source ID, target ID and relation types
    """
    # Check if multiple relation subtypes present
    if isinstance(relation_type, list):

        # Check if protein modification is a relation subtype
        if relation_type[1] in {'phosphorylation', 'glycosylation', 'ubiquitination', 'methylation'}:

            # If the object is a gene, miRNA, RNA, or protein, add protein modification
            if isinstance(v, CentralDogma):
                v = v.with_variants(pmod(KEGG_MODIFICATIONS[relation_type[1]]))

            # Add increases edge if pmod subtype is coupled with activation subtype
            if relation_type[0] == 'activation':
                graph.add_increases(
                    u, v,
                    citation=KEGG_CITATION, evidence='Extracted from KEGG',
                    subject_modifier=activity() if isinstance(u, ACTIVITY_ALLOWED_MODIFIERS) else None,
                    # Add the activity function if subject is one of the following nodes (BEL 2.0 specifications)
                    annotations={},
                )
                return

            # Add decreases edge if pmod subtype is coupled with inhibition subtype
            elif relation_type[0] == 'inhibition':
                graph.add_decreases(
                    u, v,
                    citation=KEGG_CITATION, evidence='Extracted from KEGG',
                    subject_modifier=activity() if isinstance(u, ACTIVITY_ALLOWED_MODIFIERS) else None,
                    # Add the activity function if subject is one of the following nodes (BEL 2.0 specifications)
                    annotations={},
                )
                return

            # Found multiple relationship which cannot be combined into one logic
            # (e.g., ['inhibition', 'indirect effect'])
            else:
                # Create all relationships in the list
                for relation in relation_type:
                    add_simple_edge(graph, u, v, relation)

                return

        # Found multiple relationship which cannot be combined into one logic (e.g., ['inhibition', 'indirect effect'])
        else:
            # Create all relationships in the list
            for relation in relation_type:
                add_simple_edge(graph, u, v, relation)

            return

    """Handle differently the relationships"""
    # If only one pmod relation subtype
    if relation_type in {'phosphorylation', 'glycosylation', 'ubiquitination', 'methylation'}:

        # If the object is a gene, miRNA, RNA, or protein, add protein modification
        if isinstance(v, CentralDogma):
            v = v.with_variants(pmod(KEGG_MODIFICATIONS[relation_type]))
        graph.add_increases(
            u, v,
            citation=KEGG_CITATION, evidence='Extracted from KEGG',
            subject_modifier=activity() if isinstance(u, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
        return

    # Subject activity decreases protein modification (i.e. dephosphorylation) of object
    elif relation_type == 'dephosphorylation':

        # If the object is a gene, miRNA, RNA, or protein, add protein modification
        if isinstance(v, CentralDogma):
            v = v.with_variants(pmod('Ph'))
        graph.add_decreases(
            u, v,
            citation=KEGG_CITATION, evidence='Extracted from KEGG',
            subject_modifier=activity() if isinstance(u, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
        return

    # Subject increases activity of object
    elif relation_type == 'activation':
        graph.add_increases(
            u, v,
            citation=KEGG_CITATION, evidence='Extracted from KEGG',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
        return

    # Catalytic activity of subject increases transformation of reactant(s) to product(s)
    elif relation_type in {'reversible', 'irreversible'}:
        graph.add_increases(
            u, v,
            citation=KEGG_CITATION, evidence='Extracted from KEGG',
            subject_modifier=activity('cat') if isinstance(u, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
        return

    # Subject decreases activity of object
    elif relation_type == 'inhibition':
        graph.add_decreases(
            u, v,
            citation=KEGG_CITATION, evidence='Extracted from KEGG',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
        return

    # Indirect effect and binding/association are noted to be equivalent relation types
    elif relation_type in {'indirect effect', 'binding/association', 'compound'}:
        graph.add_association(u, v, citation=KEGG_CITATION, evidence='Extracted from KEGG', annotations={})
        return

    # Subject increases expression of object
    elif relation_type == 'expression':

        # Expression object is converted to RNA abundance
        if isinstance(v, CentralDogma):
            v = v.get_rna()
        graph.add_increases(u, v, citation=KEGG_CITATION, evidence='Extracted from KEGG', annotations={})
        return

    # Subject decreases expression of object
    elif relation_type == 'repression':

        # Repression object is converted to RNA abundance
        if isinstance(v, CentralDogma):
            v = v.get_rna()
        graph.add_decreases(u, v, citation=KEGG_CITATION, evidence='Extracted from KEGG', annotations={})
        return

    elif relation_type in {'dissociation', 'hidden compound', 'missing interaction', 'state change'}:
        return

    raise ValueError(f'Unexpected relation type {relation_type} between {u} and {v}')


def get_bel_types(path, hgnc_manager, chebi_manager, flatten=None):
    """Get all BEL node and edge type statistics.

    :param str path: path to KGML file
    :param hgnc_manager: HGNC manager
    :type hgnc_manager: bio2bel_hgnc.Manager
    :param chebi_manager: ChEBI manager
    :type chebi_manager: bio2bel_chebi.Manager
    :param bool flatten: flat nodes
    :return: count of all nodes and edges in BEL graph
    :rtype: dict
    """
    bel_stats = {}

    bel_graph = kegg_to_bel(path, hgnc_manager, chebi_manager, flatten=True if flatten else False)

    bel_stats['nodes'] = bel_graph.number_of_nodes()
    bel_stats['edges'] = bel_graph.number_of_edges()

    # Get count of all BEL function types
    bel_functions_dict = count_functions(bel_graph)
    bel_stats.update(bel_functions_dict)

    # Get count of all BEL edge types
    bel_edges_dict = edge_summary.count_relations(bel_graph)
    bel_stats.update(bel_edges_dict)

    return bel_stats


def kegg_to_pickles(resource_files, resource_folder, hgnc_manager, chebi_manager, flatten=None, export_folder=None):
    """Export WikiPathways to Pickles.

    :param iter[str] resource_files: iterator with file names
    :param str resource_folder: path folder
    :param hgnc_manager: hgnc manager
    :param chebi_manager: ChEBI Manager
    :param bool flatten: flat nodes
    :param Optional[str] export_folder: export folder
    """
    if export_folder is None:
        export_folder = resource_folder

    for kgml_file in tqdm.tqdm(resource_files, desc=f'Exporting KEGG to BEL in {export_folder}'):
        _name = kgml_file[:-len('.xml')]
        _flatten = 'flatten' if flatten else 'unflatten'

        # Name of file created will be: "hsaXXX_unflatten.pickle" or "hsaXXX_flatten.pickle"
        pickle_path = os.path.join(
            export_folder if export_folder else KEGG_BEL,
            f'{_name}_{_flatten}.pickle',
        )

        # Skip not KGML files or file already exists
        if not kgml_file.endswith('.xml') or os.path.exists(pickle_path):
            continue

        bel_graph = kegg_to_bel(
            path=os.path.join(resource_folder, kgml_file),
            hgnc_manager=hgnc_manager,
            chebi_manager=chebi_manager,
            flatten=True if flatten else False,
        )

        to_pickle(bel_graph, pickle_path)
