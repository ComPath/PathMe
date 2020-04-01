# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os
import sys

import click
import networkx as nx

import pathme.kegg.cli
import pathme.reactome.cli
import pathme.wikipathways.cli
import pybel
from pybel import to_pickle
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes, remove_isolated_list_abundances
from pybel.struct.summary import count_functions
from .constants import CX_DIR, HIPATHIA_DIR, KEGG_BEL, PPI_DIR, REACTOME_BEL, SPIA_DIR, UNIVERSE_DIR, WIKIPATHWAYS_BEL
from .export_utils import export_helper, get_universe_graph, iterate_universe_graphs

logger = logging.getLogger(__name__)

commands = dict(
    kegg=pathme.kegg.cli.main,
    wikipathways=pathme.wikipathways.cli.main,
    reactome=pathme.reactome.cli.main,
)


@click.group(help='PathMe', commands=commands)
def main():
    """Run PathMe."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


@main.group()
def export():
    """Export commands."""


kegg_path_option = click.option(
    '-k', '--kegg-path',
    help='KEGG BEL folder',
    default=KEGG_BEL,
    show_default=True,
)
reactome_path_option = click.option(
    '-r', '--reactome-path',
    help='Reactome BEL folder.',
    default=REACTOME_BEL,
    show_default=True,
)
wikipathways_path_option = click.option(
    '-w', '--wikipathways-path',
    help='WikiPathways BEL folder',
    default=WIKIPATHWAYS_BEL,
    show_default=True,
)

no_flatten_option = click.option('--no-flatten', is_flag=True, help='Do not flatten complex/reactions nodes')
no_normalize_names_option = click.option('--no-normalize-names', is_flag=True, help='Do not normalize names')


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=SPIA_DIR, show_default=True)
def spia(kegg_path, reactome_path, wikipathways_path, output):
    """Export BEL Pickles to SPIA Excel."""
    click.echo(f'Results will be exported to {output}')
    export_helper(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        output=output,
    )


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=PPI_DIR, show_default=True)
def ppi(kegg_path, reactome_path, wikipathways_path, output):
    """Export BEL Pickles to PPI-like tsv file."""
    click.echo(f'Results will be exported to {output}')
    export_helper(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        output=output,
        format='ppi',
    )


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=CX_DIR, show_default=True)
@no_flatten_option
@no_normalize_names_option
def cx(kegg_path, reactome_path, wikipathways_path, output, no_flatten, no_normalize_names):
    """Export BEL Pickles to CX."""
    try:
        from pybel_cx import to_cx_file
    except ImportError:
        click.secho('Could not import pybel_cx. Use pip install pybel-cx.')
        sys.exit(1)

    click.echo(f'Results will be exported to {output}')
    for source, path, graph in iterate_universe_graphs(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=(not no_flatten),
        normalize_names=(not no_normalize_names),
    ):
        with open(os.path.join(output, f"{path.strip('.pickle')}.cx.json"), 'w') as file:
            to_cx_file(graph, file)


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=HIPATHIA_DIR, show_default=True)
@no_flatten_option
@no_normalize_names_option
def hipathia(kegg_path, reactome_path, wikipathways_path, output, no_flatten, no_normalize_names):
    """Export to hipathia directories."""
    click.echo(f'Results will be exported to {output}')
    for source, path, graph in iterate_universe_graphs(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=(not no_flatten),
        normalize_names=(not no_normalize_names),
    ):
        pybel.to_hipathia(graph, os.path.join(output, os.path.basename(path)))


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=UNIVERSE_DIR, show_default=True)
@no_flatten_option
@no_normalize_names_option
def universe(kegg_path, reactome_path, wikipathways_path, output, no_flatten, no_normalize_names):
    """Export harmonized PathMe universe."""
    logging.basicConfig(level=logging.info, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    flatten = not no_flatten
    normalize_names = not no_normalize_names

    if not flatten:
        click.echo('Complexes and Reactions will be not be flatten to single nodes')

    if not normalize_names:
        click.echo('Names will not be normalized to lower case')

    click.echo("Merging graphs to universe and harmonizing...(this might take a while)")

    # Not explode will flip the boolean coming from the cli
    universe_graph = get_universe_graph(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names,
    )
    click.echo(f'Number of isolates after getting universe: {nx.number_of_isolates(universe_graph)}')

    # Remove isolated list abundances
    remove_isolated_list_abundances(universe_graph)

    if flatten:
        # TODO: Remove node list solo de Reactome
        click.echo(f'Number of isolates after flattening: {nx.number_of_isolates(universe_graph)}')

    click.echo("Merging variants and genes")
    collapse_all_variants(universe_graph)
    collapse_to_genes(universe_graph)
    click.echo(f'Number of isolates after collapsing variants and to genes: {nx.number_of_isolates(universe_graph)}')

    universe_graph.name = 'PathMe Universe'

    click.echo(f"Export BEL graph to: {os.path.join(output, 'pathme_universe_bel_graph.bel.pickle')}")
    click.echo(universe_graph.summary_str())
    click.echo(count_functions(universe_graph))

    to_pickle(universe_graph, os.path.join(output, "pathme_universe_bel_graph.bel.pickle"))


if __name__ == '__main__':
    main()
