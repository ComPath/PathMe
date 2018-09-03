# -*- coding: utf-8 -*-

"""Tests constants."""
import logging
import os

TEST_FOLDER = os.path.dirname(os.path.realpath(__file__))
KEGG_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'kegg')
WP_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'wp')
REACTOME_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'reactome')

GLYCOLYSIS_XML = os.path.join(KEGG_TEST_RESOURCES, 'hsa00010.xml')
NOTCH_XML = os.path.join(KEGG_TEST_RESOURCES, 'hsa04330.xml')

WP22 = os.path.join(WP_TEST_RESOURCES, 'WP22.ttl')

from bio2bel.testing import TemporaryConnectionMixin
from bio2bel_hgnc import Manager as HgncManager
from bio2bel_chebi import Manager as ChebiManager
from bio2bel_kegg.manager import Manager
import tempfile

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')

pathways = os.path.join(resources_path, 'hsa.txt')
protein_pathway_url = os.path.join(resources_path, 'pathway_gene.txt')

hgnc_test_path = os.path.join(resources_path, 'hgnc_test.json')

chebi_test_path = os.path.join(resources_path, 'chebi_test.tsv.gz')


class DatabaseMixin(TemporaryConnectionMixin):
    """A test case with a populated database."""

    @classmethod
    def setUpClass(cls):
        """Create temporary file"""

        """Create temporary file"""

        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path

        # create temporary database
        cls.manager = Manager(cls.connection)

        """HGNC Manager"""

        cls.hgnc_manager = HgncManager(engine=cls.manager.engine, session=cls.manager.session)
        cls.hgnc_manager.populate(hgnc_file_path=hgnc_test_path, use_hcop=False)

        """CHEBI Manager"""

        cls.chebi_manager = ChebiManager(engine=cls.manager.engine, session=cls.manager.session)
        cls.chebi_manager._populate_compounds(url=chebi_test_path)

    @classmethod
    def tearDownClass(cls):
        """Close the connection in the manager and deletes the temporary database."""
        cls.manager.drop_all()
        cls.hgnc_manager.drop_all()
        cls.manager.session.close()
        cls.hgnc_manager.session.close()
        super().tearDownClass()