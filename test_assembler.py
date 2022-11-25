from unittest import TestCase
import pytest

from assembler import Read, DBGnode, DBGraph
import assembler

class AssemblerTestCase(TestCase):
    @pytest.mark.read
    def test_read_constructor(self):
        read = Read([">Read 1", "AGTCGTAGC"])
        self.assertEqual("Read 1", read.name)
        self.assertEqual("AGTCGTAGC", read.bases)

    @pytest.mark.read
    def test_read_str(self):
        read = Read([">Read 2", "AGTCGTAGCGTACCGTAGCCCGTTGGCAGTA"])
        self.assertEqual("Read 2: AGTCGTAGCGTACCGTAGCC...", str(read))

    def test_read_repr(self):
        read = Read([">Read 3", "AGTCGTAGCGTGTCCGACGTAGCCCGTTGGCAGTA"])
        self.assertEqual(str(read), read.__repr__())
        read = Read([">Read 3", "CCT"])
        self.assertEqual(str(read), read.__repr__())

    @pytest.mark.read
    def test_read_get_kmers(self):
        read = Read([">Read", "AGTCGTAGTC"])
        self.assertEqual({'AGTC': 2, 'GTCG': 1, 'TCGT': 1, 'CGTA': 1, 'GTAG': 1, 'TAGT': 1}, read.get_kmers(4))
        read = Read([">Read", "AGT"])
        self.assertEqual({}, read.get_kmers(4))
        read = Read([">Read", "AGTCGTACCGTGACGTACGAT"])
        self.assertEqual({'A': 5, 'G': 6, 'T': 5, 'C': 5}, read.get_kmers(1))

    @pytest.mark.read
    def test_read_eq(self):
        read1 = Read([">Read", "AGTCGTAGTC"])
        read2 = Read([">Read", "AGTCGTAGTC"])
        read3 = Read([">Read", "AGTCGTGGTC"])
        self.assertTrue(read1 == read2)
        self.assertFalse(read1 == read3)
        self.assertFalse(read1 == "test")

    @pytest.mark.dbgnode
    def test_dbgnode_get_potential_from(self):
        node = DBGnode("AGTC")
        self.assertEqual(set(["AAGT", "GAGT", "TAGT", "CAGT"]), set(node.get_potential_from()))

    @pytest.mark.dbgnode
    def test_dbgnode_get_potential_to(self):
        node = DBGnode("TAC")
        self.assertEqual(set(["ACA", "ACG", "ACT", "ACC"]), set(node.get_potential_to()))

    @pytest.mark.dbgnode
    def test_dbgnode_add_edge_to(self):
        node1 = DBGnode("TAC")
        node2 = DBGnode("ACA")
        node3 = DBGnode("ACG")
        node1.add_edge_to(node2)
        self.assertEqual(0, node1.get_edge_to_weight(node1))
        self.assertEqual(1, node1.get_edge_to_weight(node2))
        self.assertEqual(0, node1.get_edge_to_weight(node3))
        node1.add_edge_to(node2)
        self.assertEqual(0, node1.get_edge_to_weight(node1))
        self.assertEqual(2, node1.get_edge_to_weight(node2))
        self.assertEqual(0, node1.get_edge_to_weight(node3))
        node1.add_edge_to(node3)
        self.assertEqual(0, node1.get_edge_to_weight(node1))
        self.assertEqual(2, node1.get_edge_to_weight(node2))
        self.assertEqual(1, node1.get_edge_to_weight(node3))
        node1.add_edge_to(node1)
        self.assertEqual(1, node1.get_edge_to_weight(node1))
        self.assertEqual(2, node1.get_edge_to_weight(node2))
        self.assertEqual(1, node1.get_edge_to_weight(node3))

    @pytest.mark.dbgnode
    def test_dbgnode_add_edge_from(self):
        node1 = DBGnode("TAC")
        node2 = DBGnode("ATA")
        node3 = DBGnode("GTA")
        node1.add_edge_from(node2)
        self.assertEqual(0, node1.get_edge_from_weight(node1))
        self.assertEqual(1, node1.get_edge_from_weight(node2))
        self.assertEqual(0, node1.get_edge_from_weight(node3))
        node1.add_edge_from(node2)
        self.assertEqual(0, node1.get_edge_from_weight(node1))
        self.assertEqual(2, node1.get_edge_from_weight(node2))
        self.assertEqual(0, node1.get_edge_from_weight(node3))
        node1.add_edge_from(node3)
        self.assertEqual(0, node1.get_edge_from_weight(node1))
        self.assertEqual(2, node1.get_edge_from_weight(node2))
        self.assertEqual(1, node1.get_edge_from_weight(node3))
        node1.add_edge_from(node1)
        self.assertEqual(1, node1.get_edge_from_weight(node1))
        self.assertEqual(2, node1.get_edge_from_weight(node2))
        self.assertEqual(1, node1.get_edge_from_weight(node3))

    @pytest.mark.dbgraph
    def test_dbgraph_correct(self):
        kmers = {"AGT": 2, "GTC": 1, "TCA": 2, "GTG": 3}
        graph = DBGraph()
        graph.add_kmers(kmers)
        self.assertEqual(3, graph.count_edges())

    @pytest.mark.dbgraph
    def test_dbgraph_wrongkmer(self):
        kmers = {"AGT": 2, "GTC": 1, "TCAA": 2, "GTG": 3}
        graph = DBGraph()
        with pytest.raises(ValueError):
            graph.add_kmers(kmers)

    @pytest.mark.toplevel
    def test_read_fasta(self):
        expected = [Read([">1", "AGTCTAC"]), Read([">2", "TCTACCG"])]
        self.assertEqual(expected, assembler.read_fasta("data/test.fasta"))

    @pytest.mark.toplevel
    def test_build_graph(self):
        graph = assembler.build_graph("data/test.fasta", 2)
        self.assertEqual(8, graph.count_nodes())
        self.assertEqual(17, graph.count_edges())
        graph = assembler.build_graph("data/test.fasta", 4)
        self.assertEqual(6, graph.count_nodes())
        self.assertEqual(5, graph.count_edges())
        graph = assembler.build_graph("data/virus_perfectreads.fasta", 6)
        self.assertEqual(1581, graph.count_nodes())
        self.assertEqual(3644, graph.count_edges())
        graph = assembler.build_graph("data/virus_perfectreads.fasta", 12)
        self.assertEqual(2418, graph.count_nodes())
        self.assertEqual(2418, graph.count_edges())

    @pytest.mark.dbgraph
    def test_dbgraph_simplify(self):
        graph = DBGraph()
        graph.add_kmers(Read([">read", "ATGCGTAGC"]).get_kmers(3))
        self.assertEqual(7, graph.count_edges())
        self.assertEqual(7, graph.count_nodes())
        graph.simplify()
        self.assertEqual(2, graph.count_edges())
        self.assertEqual(2, graph.count_nodes())

    @pytest.mark.dbgraph
    def test_dbgraph_getFASTA(self):
        graph = DBGraph()
        graph.add_kmers(Read([">read", "ATGCGTAGC"]).get_kmers(3))
        graph.simplify()
        fasta = graph.get_FASTA()
        seqsInFasta = ["".join([y.strip() for y in x.split("\n")[1:]]) for x in fasta.split(">")[1:]]
        self.assertTrue(len(seqsInFasta) == 2)
        self.assertTrue("ATGC" in seqsInFasta, "Expected sequence 'ATGCC' in returned FASTA: \n" + fasta)
        self.assertTrue("GCGTAGC" in seqsInFasta, "Expected sequence 'GCGTAGC' in returned FASTA: \n" + fasta)

    @pytest.mark.dbgnode
    def test_dbgnode_can_extend_next(self):
        node1 = DBGnode("AGA")
        node2 = DBGnode("GGT")
        node3 = DBGnode("GTA")
        node4 = DBGnode("TAG")
        node5 = DBGnode("AGC")
        node1.add_edge_to(node2)
        node2.add_edge_to(node3)
        node3.add_edge_to(node4)
        node4.add_edge_to(node1)
        node4.add_edge_to(node5)
        node1.add_edge_from(node4)
        node2.add_edge_from(node1)
        node3.add_edge_from(node2)
        node4.add_edge_from(node3)
        node5.add_edge_from(node4)
        self.assertTrue(node1.can_extend_next())
        self.assertTrue(node2.can_extend_next())
        self.assertTrue(node3.can_extend_next())
        self.assertFalse(node4.can_extend_next())
        self.assertFalse(node5.can_extend_next())


    @pytest.mark.dbgnode
    def test_dbgnode_can_extend_prev(self):
        node1 = DBGnode("AGA")
        node2 = DBGnode("GGT")
        node3 = DBGnode("GTA")
        node4 = DBGnode("TAG")
        node5 = DBGnode("AGC")
        node1.add_edge_to(node2)
        node2.add_edge_to(node3)
        node3.add_edge_to(node4)
        node4.add_edge_to(node1)
        node4.add_edge_to(node5)
        node1.add_edge_from(node4)
        node2.add_edge_from(node1)
        node3.add_edge_from(node2)
        node4.add_edge_from(node3)
        node5.add_edge_from(node4)
        self.assertFalse(node1.can_extend_prev())
        self.assertTrue(node2.can_extend_prev())
        self.assertTrue(node3.can_extend_prev())
        self.assertTrue(node4.can_extend_prev())
        self.assertFalse(node5.can_extend_prev())

    @pytest.mark.dbgnode
    def test_dbgnode_extend_next(self):
        node1 = DBGnode("AGA")
        node2 = DBGnode("GGT")
        node3 = DBGnode("GTA")
        node1.add_edge_to(node2)
        node2.add_edge_to(node3)
        node2.add_edge_from(node1)
        node3.add_edge_from(node2)
        self.assertEqual(1, node1.get_edge_to_weight(node2))
        self.assertEqual(0, node1.get_edge_to_weight(node3))
        node1.extend_next()
        self.assertEqual(1, node1.get_edge_to_weight(node3))
        self.assertEqual(0, node1.get_edge_to_weight(node2))

    @pytest.mark.dbgnode
    def test_dbgnode_extend_prev(self):
        node1 = DBGnode("AGA")
        node2 = DBGnode("GGT")
        node3 = DBGnode("GTA")
        node1.add_edge_to(node2)
        node2.add_edge_to(node3)
        node2.add_edge_from(node1)
        node3.add_edge_from(node2)
        self.assertEqual(1, node3.get_edge_from_weight(node2))
        self.assertEqual(0, node3.get_edge_from_weight(node1))
        node3.extend_prev()
        self.assertEqual(1, node3.get_edge_from_weight(node1))
        self.assertEqual(0, node3.get_edge_from_weight(node2))
