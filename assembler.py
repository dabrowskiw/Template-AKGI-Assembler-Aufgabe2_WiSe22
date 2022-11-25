#! /usr/bin/python3

import sys


class Read:
    def __init__(self, lines):
        self.name = lines[0].strip()[1:]
        self.bases = "".join([x.strip() for x in lines[1:]]).upper()

    def get_kmers(self, kmersize):
        return None

    def __str__(self):
        return ""

    def __repr__(self):
        return ""

    def __eq__(self, other):
        return False


class DBGnode:
    def __init__(self, seq):
        pass

    def add_edge_to(self, eto):
        return None

    def add_edge_from(self, efrom):
        return None

    def get_potential_from(self):
        return None

    def get_potential_to(self):
        return None

    def get_edge_to_weight(self, other):
        return None

    def get_edge_from_weight(self, other):
        return None

    def extend_next(self):
        return None

    def extend_prev(self):
        return None

    def can_extend_next(self):
        return False

    def can_extend_prev(self):
        return False

class DBGraph:
    def __init__(self):
        pass

    def add_kmers(self, kmers):
        pass

    def count_edges(self):
        return 0

    def count_nodes(self):
        return 0

    def simplify(self):
        pass

    def get_FASTA(self):
        return ""

    def __str__(self):
        return ""


def read_fasta(readfile):
    return None


def build_graph(filename, kmersize):
    return None

if __name__ == "__main__":
    dbg = build_graph("data/virus_perfectreads.fasta", 2)
    print(dbg)
    dbg.simplify()
    print(dbg)
    print(dbg.get_FASTA())
