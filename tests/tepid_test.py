#! /usr/bin/env python

from tepid import tepid
import os
from glob import glob
import random
import string

def test():

    class Arg():

        def __init__(self):
            return

    os.mkdir("./tests/run_tests")
    os.chdir("./tests/run_tests")
    # print os.path.dirname(os.path.realpath(__file__))

    discover = Arg()
    discover.keep = True
    discover.deletions = False
    discover.insertions = False
    discover.strict = False
    discover.mask = ""
    discover.discordant = False
    discover.proc = 1
    discover.name = "run"
    discover.conc =  "../conc.bam"
    discover.split = "../split.bam"
    discover.te = "../../Annotation/Arabidopsis/TAIR9_TE.bed.gz"
    discover.prefix = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + "_"

    tepid.discover_pe(discover)
    test_ins = [a for a in open("./insertions_run.bed")]
    true_ins = [a for a in open("../insertions_test.bed")]
    test_del = [a for a in open("./deletions_run.bed")]
    true_del = [a for a in open("../deletions_test.bed")]
    assert test_ins == true_ins
    assert test_del == true_del

    empty = Arg()
    empty.keep = True
    empty.deletions = False
    empty.insertions = False
    empty.strict = False
    empty.mask = ""
    empty.discordant = False
    empty.proc = 1
    empty.name = "run"
    empty.conc =  "../conc_empty.bam"
    empty.split = "../split_empty.bam"
    empty.te = "../../Annotation/Arabidopsis/TAIR9_TE.bed.gz"
    empty.prefix = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + "_"

    assert tepid.discover_pe(empty) == 1, "Empty output test failed"


def tearDown():
    temp = glob('./*')
    for i in temp:
        os.remove(i)
    os.chdir("../..")
    os.rmdir("./tests/run_tests")
