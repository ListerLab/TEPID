#! /usr/bin/env python

from tepid import tepid
import os
from glob import glob


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

    tepid.discover_pe(discover)
    test_ins = [a for a in open("./insertions_run.bed")]
    true_ins = [a for a in open("../insertions_test.bed")]
    test_del = [a for a in open("./deletions_run.bed")]
    true_del = [a for a in open("../deletions_test.bed")]
    assert test_ins == true_ins
    assert test_del == true_del


def tearDown():
    temp = glob('./*')
    for i in temp:
        os.remove(i)
    os.chdir("../..")
    os.rmdir("./tests/run_tests")