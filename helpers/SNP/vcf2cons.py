#!/usr/bin/env python
import sys
import getopt
import subprocess
import errno
import shutil
import os, glob
import locale
import signal
import multiprocessing
import logging
import time
from collections import defaultdict
import tempfile
import re
from datetime import datetime

usage_message = '''
vcf2cons.py - Falk.Hildebrand@gmail.com
rewrite of vcf2cons.pl, including more advanced error model
'''

if sys.version > '3':
    long = int


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Options():
    def __init__(self):
        self.host_mem = 0.9


opt = Options()
cp = 0


def parse_opt(argv):
    try:
        opts, args = getopt.getopt(argv, "hm:o:i:",
                                   ["help",
                                    "out-dir="
                                    ])
    except getopt.error as msg:
        raise Usage('vcf2cons.py\n' + str(msg))
    if len(opts) == 0:
        raise Usage('vcf2cons.py\n' + usage_message)

    global opt
    opt.temp_dir = ""
    for option, value in opts:
        if option in ("-h", "--help"):
            print('vcf2cons.py\n' + usage_message)
            exit(0)
        elif option in ("-o", "--out-dir"):
            opt.out_dir = value + "/"
        elif option in ("-i"):
            opt.inF = value
        else:
            raise Usage("Invalid option %s", option)

    print opt.out_dir + "\n"
    if opt.temp_dir == "":
        opt.temp_dir = opt.out_dir + "tmp/"


def v2q_post_process(chr, seq, qual, gaps, l, reports, replaces, pos):
    #	for g in gaps:
    #		beg=0;
    #		if  g[0] > l:
    #			beg = g[0]-l
    #		end = g[0]+g[1]+l
    #		if (end > len(seq)):
    #			end = len(seq)

    print ">" + chr + " COV=" + str(reports) + " REPL=" + str(replaces) + " POS=" + ",".join(pos) + "\n" + seq + "\n"
    exit(0)


def main(argv=None):
    # return 3
    print sys.argv
    try:
        parse_opt(sys.argv[1:])
        depthStat = defaultdict(lambda: defaultdict(lambda: 0))
        # empty ints
        last_chr = ""; seq = "";  qual = 0;   last_pos = 0
        # arrays
        gaps = [];  spos = []
        last_pos = 0;  last_chr = ""
        _Q = 20;  _d = 1;   _D = 30000
        bcnts = 0; ccnts = 0;   lbcnts = 0;   lccnts = 0
        lcnt = 0;  chromL = 0
        regLen = re.compile('L=(\d+)=');

        if opt.inF is "-":
            data = sys.stdin.readlines()
        else:
            filO = open(opt.inF, "r")
            data = filO.readlines()
            filO.close
        print "easy", len(data), "lines\n"
        for line in data:
            lcnt += 1
            if line[0] == "#":
                continue
            t = line.split('\t')
            # print t[0]+t[2]+"\n"+line+"\n"
            # exit(0)
            if (last_chr == ""):
                last_chr = t[0]
                chromL = int(regLen.search(last_chr).group(1))  # ;  chromL = _; #}#die "$chromL\n";}
            if (last_chr != t[0]):
                if chromL > last_pos:
                    last_pos = chromL + 1
                if (last_pos > len(seq)):
                    seq += "n" * (last_pos - len(seq) - 1)
                    if (last_chr > 0):
                        v2q_post_process(last_chr, seq, qual, gaps, 5, lbcnts, lccnts, spos)
            t1 = int(t[1])
            print last_pos,"asdsa"
            if (t1 - last_pos) > 1:
                seq += 'n' * (t1 - last_pos - 1)
            if (t1 - last_pos < 0):
                print ("[vcf2cons] unsorted input\n" + t[1] + "-" + str(last_pos) + "\non line " + str(lcnt) + "\n");
                exit(0)
            regINDEL = re.compile('INDEL')
            regAZ = re.compile('^([A-Za-z.])(,[A-Za-z])*$')
            regAZres = regAZ.search(t[4])
            if (len(t[3]) == 1 and regINDEL.search(t[7]) != None and regAZres != None):
                # a SNP or reference
                ref = t[3]
                alt = regAZres.group(1)
                print ref, " ", alt, "\n"
                exit(0)
                q = 0
                b = ref
                dep = 0
                regDP = re.compile('DP=(\d+)')
                regDPres = regDP.search(t[7])
                if (regDPres != None):
                    dep = int(regDPres.group(1))
                    print t[7] + " " + str(dep) + "\n"
                    exit(0)

                if (dep <= 1):
                    b = "N"
                if ("," in alt or len(alt) > 1):
                    print line + "\n"
                    exit(3)
                regAF = re.compile('AF=(-?\d[\.\d+]?)')
                regAFres = regAF.search(t[7])

                if (regAFres != None):
                    q = regAFres.group(1)
                    if (q > 0.0501 and alt != '.' and dep > 1):
                        b = alt
                        ccnts += 1
                        lccnts += 1
                        spos.append(t[1])
                        depthStat[dep]['alt'] += 1
                        depthStat[dep]['norm'] -= 1

                depthStat[dep]['norm'] -= 1

                regQA=re.compile('QA=(\d+)')
                regQAres=regQA.search(t[7])
                if (regQAres != None and int(regQAres.group(1)) > _Q and dep >= _d and dep < _D):
                    b = b.upper()
                else :
                    b = b.lower()
                seq.append(b)
                bcnts +=1; lbcnts +=1
            elif (t[4] not in '.'):
                print line+"\n";exit(33)
                gaps.append(t[1],len(t[3]))
            last_pos = int(t[1])

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
    return 2

if __name__ == "__main__":
    sys.exit(main())
