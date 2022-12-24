#!/usr/bin/env python

# encoding: utf-8


import os
import sys
import random
import logging
from operator import itemgetter
import string
import shutil
import re
import argparse
import numpy as np
from subprocess import call as subpcall
from optparse import OptionParser

#additional modules

import pybedtools #useful to sort inside python
import pyfaidx # fastest way to deal with .fasta in python

chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
          "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

def calculate_nums(fai, num, working_dir, mean, std):

    fai_file = os.path.abspath(fai)
    chrom_length_map = {}
    total_length = 0
    for line in open(fai_file):
        tmp = line.strip().split("\t")
        if tmp[0] in chroms:
            chrom_length_map[tmp[0]] = int(tmp[1])
            total_length += int(tmp[1])

    writer = open(working_dir + 'sim_chrom.config.txt', 'w')

    for k,v in chrom_length_map.items():
        count = (v / float(total_length)) * num
        out_str = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(k, int(count), mean, std, chrom_length_map[k])
        writer.write(out_str)
    writer.close()


def random_regions(working_dir):

    exclude_regions = working_dir + 'sv_repeat_telomere_centromere.bed'
    total_num = 0
    for line in open(working_dir + 'sim_chrom.config.txt', "r"):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        sv_num = closest_number(int(tmp[1]), 5)
        sv_mean = tmp[2]
        sv_std = tmp[3]
        chrom_len = tmp[4]

        chrom_dim_tsv = working_dir + chrom + ".dim.tsv"

        write_dim = open(chrom_dim_tsv, 'w')
        write_dim.write("{0}\t{1}\n".format(chrom, chrom_len))
        write_dim.close()

        tmp_file = working_dir + 'tmp.txt'
        exclude_file = working_dir + chrom + ".exclude.bed"
        write_exclude = 'echo %s > %s && grep -w -f %s %s | sortBed > %s && rm %s'%(chrom, tmp_file, tmp_file, exclude_regions, exclude_file, tmp_file)
        subpcall(write_exclude, shell=True)

        sv_out_bed = working_dir + chrom + ".random.bed"

        command = 'Rscript randomregion.r -d %s -n %s -l %s -s %s -x %s -v "deletion,tandem duplication,inverted tandem duplication,translocation cut-paste,translocation copy-paste"' \
                  ' -r "20:20:20:20:20" | sortBed > %s'%(os.path.abspath(chrom_dim_tsv), sv_num, sv_mean, sv_std, os.path.abspath(exclude_file), os.path.abspath(sv_out_bed))
        total_num += sv_num
        print("{0}, regions: {1}".format(chrom, sv_num))
        subpcall(command, shell=True)

    subpcall('cd %s && rm *.tsv *.exclude.bed'%(working_dir), shell=True)

    print("Total SV created: ", total_num)


def closest_number(n, m):
    # Find the quotient
    q = int(n / m)

    # 1st possible closest number
    n1 = m * q

    # 2nd possible closest number
    if ((n * m) > 0):
        n2 = (m * (q + 1))
    else:
        n2 = (m * (q - 1))

        # if true, then n1 is the required closest number
    if (abs(n - n1) < abs(n - n2)):
        return n1

        # else n2 is the required closest number
    return n2

def sim_genome(genome, bedfile, output, chrom):

    immutable_ref = pyfaidx.Fasta(
        os.path.abspath(genome))  # load referene, that will be used to modify real .fasta
    classic_chrs = immutable_ref.keys()  # allowed chromosomes
    possible_variants = ['deletion', 'insertion', 'inversion', 'tandem duplication', 'inverted tandem duplication',
                         'SNP', 'tandem repeat expansion', 'tandem repeat contraction', 'perfect tandem repetition',
                         'approximate tandem repetition', 'translocation cut-paste', 'translocation copy-paste',
                         'interspersed duplication', 'reciprocal translocation', 'dispersed inverted duplication', 'dispersed duplication']  # allowed variants
    valid_dna = 'ACGT'  # allowed nucleotides
    haplopattern = re.compile("^h[0-9]+$")  # allowed haplotypes for inter-haplotypes SVs are h1,h2,h3 ...

    bedlist = []

    for bed in bedfile[0]:

        if bed not in bedlist:
            bedlist.append(bed)  # remove possible duplicates but mantain input order that should be lost in set

    d = dict()  # initialize one empty dict for each haplotype

    logging.info('Organizing SVs')

    created_svs = {}

    for i, bed in enumerate(bedlist):

        if not os.path.exists(os.path.abspath(bed)):

            logging.error('.bed file ' + os.path.abspath(bed) + ' does not exist')
            exitonerror(os.path.abspath(output))

        else:

            bedh = pybedtools.BedTool(os.path.abspath(bed))

            try:

                srtbedh = bedh.sort()

            except:

                logging.error('Incorrect .bed format for .bed file ' + os.path.abspath(
                    bed) + '. Are you sure all the fields are tab-delimited?')
                exitonerror(os.path.abspath(output))

            logging.info('Organizing SVs for ' + os.path.abspath(bed))

            d["h{0}".format(i + 1)] = dict()

            for entries in srtbedh:
                sv_key = "{0}\t{1}\t{2}".format(entries[0], entries[1], entries[2])
                sv_value = (entries[0], entries[1], entries[2], entries[3], entries[4])
                created_svs[sv_key] = sv_value

                if str(entries[0]) not in classic_chrs:  # exclude errors in the first field

                    logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file ' + os.path.abspath(bed))
                    exitonerror(os.path.abspath(output))

                try:  # exclude errors in the second field

                    int(entries[1])

                except:

                    logging.error(
                        'Cannot convert ' + str(entries[1]) + ' to integer number in .bed file  ' + os.path.abspath(
                            bed) + '. Start must be an integer')
                    exitonerror(os.path.abspath(output))

                try:  # exclude errors in the third field

                    int(entries[2])

                except:

                    logging.error(
                        'Cannot convert ' + str(entries[2]) + ' to integer number in .bed file  ' + os.path.abspath(
                            bed) + '. End must be an integer')
                    exitonerror(os.path.abspath(output))

                if (int(entries[2]) - int(entries[1]) == 0):
                    logging.error('Start ' + str(entries[1]) + ' and end ' + str(
                        entries[2]) + ' cannot have the same value in .bed ' + os.path.abspath(bed) + '.')
                    exitonerror(os.path.abspath(output))

                if str(entries[3]) not in possible_variants:  # exclude errors in the fourth field

                    logging.error(str(entries[3]) + ' is not a valid variant in .bed ' + os.path.abspath(bed))
                    exitonerror(os.path.abspath(output))

                # exclude errors in the sixth field

                try:

                    int(entries[5])

                except:

                    logging.error(
                        'Cannot convert ' + str(entries[5]) + ' to integer number in .bed file ' + os.path.abspath(
                            bed) + '. Length of random sequence at breakpoint must be an integer')

                else:  # everything fine for now


                    if str(entries[3]) == 'SNP':

                        if str(entries[4]) not in valid_dna:  # information is just a valid nucleotide

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(
                                entries[3]) + '. Must be a valid DNA base included in A,C,T,G')
                            exitonerror(os.path.abspath(output))

                        if (int(entries[5])) != 0:
                            logging.warning('Incorrect length of random sequence at breakpoint ' + str(
                                entries[5]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(
                                entries[3]) + '. Coherced to 0')

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [
                                (int(entries[1]), int(entries[2]), str(entries[3]), str(entries[3]), str(entries[4]), '')]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append(
                                (int(entries[1]), int(entries[2]), str(entries[3]), str(entries[3]), str(entries[4]), ''))


                    elif str(entries[3]) == 'inversion':  # information must be None

                        if str(entries[4]) != 'None':
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[3]) + '. Must be None')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]),
                                                                         str(entries[3]), str(entries[3]), str(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             str(entries[3]), str(entries[3]), str(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))


                    elif str(entries[3]) == 'deletion':  # information must be None

                        if str(entries[4]) != 'None':
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[3]) + '. Must be None')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[3]), str(entries[4]), ''.join(random.choices(valid_dna, k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             str(entries[3]), str(entries[3]), str(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))


                    elif str(entries[3]) == 'insertion':  # information must be a valid DNA sequence

                        if not (all(i in valid_dna for i in entries[4].upper())):  # validate sequence

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(
                                entries[3]) + '. Must be a valid DNA string with A,C,T,G characters')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]),
                                                                         str(entries[3]), str(entries[3]), str(entries[4]).upper(),
                                                                         ''.join(random.choices(valid_dna,
                                                                                                k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             str(entries[3]), str(entries[3]), str(entries[4]).upper(),
                                                                             ''.join(random.choices(valid_dna, k=int(
                                                                                 entries[5])))))


                    elif str(entries[3]) == 'tandem duplication':  # information must be an integer

                        try:

                            int(entries[4])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[3]) + '. Must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]),
                                                                         str(entries[3]), str(entries[3]), int(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             str(entries[3]), str(entries[3]), int(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))



                    elif str(entries[3]) == 'inverted tandem duplication':  # information must be an integer

                        try:

                            int(entries[4])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[3]) + '. Must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]),
                                                                         str(entries[3]), str(entries[3]), int(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             str(entries[3]), str(entries[3]), int(entries[4]), ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))


                    elif str(entries[3]) == 'reciprocal translocation':

                        entr_4 = re.split('[:]', str(entries[4]))

                        if len(entr_4) != 5:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2')
                            exitonerror(os.path.abspath(output))

                        if not haplopattern.match(str(entr_4[0])):
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Haplotype must be hN, with N being any integer.')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[1]) not in classic_chrs:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Chromosome must be a valid chromosome')
                            exitonerror(os.path.abspath(output))

                        try:

                            int(entr_4[2])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Breakpoint must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[3]) not in ['forward', 'reverse'] or str(entr_4[4]) not in ['forward', 'reverse']:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Orientation1 and orientation2 must be forward or reverse')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys() and str(entr_4[4]) == 'forward':

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'del-ins', str(entries[3]),
                                                                         immutable_ref[str(entr_4[1])][
                                                                         (int(entr_4[2]) - 1) + 1:(int(entr_4[2]) - 1) + 1 + (int(entries[2]) - int(entries[1]))].seq,
                                                                         ''.join(random.choices(valid_dna,k=int(entries[5]))))]


                        elif str(entries[0]) not in d["h{0}".format(i + 1)].keys() and str(entr_4[4]) == 'reverse':

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'del-invins', str(entries[3]),
                                                                         immutable_ref[str(entr_4[1])][
                                                                         (int(entr_4[2]) - 1) + 1:(int(
                                                                             entr_4[2]) - 1) + 1 + (
                                                                                                  int(entries[2]) - int(
                                                                                                      entries[1]))].seq,
                                                                         ''.join(random.choices(valid_dna,
                                                                                                k=int(entries[5]))))]


                        elif str(entries[0]) in d["h{0}".format(i + 1)].keys() and str(entr_4[4]) == 'forward':

                            d["h{0}".format(i + 1)][str(entries[0])].append((
                                                                            int(entries[1]), int(entries[2]), 'del-ins', str(entries[3]),
                                                                            immutable_ref[str(entr_4[1])][
                                                                            (int(entr_4[2]) - 1) + 1:(int(
                                                                                entr_4[2]) - 1) + 1 + (int(
                                                                                entries[2]) - int(entries[1]))].seq,
                                                                            ''.join(random.choices(valid_dna,
                                                                                                   k=int(entries[5])))))


                        elif str(entries[0]) in d["h{0}".format(i + 1)].keys() and str(entr_4[4]) == 'reverse':

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             'del-invins', str(entries[3]),
                                                                             immutable_ref[str(entr_4[1])][
                                                                             (int(entr_4[2]) - 1) + 1:(int(
                                                                                 entr_4[2]) - 1) + 1 + (int(
                                                                                 entries[2]) - int(entries[1]))].seq,
                                                                             ''.join(random.choices(valid_dna, k=int(
                                                                                 entries[5])))))

                        if str(entr_4[0]) not in d.keys():
                            d[entr_4[0]] = dict()

                        if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) + 1, int(entr_4[2]) + 1 + (
                            int(entries[2]) - int(entries[1])), 'del-ins', str(entries[3]), immutable_ref[str(entries[0])][
                                                                           int(entries[1]) - 1:int(entries[2])].seq,
                                                                  ''.join(
                                                                      random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) + 1, int(entr_4[2]) + 1 + (
                            int(entries[2]) - int(entries[1])), 'del-invins', str(entries[3]), immutable_ref[str(entries[0])][
                                                                              int(entries[1]) - 1:int(entries[2])].seq,
                                                                  ''.join(
                                                                      random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) + 1, int(entr_4[2]) + 1 + (
                            int(entries[2]) - int(entries[1])), 'del-ins', str(entries[3]), immutable_ref[str(entries[0])][
                                                                           int(entries[1]) - 1:int(entries[2])].seq,
                                                                      ''.join(random.choices(valid_dna,
                                                                                             k=int(entries[5])))))


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) + 1, int(entr_4[2]) + 1 + (
                            int(entries[2]) - int(entries[1])), 'del-invins', str(entries[3]), immutable_ref[str(entries[0])][
                                                                              int(entries[1]) - 1:int(entries[2])].seq,
                                                                      ''.join(random.choices(valid_dna,
                                                                                             k=int(entries[5])))))


                    elif str(entries[3]) == 'translocation cut-paste':

                        entr_4 = re.split('[:]', str(entries[4]))

                        if len(entr_4) != 4:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(
                                entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation.')
                            exitonerror(os.path.abspath(output))

                        if not haplopattern.match(str(entr_4[0])):
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Haplotype must be hN, with N being any integer.')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[1]) not in classic_chrs:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Chromosome must be a valid chromosome.')
                            exitonerror(os.path.abspath(output))

                        try:

                            int(entr_4[2])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Breakpoint must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[3]) not in ['forward', 'reverse']:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Orientation must be forward or reverse')
                            exitonerror(os.path.abspath(output))

                        if str(entries[0]) not in d["h{0}".format(i + 1)].keys():

                            d["h{0}".format(i + 1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', str(entries[3]),
                                                                         'None', ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]

                        else:

                            d["h{0}".format(i + 1)][str(entries[0])].append((int(entries[1]), int(entries[2]),
                                                                             'deletion', str(entries[3]), 'None', ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))

                        if str(entr_4[0]) not in d.keys():
                            d[entr_4[0]] = dict()

                        if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                                                  immutable_ref[str(entries[0])][
                                                                  int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                                                  immutable_ref[str(entries[0])][
                                                                  int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                                                      immutable_ref[str(entries[0])][
                                                                      int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                                                     immutable_ref[str(entries[0])][
                                                                     int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                                                         random.choices(valid_dna, k=int(entries[5])))))


                    elif str(entries[3]) == 'dispersed duplication' or str(entries[3]) == 'dispersed inverted duplication':
                        entr_4 = re.split('[:]', str(entries[4]))

                        if len(entr_4) != 4:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(
                                entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation.')
                            exitonerror(os.path.abspath(output))

                        if not haplopattern.match(str(entr_4[0])):
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Haplotype must be hN, with N being any integer.')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[1]) not in classic_chrs:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Chromosome must be a valid chromosome.')
                            exitonerror(os.path.abspath(output))

                        try:

                            int(entr_4[2])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Breakpoint must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[3]) not in ['forward', 'reverse']:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Orientation must be forward or reverse')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[0]) not in d.keys():
                            d[entr_4[0]] = dict()

                        if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])] = [
                                (int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                 immutable_ref[str(entries[0])][
                                 int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                    random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])] = [
                                (int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                 immutable_ref[str(entries[0])][
                                 int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                    random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])].append(
                                (int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                 immutable_ref[str(entries[0])][int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                    random.choices(valid_dna, k=int(entries[5])))))


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])].append(
                                (int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                 immutable_ref[str(entries[0])][int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                    random.choices(valid_dna, k=int(entries[5])))))
                    else:  # is a translocation copy paste or interspersed duplication, they are the same


                        entr_4 = re.split('[:]', str(entries[4]))

                        if len(entr_4) != 4:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(
                                entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation.')
                            exitonerror(os.path.abspath(output))

                        if not haplopattern.match(str(entr_4[0])):
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Haplotype must be hN, with N being any integer.')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[1]) not in classic_chrs:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Chromosome must be a valid chromosome.')
                            exitonerror(os.path.abspath(output))

                        try:

                            int(entr_4[2])

                        except:

                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Breakpoint must be an integer')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[3]) not in ['forward', 'reverse']:
                            logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(
                                bed) + ' for variant ' + str(entries[
                                                                 3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Orientation must be forward or reverse')
                            exitonerror(os.path.abspath(output))

                        if str(entr_4[0]) not in d.keys():
                            d[entr_4[0]] = dict()

                        if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                                                  immutable_ref[str(entries[0])][
                                                                  int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                                                  immutable_ref[str(entries[0])][
                                                                  int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5]))))]


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) - 1, int(entr_4[2]), 'insertion', str(entries[3]),
                                                                      immutable_ref[str(entries[0])][int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                random.choices(valid_dna, k=int(entries[5])))))


                        elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

                            d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2]) - 1, int(entr_4[2]), 'invinsertion', str(entries[3]),
                                                                     immutable_ref[str(entries[0])][int(entries[1]) - 1:int(entries[2])].seq, ''.join(
                                                                         random.choices(valid_dna, k=int(entries[5])))))

    logging.info('SVs organized')
    logging.info('Generating .fasta haplotypes with SVs')

    # created_svs_haps = {}
    for dicts in d.keys():
        logging.info('Generating SVs for ' + str(dicts))
        parse_dict(classic_chrs, immutable_ref, d[dicts], os.path.abspath(output + '/' + chrom + '.' + str(dicts) + '.fa'))


    logging.info('Haplotypes with SVs generated')
    logging.info('Done')
    print('Done')




def parse_dict(chromosomes, fasta, dictionary, output_fasta):
    trans = str.maketrans('ATGC', 'TACG')
    skipped = 0

    created_svs = []

    for chrs in chromosomes:

        chrom = fasta[chrs]
        seq = chrom[:len(chrom)].seq

        if chrs not in dictionary.keys():  # chromosome not there, write unchanged

            write_unmodified_chromosome(chrs, seq, output_fasta)

        else:

            alterations_list = dictionary[chrs]

            if not len(alterations_list) == 1:  # else is already sorted, as it has length 1

                alterations_list = sorted(alterations_list, key=itemgetter(0, 1))

                new_alterations_list = []

                l = 0

                # Skip overlapped SVs in the alteration list
                while l < len(alterations_list):

                    start, end, typ, g_type, info, seqatbreak = alterations_list[l]
                    seqatbreak = ""  # # NOTE: we assign seqatbreak to empty to abtain single base level resolution

                    if new_alterations_list == []:

                        new_alterations_list.append((start, end, typ, g_type, info, seqatbreak))

                    else:

                        if (new_alterations_list[-1][0] <= start <= new_alterations_list[-1][1]) or (
                                new_alterations_list[-1][0] <= end <= new_alterations_list[-1][1]):

                            logging.warning('Variant with coordinates ' + chrs + ':' + str(start) + '-' + str(
                                end) + ' overlaps with variant ' + chrs + ':' + str(
                                new_alterations_list[-1][0]) + '-' + str(
                                new_alterations_list[-1][1]) + ' in the current haplotype. Skipped')
                            skipped += 1

                        else:

                            new_alterations_list.append((start, end, typ, g_type, info, seqatbreak))

                    l += 1

            else:

                new_alterations_list = alterations_list

            i = 0

            while i < len(new_alterations_list):

                start, end, typ, g_type, info, seqatbreak = new_alterations_list[i]

                seqatbreak = ""     # # NOTE: we assign seqatbreak to empty to abtain single base level resolution
                if start == 0:
                    start += 1

                if i == 0:  # first entry for the cromosome, write until the first variant start

                    seq_until_start = seq[:start - 1]

                    write_start_sequence(chrs, seq_until_start + seqatbreak, output_fasta)

                if typ == 'inversion':  # inverted sequence

                    alt_seq = Reverse(seq, start, end).translate(trans)

                    write_sequence_between(alt_seq, output_fasta)
                    created_svs.append((chrs, start, end, g_type, typ, info))

                elif typ == 'deletion':  # write nothing; deletions and insertions are also valid for translocation, as they are translated before intro insertions and deletions

                    alt_seq = ''

                    write_sequence_between(alt_seq, output_fasta)
                    created_svs.append((chrs, start, end, g_type, typ, info))

                elif typ == 'insertion':  # write specified insertion; deletions and insertions are also valid for translocation, are they are translated before intro insertions and deletions

                    alt_seq = info
                    r_type = typ
                    if g_type == "translocation copy-paste":
                        r_type = 'dispersed duplication'
                    created_svs.append((chrs, start, end + len(info), g_type, r_type, info))

                    write_sequence_between(seq[start - 1:end] + alt_seq, output_fasta)

                elif typ == 'invinsertion':

                    r_type = typ
                    if g_type == "translocation copy-paste":
                        r_type = 'inverted duplication'

                    alt_seq = info[::-1].translate(trans)
                    created_svs.append((chrs, start, end + len(info), g_type, r_type, info))

                    write_sequence_between(seq[start - 1:end] + alt_seq, output_fasta)


                elif typ == 'del-ins':

                    write_sequence_between(info, output_fasta)
                    created_svs.append((chrs, start, end, g_type, typ, info))

                elif typ == 'del-invins':

                    alt_seq = info[::-1].translate(trans)

                    write_sequence_between(alt_seq, output_fasta)
                    created_svs.append((chrs, start, end, g_type, typ, info))


                elif typ == 'tandem duplication':
                    created_svs.append((chrs, start, end, g_type, typ))
                    write_sequence_between(seq[start - 1:end] * info, output_fasta)


                elif typ == 'inverted tandem duplication':
                    created_svs.append((chrs, start, end, g_type, typ, info))
                    write_sequence_between(
                        seq[start - 1:end] + (Reverse(seq, start, end).translate(trans) * (info - 1)),
                        output_fasta)  # first part is not inverted, duplicated part it is

                elif typ == 'SNP':
                    created_svs.append((chrs, start, end, g_type, typ))
                    until_end = seq[start - 1:end - 1]

                    write_sequence_between(until_end + info, output_fasta)


                if i == len(new_alterations_list) - 1:

                    write_end_sequence(seq[end:], output_fasta)  # end not included, as it was included in the variant


                elif i < len(new_alterations_list) - 1:

                    nextstart = new_alterations_list[i + 1][0]
                    thisend = end

                    write_sequence_between(seq[thisend:nextstart - 1] + seqatbreak, output_fasta)

                i += 1

    if skipped > 0:
        logging.warning('Skipped ' + str(skipped) + ' SVs for the current haplotype as they overlapped others')

    return created_svs

def random_select_comps(comps, entry, ratio, sv_type):

    # print(comps)
    # Adding at least one event
    # num_comps = random.randint(1, 2)
    # num_comps = 2
    num_comps = len(comps)
    # print(num_comps)
    # nested sv list
    sv_list = [0 for i in range(0, num_comps + 1)]
    # start pos of each single event
    sv_starts = [0 for i in range(0, num_comps + 1)]
    # pos index, 1,2,3,... select index to insert sv
    pos_idx = [i for i in range(0, num_comps + 1)]
    # nested event name
    struct_name = [0 for i in range(0, num_comps + 1)]

    entry_idx_value = 1



    if num_comps == 2:
        sv_list[1] = entry
        pos_idx.remove(1)
        struct_name[1] = str(entry[3])
    else:
        entry_idx = random.randrange(len(pos_idx))
        entry_idx_value = pos_idx[entry_idx]
        pos_idx[entry_idx_value] = pos_idx[-1]
        del pos_idx[-1]
        sv_list[entry_idx_value] = entry
        struct_name[entry_idx_value] = entry[3]

    # add source event start pos
    sv_starts[entry_idx_value] = entry[1]

    for i in range(0, num_comps):
        added_oper = comps[0]

        if sv_type == "dispersed duplication" or sv_type == "dispersed inverted duplication" or sv_type == "tandem duplication" or sv_type == "inverted tandem duplication":

            if sv_type == 'tandem duplication' and num_comps == 1 and comps == ['tandem duplication']:
                size = int(entry[2]) - int(entry[1])
                pass
            else:

                # if num_comps == 1:
                #     lower_size = (int(entry[2]) - int(entry[1])) * 1.5
                #     upper_size = (int(entry[2]) - int(entry[1])) * 3
                # else:
                #     if added_oper == 'deletion':
                #         lower_size = (int(entry[2]) - int(entry[1])) * 3
                #         upper_size = (int(entry[2]) - int(entry[1])) * 5
                #     else:
                #         lower_size = (int(entry[2]) - int(entry[1])) * 1.5
                #         upper_size = (int(entry[2]) - int(entry[1])) * 3
                # size = random.randint(int(lower_size), int(upper_size))

                if num_comps == 1:
                    lower_size = (int(entry[2]) - int(entry[1])) * 3
                    upper_size = (int(entry[2]) - int(entry[1])) * 5
                else:
                    if added_oper == 'deletion':
                        lower_size = (int(entry[2]) - int(entry[1])) * 3
                        upper_size = (int(entry[2]) - int(entry[1])) * 5
                    else:
                        lower_size = (int(entry[2]) - int(entry[1])) * 2
                        upper_size = (int(entry[2]) - int(entry[1])) * 3
                size = random.randint(int(lower_size), int(upper_size))
        else:
            upper_size = (int(entry[2]) - int(entry[1])) * ratio
            if upper_size < 100:
                size = 100
            else:
                size = random.randint(100, int(upper_size))

        # added_oper = comps[random.randint(0, len(comps) - 1)]

        comps.remove(added_oper)
        comp_idx = random.randrange(len(pos_idx))
        comp_idx_value = pos_idx[comp_idx]

        pos_idx[comp_idx] = pos_idx[-1]
        del pos_idx[-1]
        sv_starts[comp_idx_value] = size

        comp = [added_oper, comp_idx_value, size]
        sv_list[comp_idx_value] = comp
        struct_name[comp_idx_value] = added_oper

    # print(sv_list)
    # print(entry_idx_value)
    return sv_list, entry_idx_value

def refine_trans_insert_pos(source_start, source_end, info_token):
    insert_pos = int(info_token[2])

    new_insert_pos = insert_pos

    min_dist = min(abs(insert_pos - source_start), abs(insert_pos - source_end))

    if min_dist > 20000:
        if insert_pos > source_end:
            new_insert_pos = random.randrange(source_end, source_end + 20000)
        elif insert_pos < source_start:
            new_insert_pos = random.randrange(source_start - 20000, source_start)

    new_info = '{0}:{1}:{2}:{3}'.format(info_token[0], info_token[1], new_insert_pos, info_token[3])

    return new_info


def create_mixed_bed(mixed_sv):
    starts = []
    ends = []
    chrom = ""
    var_type = ""
    comps_str = ""
    for sv in mixed_sv:
        chrom = sv[0]
        starts.append(sv[1])
        ends.append(sv[2])
        var_type += "{0};".format(sv[3])
        comp = "<{0},{1},{2}>".format(sv[3], sv[1], sv[2])
        if "translocation" in sv[3]:
            comp =  "<{0},{1},{2},{3}>".format(sv[3], sv[1], sv[2], sv[4])
        comps_str += comp

    return (chrom, min(starts), max(ends), var_type[:-1], comps_str)

import pysam
def fetch_ref_seq(ref_path, chr, start, end):
    ref = pysam.FastaFile(ref_path)
    if 'chr' not in chr:
        chr = 'chr' + chr

    ref_cutted = ref.fetch(chr, start, end)

    # print(chr, start, end, ref_path, ref)
    ref_cutted_remove_N = ref_cutted.replace('N', 'A')
    return ref_cutted_remove_N

def add_flank_region(bed, output, chrom, ratio, chrom_path, local=True):
    bedh = pybedtools.BedTool(os.path.abspath(bed))

    # choose some del to remain del rather change to del-in
    all_dels = []
    all_ddups = []
    all_iddups = []

    all_tdups = []
    all_itdups = []

    for entries in bedh:
        if str(entries[3]) == 'deletion':
            all_dels.append(entries)
        elif str(entries[3]) == 'tandem duplication':
            all_tdups.append(entries)
        elif str(entries[3]) == 'inverted tandem duplication':
            all_itdups.append(entries)
        elif str(entries[3]) == 'dispersed duplication':
            all_ddups.append(entries)
        elif str(entries[3]) == 'dispersed inverted duplication':
            all_iddups.append(entries)

    import random
    random.shuffle(all_dels)
    random.shuffle(all_ddups)
    random.shuffle(all_iddups)
    random.shuffle(all_tdups)
    random.shuffle(all_itdups)

    del_invs = all_dels

    tdup_del = all_tdups[0: int(len(all_tdups) / 2)]
    tdup_del_inv = all_tdups[int(len(all_tdups) / 2): ]
    tdup_mutil = []

    itdup_del = all_itdups[0: int(len(all_itdups) / 3)]
    itdup_del_inv = all_itdups[int(len(all_itdups) / 3): 2 * int(len(all_itdups) / 3)]


    ddup_del = all_ddups

    iddup_del = all_iddups[0: int(len(all_iddups) / 2)]

    # print("del_inv ", len(del_invs))
    # print("tdup_del ", len(tdup_del))
    # print("tdup_del_inv ", len(tdup_del_inv))
    # print("tdup_mutil ", len(tdup_mutil))
    # print("itdup_del ", len(itdup_del))
    # print("itdup_del_inv ", len(itdup_del_inv))
    # print("ddup_del ", len(ddup_del))
    # print("iddup_del ", len(iddup_del))

    # print(not_change_dels)
    new_sv_list = []
    mixed_sv_list = []

    srtbedh = bedh.sort()
    for entries in srtbedh:
        # Make simple SVs more complex
        mixed_sv = []
        sv = (entries[0], int(entries[1]), int(entries[2]), entries[3], entries[4], int(entries[5]))

        # create del-inv, del-dup
        if str(entries[3]) == 'deletion':
            # operations = ['inversion', 'tandem duplication', 'inverted tandem duplication']

            if entries not in del_invs:
                new_sv_list.append(sv)
                mixed_sv_list.append(sv)
            else:
                # print(entries)
                operations = ['inversion']
                structs, entry_idx = random_select_comps(operations, sv, ratio, "deletion")

                for comp in structs:

                    if len(comp) == 3:

                        col5 = None
                        # duplication operation
                        if 'tandem duplication' in comp[0]:
                            col5 = 2



                        if comp[1] < entry_idx:
                            added_sv = (str(entries[0]), int(entries[1]) - comp[2] + 1, int(entries[1]) - 1, comp[0], col5, 0)
                            # print('1', added_sv)

                        else:
                            added_sv = (str(entries[0]), int(entries[2]) + 1, int(entries[2]) + comp[2] + 1, comp[0], col5, 0)
                            # print('2', added_sv)

                        new_sv_list.append(added_sv)
                        mixed_sv.append(added_sv)

                    else:

                        new_sv_list.append(comp)
                        mixed_sv.append(comp)

                mixed_sv_list.append(create_mixed_bed(mixed_sv))


        elif str(entries[3]) == 'inverted tandem duplication':

            if entries in itdup_del:
                operations = ['deletion']
            elif entries in itdup_del_inv:
                operations = ['deletion', 'inversion']
            else:
                new_sv_list.append(sv)
                mixed_sv_list.append(sv)
                continue
            # operations = ['deletion', 'inversion', 'tandem duplication']
            # operations = ['deletion', 'inversion']
            new_sv_list.append(sv)
            mixed_sv.append(sv)
            structs, entry_idx = random_select_comps(operations, sv, ratio, "inverted tandem duplication")
            for comp in structs:

                if len(comp) == 3:
                    col5 = None
                    # duplication operation
                    if 'tandem duplication' in comp[0]:
                        col5 = 2

                    entries = new_sv_list[-1]

                    added_sv = (
                        str(entries[0]), int(entries[2]) + 1, int(entries[2]) + comp[2] + 1, comp[0], col5, 0)
                    # if comp[1] < entry_idx:
                    #     added_sv = (
                    #     str(entries[0]), int(entries[1]) - comp[2] + 1, int(entries[1]) - 1, comp[0], col5, 0)
                    # else:
                    #     added_sv = (
                    #     str(entries[0]), int(entries[2]) + 1, int(entries[2]) + comp[2] + 1, comp[0], col5, 0)

                    new_sv_list.append(added_sv)
                    mixed_sv.append(added_sv)
                # else:
                #     new_sv_list.append(comp)
                #     mixed_sv.append(comp)
            mixed_sv_list.append(create_mixed_bed(mixed_sv))

        elif str(entries[3]) == 'tandem duplication':

            if entries in tdup_del:
                operations = ['deletion']
            elif entries in tdup_del_inv:
                operations = ['deletion', 'inversion']
            elif entries in tdup_mutil:
                operations = ['tandem duplication']
            else:
                new_sv_list.append(sv)
                mixed_sv_list.append(sv)
                continue

            # operations = ['deletion', 'inversion', 'inverted tandem duplication']
            # operations = ['deletion', 'inversion']
            # operations = ['tandem duplication']

            operations_bak = operations.copy()
            new_sv_list.append(sv)
            mixed_sv.append(sv)

            structs, entry_idx = random_select_comps(operations, sv, ratio, "tandem duplication")
            for comp in structs:
                if len(comp) == 3:
                    # print('-', comp)

                    col5 = None
                    # duplication operation
                    if 'tandem duplication' in comp[0]:
                        col5 = 2

                    entries = new_sv_list[-1]
                    if operations_bak == ['tandem duplication']:
                        ins_seq = fetch_ref_seq(chrom_path, str(entries[0]), int(entries[1]), int(entries[2]))
                        added_sv = (str(entries[0]), int(entries[2]) + 1, int(entries[2]) + 2, 'insertion', ins_seq, 0)

                    else:
                        added_sv = (str(entries[0]), int(entries[2]) + 1, int(entries[2]) + comp[2] + 1, comp[0], col5, 0)

                    # if comp[1] < entry_idx:
                    #     added_sv = (
                    #     str(entries[0]), int(entries[1]) - comp[2] + 1, int(entries[1]) - 1, comp[0], col5, 0)
                    # else:
                    #     added_sv = (
                    #     str(entries[0]), int(entries[2]) + 1, int(entries[2]) + comp[2] + 1, comp[0], col5, 0)

                    new_sv_list.append(added_sv)
                    mixed_sv.append(added_sv)
                # else:
                #     new_sv_list.append(comp)
                #     mixed_sv.append(comp)
            mixed_sv_list.append(create_mixed_bed(mixed_sv))

        # Translocation are adjusted within 100k
        elif str(entries[3]) == "translocation copy-paste" or str(entries[3]) == "dispersed duplication" or str(entries[3]) =="dispersed inverted duplication":

            trans_info = re.split('[:]', str(entries[4]))
            new_info = str(entries[4])
            if local:
                new_info = refine_trans_insert_pos(sv[1], sv[2], trans_info)
            sv = (entries[0], int(entries[1]), int(entries[2]), entries[3], new_info, int(entries[5]))

            if entries in ddup_del:
                operations = ['deletion']
            elif entries in iddup_del:
                operations = ['deletion']
            else:
                new_sv_list.append(sv)
                mixed_sv_list.append(sv)
                continue

            new_sv_list.append(sv)
            mixed_sv.append(sv)

            new_info_tokens = re.split('[:]', new_info)
            new_insert_pos = int(new_info_tokens[2])

            source_start = int(entries[1])
            source_end = int(entries[2])

            # print(sv, new_insert_pos, source_start, source_end)
            # operations = ['deletion', 'inversion', "inverted tandem duplication", "tandem duplication"]



            # operations = ['deletion']

            structs, entry_idx = random_select_comps(operations, sv, ratio, "dispersed duplication")
            for comp in structs:
                if len(comp) == 3:
                    # print(comp)
                    col5 = None
                    # duplication operation
                    if 'tandem duplication' in comp[0]:
                        col5 = 2
                    if new_insert_pos < source_start:
                        # between source and insert pos
                        # if comp[1] < entry_idx:
                        #     added_sv = (str(entries[0]), new_insert_pos + 50, new_insert_pos + comp[2] + 50, comp[0], col5, 0)
                        # else:
                        #     added_sv = (str(entries[0]), new_insert_pos - comp[2] - 50, new_insert_pos - 50, comp[0], col5, 0)
                        added_sv = (str(entries[0]), new_insert_pos - comp[2] - 50, new_insert_pos - 50, comp[0], col5, 0)


                    elif new_insert_pos > source_end:
                        # if comp[1] < entry_idx:
                        #     added_sv = (str(entries[0]), new_insert_pos - comp[2] - 50, new_insert_pos - 50,  comp[0], col5, 0)
                        # else:
                        #     added_sv = (str(entries[0]), new_insert_pos + 50, new_insert_pos + 50 + comp[2], comp[0], col5, 0)
                        added_sv = (str(entries[0]), new_insert_pos + 50, new_insert_pos + 50 + comp[2], comp[0], col5, 0)

                    new_sv_list.append(added_sv)
                    mixed_sv.append(added_sv)


            mixed_sv_list.append(create_mixed_bed(mixed_sv))
        #
        # # Translocation are adjusted within 100k
        # elif str(entries[3]) == "translocation cut-paste":
        #
        #     trans_info = re.split('[:]', str(entries[4]))
        #     new_info = str(entries[4])
        #     if local:
        #         new_info = refine_trans_insert_pos(sv[1], sv[2], trans_info)
        #
        #     sv = (entries[0], int(entries[1]), int(entries[2]), entries[3], new_info, int(entries[5]))
        #
        #     source_start = int(entries[1])
        #     source_end = int(entries[2])
        #
        #     new_sv_list.append(sv)
        #     mixed_sv.append(sv)
        #
        #     new_info_tokens = re.split('[:]', new_info)
        #     new_insert_pos = int(new_info_tokens[2])
        #
        #     operations = ['inversion', "inverted tandem duplication", "tandem duplication"]
        #
        #     structs, entry_idx = random_select_comps(operations, sv, ratio)
        #
        #     for comp in structs:
        #         if len(comp) == 3:
        #             col5 = None
        #             # duplication operation
        #             if 'tandem duplication' in comp[0]:
        #                 col5 = 2
        #
        #             if new_insert_pos < source_start:
        #                 # between source and insert pos
        #                 if comp[1] < entry_idx:
        #                     added_sv = (
        #                     str(entries[0]), new_insert_pos + 50, new_insert_pos + comp[2] + 50, comp[0], col5, 0)
        #                 else:
        #                     added_sv = (str(entries[0]), new_insert_pos - comp[2] - 50, new_insert_pos - 50, comp[0], col5, 0)
        #
        #             elif new_insert_pos > source_end:
        #                 if comp[1] < entry_idx:
        #                     added_sv = (str(entries[0]), new_insert_pos - comp[2] - 50, new_insert_pos - 50,  comp[0], col5, 0)
        #                 else:
        #                     added_sv = (str(entries[0]), new_insert_pos + 50, new_insert_pos + 50 + comp[2], comp[0], col5, 0)
        #
        #             new_sv_list.append(added_sv)
        #             mixed_sv.append(added_sv)
        #
        #     mixed_sv_list.append(create_mixed_bed(mixed_sv))

        else:
            new_sv_list.append(sv)
            mixed_sv_list.append(sv)

    writer = open(bed.replace('.bed', '') + ".added_SVs.bed", 'w')
    for sv in new_sv_list:

        out_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(sv[0], sv[1], sv[2], sv[3], sv[4], sv[5])
        writer.write(out_str)


    mixed_writer = open(bed.replace('.bed', '') + ".added_SVs.nested.bed", 'w')
    for sv in mixed_sv_list:
        if sv[3] == 'inversion;deletion':
            svtype = 'deletion;inversion'
        else:
            svtype = sv[3]

        out_str = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(sv[0], sv[1], sv[2], svtype, sv[4])
        mixed_writer.write(out_str)


def exitonerror(output):
    print('An error occured. Check .log file at ' + os.path.abspath(output + '/VISOR_HACk.log') + ' for more details.')
    sys.exit(1)


def write_unmodified_chromosome(chromosome, seq, output_fasta):
    with open(os.path.abspath(output_fasta), 'a') as faout:
        faout.write('>' + chromosome + '\n' + seq + '\n')


def write_start_sequence(chromosome, seq, output_fasta):
    with open(os.path.abspath(output_fasta), 'a') as faout:
        faout.write('>' + chromosome + '\n' + seq)


def write_sequence_between(seq, output_fasta):
    with open(os.path.abspath(output_fasta), 'a') as faout:
        faout.write(seq)


def write_end_sequence(seq, output_fasta):
    with open(os.path.abspath(output_fasta), 'a') as faout:
        faout.write(seq + '\n')


def Change_Random_Char(word):
    length = len(word)
    word = list(word)
    k = random.sample(range(0, length), 1)
    k.sort()
    nuc_list = ['A', 'T', 'C', 'G']

    for index in k:
        add_rem = word[index]
        nuc_list.remove(add_rem)
        word[index] = ''.join(random.sample(nuc_list, k=1))
        nuc_list.append(add_rem)

    return (''.join(word))


def Delete_Random_Char(word):
    index = random.randint(0, len(word) - 1)
    word = word[:index] + word[index + 1:]

    return word


def Insert_Random_Char(word):
    nucs = ['A', 'T', 'C', 'G']

    index = random.randint(0, len(word) - 1)
    word = word[:index] + random.choice(nucs) + word[index:]

    return word


def Reverse(sequence, start, end):
    new_seq = sequence[start - 1:end][::-1]

    return new_seq


def PTR(infofield, sequence, start, end):  # new ptr

    info = re.split('[:]', infofield)
    motif, length = str(info[0]), int(info[1])
    new_seq = sequence[start - 1:end] + motif * length

    return new_seq


def ATR(infofield, sequence, start, end):  # new atr

    info = re.split('[:]', infofield)
    motif, length, altnum = str(info[0]), int(info[1]), int(info[2])
    new_seq = motif * length

    alterations = ['insertion', 'deletion', 'substitution']

    counter = 0

    while counter < altnum:

        alt_type = random.choice(alterations)

        if alt_type == 'substitution':

            new_seq = Change_Random_Char(new_seq)

        elif alt_type == 'deletion':

            new_seq = Delete_Random_Char(new_seq)

        else:  # insertion

            new_seq = Insert_Random_Char(new_seq)

        counter += 1

    return sequence[start - 1:end] + new_seq


def EXPTR(infofield, sequence, start, end):  # expand tr

    info = re.split('[:]', infofield)
    motif, num = str(info[0]), int(info[1])
    trseq = sequence[start - 1:end]

    firstbase = trseq[0]
    exprep = trseq[1:] + motif * num

    new_seq = firstbase + exprep

    return new_seq


def CTRTR(infofield, sequence, start, end):  # contract tr

    info = re.split('[:]', infofield)
    motif, num = str(info[0]), int(info[1])

    trseq = sequence[start - 1:end]

    firstbase = trseq[0]
    rep = trseq[1:]
    newind = len(motif) * num
    delrep = rep[newind:]
    new_seq = firstbase + delrep

    return new_seq

def run():
    parser = argparse.ArgumentParser(description='Simulate complex structure variants')
    subparsers = parser.add_subparsers(title='modules', dest='command',)  # two submodules

    parser_num = subparsers.add_parser('config', help='Create configuration file for simulation')
    parser_num.add_argument("-f", dest="fa", help="Chromosome file")
    parser_num.add_argument("-n", dest="num", type=int, help="Total number of SVs to simulate")
    parser_num.add_argument("-w", dest="dir", help="Working directory")
    parser_num.add_argument("-l", dest="mean", help="Average length of simulated SVs")
    parser_num.add_argument("-s", dest="std", help="Standard deviation of SV length")

    parser_random = subparsers.add_parser('random', help='Create random regions for each chromosome')
    parser_random.add_argument("-w", dest="dir", help="Working directory contains the simulation config file")


    parser_add = subparsers.add_parser('add',help='Add random simple events to current events in BED')

    parser_add.add_argument("-i", dest="input", help="Original generated from VISOR script")
    parser_add.add_argument("-w", dest="dir", help="Working dir")
    parser_add.add_argument("-c", dest="chrom", help="Chromosome to add events")
    parser_add.add_argument("-p", dest="chrom_path", help="path of chrom fa")
    parser_add.add_argument("-e", dest="local", default=True, action='store_false', help="Enable local translocation copy-paste(default=True)")
    parser_add.add_argument("-r", dest="ratio", type=float, default=0.5, help="Length ratio of SVs in flanking region (default=0.5)")

    parser_sim = subparsers.add_parser('sim_clone', help='Create abnormal genome')

    parser_sim.add_argument('-g', dest='genome', help='Template reference genome', metavar='FASTA', required=True)
    parser_sim.add_argument('-c', dest='chrom', help='Chromosome to modify', type=str, required=True)
    parser_sim.add_argument('-bed', dest='bedfile', help='One or more BED files (one for each haplotype)', metavar='BED', nargs='+', action='append', required=True)
    parser_sim.add_argument('-o', dest='output', help='Output folder', metavar='FOLDER', required=True)

    args = parser.parse_args()

    if args.command == 'config':
        calculate_nums(args.fa, args.num, args.dir, args.mean, args.std)

    elif args.command == 'random':
        random_regions(args.dir)

    elif args.command == 'add':
        add_flank_region(args.input, args.dir, args.chrom, args.ratio, args.chrom_path)

    elif args.command == 'sim_clone':

        genome_file = args.genome
        bed_file = args.bedfile
        output = args.output
        chrom = args.chrom
        sim_genome(genome_file, bed_file, output, chrom)

script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('simulate.py         Last Update:2019-11-03\n')
    print('Simulate complex structural variants\n')
    print('Usage:')
    print('simulate.py [commands] <parameters>\n')
    print('commands:')
    print('config: Create configuration file for simulation')
    print('random: Create random regions for each chromosome')
    print('add: Add random simple events to current events in BED')
    print('sim_clone: Create abnormal genome')

    print("=======================================================")
else:
    run()