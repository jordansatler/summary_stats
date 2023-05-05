#!/usr/bin/env python

"""
Calculate a set of summary statistics for multilocus sequence data.
Summary stats calculated include nucleotide diversity (pi),
number of segregating sites (SS), Watterson's theta (per locus),
Watterson's theta (per site), and Tajima's D. Requires a traits file
(individuals are assigned to species) and a folder of nexus formatted loci.

usage:
    python sumstats.py traits.file /path/to/loci/
"""

import os
import sys
import glob
import numpy as np
import shutil
import dendropy
from dendropy.calculate import popgenstat
from Bio import SeqIO

def traits(infile):
    """read traits file and assign individuals to OTUs"""
    with open(infile, 'r') as t:
        return {i.split()[0]:i.split()[1] for i in t}

def locus_files(files):
    """get list of nexus files"""
    return glob.glob(files + "/*.nex*")

def read_file(traits, locus):
    """read nexus file and filter if meets sampling requirement"""
    # record if OTU is represented in locus
    partitions = {sp:{} for sp in set(traits.values())}
    with open(locus, 'r') as l:
        records = [rec for rec in SeqIO.parse(locus, "nexus") 
                               if rec.id in traits]
        for i in records:
            partitions[traits[i.id]][i.id] = str(i.seq)
        return partitions

def get_character_matrix(partitions):
    """convert dicts to character matrices and calculate pop get stats"""
    sumstats = {}
    for sp, data in partitions.items():
        # require at least two alleles per species per locus
        if len(data) < 2:
            continue

        newD = dendropy.DnaCharacterMatrix.from_dict(data)
        # test if Tajima's D can be calculated
        try:
            popgenstat.tajimas_d(newD)
        except ZeroDivisionError:
            # print("Tajima's D gave a zero division error!")
            sumstats[sp] = [popgenstat.nucleotide_diversity(newD),
                            popgenstat.num_segregating_sites(newD),
                            popgenstat.wattersons_theta(newD),
                            popgenstat.wattersons_theta(newD) / len(newD[0])]
            continue
        # return all four SS if Tajima's D calculated
        sumstats[sp] = [popgenstat.nucleotide_diversity(newD),
                        popgenstat.num_segregating_sites(newD),
                        popgenstat.wattersons_theta(newD),
                        popgenstat.wattersons_theta(newD) / len(newD[0]),
                        popgenstat.tajimas_d(newD)]
    return sumstats

def put_ss_together(ssList):
    """place summary statistics into separate lists"""
    vals = {}
    for k, v in ssList.items():
        pi = [l[0] for l in v]
        S = [l[1] for l in v]
        Theta_per_locus = [l[2] for l in v]
        Theta_per_site = [l[3] for l in v]
        D = [l[4] for l in v if len(l) == 5]
        # add SS to new dictionary
        vals[k] = [pi, S, Theta_per_locus, Theta_per_site, D]
    return vals

def write_ss(ss_pop):
    """write summary statistics to a file"""
    with open("results_sumstats.txt", 'w') as out:
        out.write("### Infiles ###\n\nLoci:\t{0}\nTraits:\t{1}\n".format(sys.argv[2], 
                                                                           sys.argv[1]))
        out.write("\n### Species ###")
        for sp, vals in sorted(ss_pop.items()):
            out.write('\n\n*** ' + sp + ' ***\n\n')
            for i in range(len(vals)):
                # make sure summary stat was calculated
                if len(vals[i]) == 0:
                    continue

                if i == 0:
                    out.write("pi\tN\t{0}\tmean\t{1:.4f}\tsd\t{2:.4f}\n".format(
                               len(vals[i]), np.mean(vals[i]), np.std(vals[i])))
                elif i == 1:
                    out.write("SS\tN\t{0}\tmean\t{1:.4f}\tsd\t{2:.4f}\n".format(
                               len(vals[i]), np.mean(vals[i]), np.std(vals[i])))
                elif i == 2:
                    out.write("Watterson's theta_locus\tN\t{0}\tmean\t{1:.4f}\tsd\t{2:.4f}\n".format(
                               len(vals[i]), np.mean(vals[i]), np.std(vals[i])))
                elif i == 3:
                    out.write("Watterson's theta_site\tN\t{0}\tmean\t{1:.4f}\tsd\t{2:.4f}\n".format(
                               len(vals[i]), np.mean(vals[i]), np.std(vals[i])))
                else:
                    out.write("Tajima's D\tN\t{0}\tmean\t{1:.4f}\tsd\t{2:.4f}\n".format(
                               len(vals[i]), np.mean(vals[i]), np.std(vals[i])))

def move_file():
    """move fasta files to results folder"""
    if os.path.exists('../results/'):
        print("\nresults folder already exists")
        print("leaving new results file in script folder")
    else:
        os.mkdir('../results/')
        shutil.move("results_sumstats.txt", "../results/results_sumstats.txt")

def main():
    if len(sys.argv) != 3:
        print("python sumstats.py traits.file /path/to/loci/")
        sys.exit()

    tr = traits(sys.argv[1])
    loci = locus_files(sys.argv[2])

    # store results
    res = {sp:[] for sp in set(tr.values())}
    for index, i in enumerate(loci):
        print("Locus {0}\tPercent complete {1:.2f}%".format(index + 1,
                                                           ((index + 1) / float(len(loci)) * 100)))
        partition = read_file(tr, i)
        results = get_character_matrix(partition)

        # add results for locus to total
        for k, v in results.items():
            res[k].append(v)

    by_stat = put_ss_together(res)
    write_ss(by_stat)
    move_file()

if __name__ == '__main__':
    main()
