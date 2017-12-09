#!/usr/bin/env python3

'''
This script is used to cluster branches of the phylogeny with k-means.
Per-branch lambda estimations are used as input (see "per_branch" directory).

For convenience, this script also creates a new CAFE script (and a script
for submission to a computational cluster with SGE qsub) for the branch
combinations determined with kmeans. 

usage:
grep -A 1 "Lambda Search Result" ../per_branch/*.sh.out | grep "Score" | ../kmeans.py
'''


import sys
import os
import re
import argparse
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import scale


cafe_file = """#!shell
version
date

tree ((((nvi:162,(ame:97,(pca:76,((((((aec:9,ace:9):9,sin:18):9,pba:27):10,cfl:37):9,lhu:46):9,hsa:55):21):21):65):182,(tca:327,(aae:158,dme:158):169):17):29,rpr:373):14,((((mna:121,cse:121):15,zne:136):37,bge:173):75,lmi:248):139)

load -filter -r 10000 -i ../mcupgma_famsizes_p80_d100_s35_noIRs_noEdan_plusManuallyAnnotatedFamilies_noTEs_withPolistes.tsv -t 6 -l k{k:02d}_r{i}_CAFE.out -p 0.05

lambda -s -t {lambdatree}

report k{k:02d}_r{i}

date
"""


def stdin2lambda_df():
    best_rows = {}
    for raw_line in sys.stdin:
        values = raw_line.strip().split()
        lambda1, lambda2 = (float(l) for l in values[2].split(','))
        branch = os.path.basename(values[0]).split('_')[0]
        likelihood = float(values[-1])
        print(branch, lambda1, lambda2, likelihood)
        if branch not in best_rows:
            best_rows[branch] = [branch, lambda1, lambda2, likelihood]
        else:
            if likelihood < best_rows[branch][-1]:
                best_rows[branch] = [branch, lambda1, lambda2, likelihood]
    df_raw = sorted(list(best_rows.values()), key=lambda x: int(x[0][1:]))
    return pd.DataFrame(df_raw, columns=('branch', 'foreground-lambda', 'background-lambda', 'likelihood'))


cafe_qsubber = """#!/bin/bash
#$ -N {cafe_sh}
#$ -o {wd}/{cafe_sh}.out
#$ -e {wd}/{cafe_sh}.err
#$ -pe smp 6
#$ -V
#$ -S /bin/bash
#$ -wd {wd}/


echo "CAFE job {cafe_sh} started on"; hostname

./cafe.linux.x86_64 {cafe_sh}
"""


def sort_groups(lambdatree):
    i, j = 0, 0
    old2new = {}
    out_d = {}
    out_chars = []
    for char in re.split(r'[\(\)\,]', lambdatree):
        if char.strip():
            i += 1
            branch = 'b' + str(i)
            old2new[branch] = char
            if char not in old2new:
                j += 1
                old2new[char] = str(j)
            out_d[branch] = old2new[char]
    return out_d


def main(k, r):
    k1, k2 = k.split(':')
    kstart, kend = int(k1), int(k2) + 1
    i1, i2 = r.split(':')
    istart, iend = int(i1), int(i2) + 1

    df = stdin2lambda_df()  # pd.read_table('lambdas.tsv', index_col=False)#, header=None, names=headers)

    print(df)

    df['branch'] = df['branch'].apply(lambda x: 'b' + ''.join([l for l in x if l.isdigit()]))
    df = df.sort_values('branch')

    relevant_cols = scale(df[['foreground-lambda', 'background-lambda']])

    #tree = '(((({b1},({b2},(((((({b3},{b4}){b5},{b6}){b7},{b8}){b9},{b10}){b11},{b12}){b13},{b14}){b15}){b16}){b17},({b18},({b19},{b20}){b21}){b22}){b23},{b24}){b25},(((({b26},{b27}){b28},{b29}){b30},{b31}){b32},{b33}){b34})'
    tree = '(((({b1},({b2},({b3},(((((({b4},{b5}){b6},{b7}){b8},{b9}){b10},{b11}){b12},{b13}){b14},{b15}){b16}){b17}){b18}){b19},({b20},({b21},{b22}){b23}){b24}){b25},{b26}){b27},(((({b28},{b29}){b30},{b31}){b32},{b33}){b34},{b35}){b36})'

    for k in range(kstart, kend):
        X = np.matrix(relevant_cols)
        kmeans = KMeans(n_clusters=k).fit(X)
        df['cluster'] = [c + 1 for c in kmeans.labels_]
        print('\nk-means clusters for k =', k)
        #print(df)
        branch2lambda = {}
        for row in df.itertuples():
            branch2lambda[row.branch] = row.cluster
        lambdatree_unsorted = tree.format(**branch2lambda)
        lambda_d = sort_groups(lambdatree_unsorted)
        lambdatree = tree.format(**lambda_d)
        print('lambda -s -t ' + lambdatree)

        for i in range(istart, iend):
            cwd = os.getcwd()
            shname = 'k{k:02d}_r{i}_CAFE.sh'.format(k=k, i=i)
            qsubname = 'k{k:02d}_r{i}_CAFE_qsub.sh'.format(k=k, i=i)
            with open(shname, 'w') as shfile:
                shfile.write(cafe_file.format(
                    k=k,
                    i=i,
                    lambdatree=lambdatree,
                ))
            with open(qsubname, 'w') as qfile:
                qfile.write(cafe_qsubber.format(
                    wd=cwd,
                    cafe_sh=shname,
                    lambdatree=lambdatree,
                ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='CAFE k-means qsub script generator! usage: grep -A 1 "Lambda Search Result" ../per_branch/b*.sh.out | grep "Score" | ../kmeans.py -k 2:20 -i 1:3')
    parser.add_argument('-k', help='the range of k-values to try out for k-means '
                        'clustering, example: 2:10', default='2:20')
    parser.add_argument('-r', help='the number of times each script should be run, '
                        'example: 1:3 to generate three runs with IDs 1,2,3', default='1:5')
    args = parser.parse_args()
    main(**vars(args))
