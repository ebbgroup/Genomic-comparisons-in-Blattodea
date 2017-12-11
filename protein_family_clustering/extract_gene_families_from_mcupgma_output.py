#!/usr/bin/python3.4


import argparse
import sys
import gzip
import numpy as np


sys.setrecursionlimit(2000)



# gene ID prefix mapped to short species name
ID2SP = {
    'SMAR': 'sma',
    'EDAN': 'eda',
    'RPR': 'rpr',
    'HSAL': 'hsa',
    'LH': 'lhu',
    'CFLO': 'cfl',
    'PB': 'pba',
    'SINV': 'sin',
    'AECH': 'aec',
    'ACEP': 'ace',
    'GB': 'ame',
    'Nasvi': 'nvi',
    'FBpp': 'dme',
    'AAE': 'aae',
    'XM_': 'tca',
    'NM_': 'tca',
    'Znev_': 'zne',
    'Csec_': 'cse',
    'Mnat_': 'mna',
    'Bger_': 'bge',
    'LOCMI': 'lmi',
    'rna': 'pca',

    'Acep_': 'ace',
    'l(2)efl-PA': 'ace',  # why...
}

SPECIES = set(ID2SP.values())
OK_SPECIES = sorted(SPECIES - {'sma', 'eda'})  # exclude outgroups from clusters


def parse_upgma_tree(infile):
    print('parsing MC-UPGMA clustering (step 1)')
    tree_d = {}
    with open(infile, 'r') as tree:
        for row in tree:
            values = row.rstrip().split()
            id1, id2, evalue, cluster_id = \
                int(values[0]), int(values[1]), float(values[2]), int(values[3])
            assert cluster_id not in tree_d, 'duplicate cluster definition'
            tree_d[cluster_id] = (id1, id2)
    return tree_d


def parse_numeric2id_map(infile):
    print('parsing ID map')
    map_d = {}
    with gzip.open(infile, 'r') as tsvmap:
        for row in tsvmap:
            values = row.decode().rstrip().split()
            num_id = int(values[0])
            str_id = values[1]
            assert num_id not in map_d, 'duplicate id in map'
            for id_pfx, species in ID2SP.items():
                if str_id.startswith(id_pfx):
                    species_and_id = '{}|{}'.format(species, str_id)
                    break
            else:
                raise Exception('Couldn\'t determine species for {}'.format(str_id))
            map_d[num_id] = species_and_id
    return map_d


def resolve_cluster(highest_id, num2txt, cluster2ids, cluster2genes):
    if highest_id in cluster2genes:
        genes_in_cluster = cluster2genes[highest_id]
    else:
        genes_in_cluster = set()
        def _deeper(cluster_id):
            #print('\nthe ID', cluster_id, 'defines:')
            if cluster_id in num2txt:
                #print(num2txt[cluster_id])
                genes_in_cluster.add(num2txt[cluster_id])
            else:
                subclusters = cluster2ids[cluster_id]
                #print('{0} and {1}'.format(*subclusters))
                for sub_id in cluster2ids[cluster_id]:
                    _deeper(sub_id)
        _deeper(highest_id)
        cluster2genes[highest_id] = genes_in_cluster
    return genes_in_cluster




def main(mcupgma_tree, numeric2text, protolevel, max_diff, max_mean_size):

    num2txt = parse_numeric2id_map(numeric2text)
    cluster2ids = parse_upgma_tree(mcupgma_tree)
    cluster2genes = {}

    print('parsing MC-UPGMA clustering (step 2)')
    with open(mcupgma_tree, 'r') as tree:
        for line in reversed(list(tree)):
            values = line.rstrip().split()
            clusters_protolevel = float(values[2])
            clust_id = int(values[3])
            if clusters_protolevel <= protolevel:
                genes_in_cluster = resolve_cluster(
                    clust_id, num2txt, cluster2ids, cluster2genes
                )

    genes_already_printed = set()

    print('writing clusters to files if their size is ok...')
    clusters_sorted_by_size = sorted(
        cluster2genes, key=lambda k: len(cluster2genes[k]), reverse=True)

    # get some descriptive output file names:
    params = 'p{}'.format(round(protolevel))
    if max_diff != float("inf"):
        params += '_d{}'.format(round(max_diff))
    if max_mean_size != float("inf"):
        params += '_s{}'.format(round(max_mean_size))
    orthofile_n = 'mcupgma_ortho_{}.txt'.format(params)
    famfile_n = 'mcupgma_famsizes_{}.tsv'.format(params)

    with open(orthofile_n, 'w') as orthofile, open(famfile_n, 'w') as famfile:

        fam_header = ['Description', 'ID'] + OK_SPECIES
        famfile.write('\t'.join(fam_header) + '\n')
        big_fams = 0

        for clust_id in clusters_sorted_by_size:
            gene_ids = cluster2genes[clust_id]

            speciesnames = [g.split('|')[0] for g in gene_ids]
            famsizes = [speciesnames.count(s) for s in OK_SPECIES]
            mean_famsize = np.mean(famsizes)
            sizediff = max(famsizes) - min(famsizes)

            if mean_famsize > max_mean_size:
                # cluster is too big! split it into smaller subclusters
                print('Splitting {} because mean famsize {} is too big'.format(clust_id, mean_famsize))
                continue

            if sizediff > max_diff:
                # cluster is too big! split it into smaller subclusters
                print('Splitting {} because size difference {} is too big'.format(clust_id, sizediff))
                continue

            if gene_ids & genes_already_printed:
                # a higher-level cluster was already printed
                continue

            orthofile.write('C_{}: {}\n'.format(clust_id, ' '.join(gene_ids)))

            if sum(famsizes) >= 100:
                big_fams += 1

            famfile.write('n/a\tC_{}\t{}\n'.format(clust_id, '\t'.join([str(f) for f in famsizes])))

            genes_already_printed |= gene_ids

    print('Number of families with over 100 genes:', str(big_fams))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parses MC-UPGMA result tree and converts it to a usable format')
    parser.add_argument('mcupgma_tree', help='path to the MC-UPGMA output tree file')
    parser.add_argument('numeric2text', help='path to the tsv that maps numerical gene IDs '
                        'to regular gene IDs. generated by 1_get_edges_input_for_mcupgma.sh')
    parser.add_argument('-p', '--protolevel', default=100.0, type=float,
                        help='smaller ProtoLevel cutoff = smaller clusters (0-100; default: 100)')
    parser.add_argument('-d', '--max_diff', default=float('inf'), type=float,
                        help='split clusters into subclusters if the family size '
                        'differs too much between the species with the largest family '
                        'and the species with the smallest family. This parameter '
                        'specified the maximum tolerated difference (default: never split)')
    parser.add_argument('-s', '--max_mean_size', default=float('inf'), type=float,
                        help='split clusters into subclusters if the mean family size exceeds '
                        'this threshold (default: never split)')
    args = parser.parse_args()
    main(**vars(args))
