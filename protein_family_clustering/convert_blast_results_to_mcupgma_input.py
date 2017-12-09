#!/usr/bin/python3.4


import argparse
import gzip
import io


gene2id = {}
edge2eval = {}


def main(blastfile, blast_evalue_cutoff):
    id_cnt = 0
    with open(blastfile, 'r') as infile, gzip.open('edges.gz', 'wb') as out_raw, gzip.open('numeric2text_protID.tsv.gz', 'wb') as idmap_raw:
        with io.TextIOWrapper(out_raw, encoding='utf-8') as outfile, io.TextIOWrapper(idmap_raw, encoding='utf-8') as idmap:
            for i, row in enumerate(infile):
                if i % 10000 == 0:
                    print('Parsing BLAST file row {}'.format(i), end='\r')
                values = row.split()
                queryId = values[0].split('|')[1]
                subjectId = values[1].split('|')[1]

                # Self edges are not allowed, i.e cluster_id1 == cluster_id2 is illegal.
                if queryId == subjectId:
                    #print('Skipped self edge', queryId, subjectId)
                    continue

                blast_evalue_str = values[10]
                blast_evalue = float(blast_evalue_str)

                # filter hits with poor e-values:
                if blast_evalue > blast_evalue_cutoff:
                    continue

                # assign a unique numerical ID to each gene ID (if it doesnt have one already)
                for prot_id in (queryId, subjectId):
                    if prot_id not in gene2id:
                        id_cnt += 1
                        gene2id[prot_id] = id_cnt
                        idmap.write('{}\t{}\n'.format(id_cnt, prot_id))

                # It is illegal to include both edges i<->j and j<->i and it will lead to unexpected clustering failures. 
                edge_id = tuple(sorted([gene2id[queryId], gene2id[subjectId]]))
                # only take the lowest hit:
                if edge_id not in edge2eval or blast_evalue < edge2eval[edge_id]:
                    edge2eval[edge_id] = blast_evalue


            print('\nWriting output...\n')
            i = 0
            for edge_tuple, best_evalue in edge2eval.items():
                i += 1
                if i % 1000 == 0:
                    print('Writing graph edge pair {}'.format(i), end='\r')
                edge1, edge2 = edge_tuple
                assert edge1 < edge2
                assert best_evalue < blast_evalue_cutoff
                assert best_evalue >= 0
                outline = '{}\t{}\t{}\n'.format(edge1, edge2, best_evalue)
                outfile.write(outline)
    
    print('Maximum ID is', id_cnt)
    print('Wrote output file to edges.tsv. Gzip this file and call it "x.edges.gz" '
          'so it can be used as an input for MC-UPGMA.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='converts blast all vs. all result to edges.gz (mc-upgma input format) '
            'NOTE: output is not actually gzipped, this has to be done manually')
    parser.add_argument('blastfile', help='a BLAST all vs all output table')
    parser.add_argument('-c', '--blast_evalue_cutoff', help='default: 10.0', type=float, default=10.0)
    args = parser.parse_args()
    main(**vars(args))
