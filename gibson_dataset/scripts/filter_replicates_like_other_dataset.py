'''Filter the `replicates` dataset so that it has identical taxas to 
another dataset. We cannot pass the `replicates` dataset into a consistency
filter because there is no such thing as a consecutive timepoint in that 
dataset
'''
import mdsine2 as md2
import argparse
import logging

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--replicate-dataset', '-r', type=str, dest='replicate',
        help='Location of the replicate dataset')
    parser.add_argument('--like-other', '-l', type=str, dest='other',
        help='Location of the other dataset to filter like')
    parser.add_argument('--output-basepath', '-o', type=str, dest='path',
        help='Location to save the output')
    args = parser.parse_args()
    md2.config.LoggingConfig(level=logging.INFO)

    replicates = md2.Study.load(args.replicate)
    other = md2.Study.load(args.other)

    to_delete = []
    for taxa in replicates.taxas:
        if taxa.name not in other.taxas:
            to_delete.append(taxa.name)
    replicates.pop_taxas(to_delete)
    replicates.save(args.path)
    
    