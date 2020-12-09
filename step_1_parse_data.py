'''Parse the input files

mdsine2 version: 4.0.6
Author : David Kaplan
Date: 11/30/20

Input tables
------------
Metadata
    Holds the metadata for each sample. Columns:
    sampleID : str
        Name of the sample
    subject : str
        Name of the subject this sample belongs to
    time : float
        Timepoint of the sample

Perturbations
    These are the perturbations for each subject. Columns
    name : str
        Name of the perturbation
    start, end : float
        Start and end of the perturbation
    subject : str
        This is the subject this perturbation corresponds to. Subject must
        be contained in the metadata file as well.

qPCR
    These are the qPCR measurements for each sample. Columns:
    sampleID : str
        Name of the sampe
    measurement1, ... : float
        Rest of the columns are the replicate measurements

Counts
    These are the counts for each taxa. Columns:
    name : str
        Name of the taxa
    `sampleID`s
        Each sampleID has its own count

Taxonomy
    This is the taxonomy name for each taxa in counts. Columns:
    name : str
        Name of the taxa. This corresponds to the name in `counts`
    sequence : str
        Sequence associated with the taxa
    kingdom : str
        Kingdom taxonomic classification
    phylum : str
        Phylum taxonomic classification
    class : str
        Class taxonomic classification
    order : str
        Order taxonomic classification
    family : str
        Family taxonomic classification
    genus : str
        Genus taxonomic classification
    species : str
        Species taxonomic classification
'''
import argparse
import mdsine2 as md2
import logging

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--name', '-n', type=str, dest='name',
        help='Name of the dataset')
    parser.add_argument('--taxonomy', '-t', type=str, dest='taxonomy',
        help='This is the table showing the sequences and the taxonomy for each' \
            ' ASV or OTU')
    parser.add_argument('--metadata', '-m', type=str, dest='metadata',
        help='This is the metadata table')
    parser.add_argument('--reads', '-r', type=str, dest='reads',
        help='This is the reads table', default=None)
    parser.add_argument('--qpcr', '-q', type=str, dest='qpcr',
        help='This is the qPCR table', default=None)
    parser.add_argument('--perturbations', '-p', type=str, dest='perturbations',
        help='This is the perturbation table', default=None)
    parser.add_argument('--sep', '-s', type=str, dest='sep',
        help='This is the separator for the tables', default='\t')
    parser.add_argument('--outfile', '-o', type=str, dest='outfile',
        help='This is where you want to save the parsed dataset')
    args = parser.parse_args()

    md2.config.LoggingConfig(level=logging.INFO)
    study = md2.dataset.parse(name=args.name, metadata=args.metadata, taxonomy=args.taxonomy,
        reads=args.reads, qpcr=args.qpcr, perturbations=args.perturbations, sep=args.sep)
    study.save(args.outfile)