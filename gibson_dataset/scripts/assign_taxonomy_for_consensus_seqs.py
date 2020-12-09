'''Reassign the taxonomy based on the consensus sequences that were 
computed with the script `preprocess.py` with RDP

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6
'''
import pandas as pd
import mdsine2 as md2
import os
import argparse
import logging

def parse_rdp(fname, confidence_threshold):
    '''Parse the taxonomic assignment document from RDP with a confidence
    threshold `confidence_threshold`

    Parameters
    ----------
    fname : str
        This is the name of the taxonomic assignment document from RDP
    confidence_threshold : float
        This is the minimum confidence needed for us to include it in the
        classification

    Returns
    -------
    dict( key1 -> dict ( key2 -> value ) )
        key1 : str
            OTU name
        key2 : str
            taxonomic level
        value : str
            taxonomic name
    '''
    f = open(fname, 'r')
    txt = f.read()
    f.close()

    columns = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    data = []
    index = []

    for i, line in enumerate(txt.split('\n')):
        if 'OTU' != line[:3]:
            continue

        splitting = line.split('%;')
        otuname = splitting[0].split(';')[0]
        splitting = splitting[1:]

        index.append(otuname)

        temp = []
        taxaed = []

        for tax_idx in range(len(splitting)):
            tax_key = columns[tax_idx]
            tax, confidence = splitting[tax_idx].replace('%', '').split(';')
            if float(confidence) > confidence_threshold:
                temp.append(tax)
                taxaed.append(tax_key)
            else:
                break

        for tax_key in columns:
            if tax_key not in taxaed:
                temp.append(md2.pylab.base.DEFAULT_TAXA_NAME)
        data.append(temp)

    df = pd.DataFrame(data, columns=columns, index=index)
    return df

def parse_silva_NOT_USED(fname):
    '''Parse the taxonomic assignment table from Silva

    Parameters
    ----------
    fname : str
        This is the name of the taxonomic assignment table from Silva

    Returns
    -------
    dict( key1 -> dict ( key2 -> value ) )
        key1 : str
            OTU name
        key2 : str
            taxonomic level
        value : str
            taxonomic name
    '''
    tbl = pd.read_csv(fname, sep='\t', index_col=0)
    for otuname in tbl.index:
        if 'OTU' not in otuname:
            continue

        d_silva[otuname] = {}
        
        taxas = tbl['Taxon'][otuname].split('; ')
        for tax in taxas:
            if 'd__' in tax:
                key = 'tax_kingdom'
                tax = tax.replace('d__', '')
            elif 'p__' in tax:
                key = 'tax_phylum'
                tax = tax.replace('p__', '')
            elif 'c__' in tax:
                key = 'tax_class'
                tax = tax.replace('c__', '')
            elif 'o__' in tax:
                key = 'tax_order'
                tax = tax.replace('o__', '')
            elif 'f__' in tax:
                key = 'tax_family'
                tax = tax.replace('f__', '')
            elif 'g__' in tax:
                key = 'tax_genus'
                tax = tax.replace('g__', '')
            elif 's__' in tax:
                key = 'tax_species'
                tax = tax.replace('s__', '')

                # replace the {genus}_ prefix of the species
                tax = tax.replace(d_silva[otuname]['tax_genus'] + '_', '')

            if 'uncultured' in tax:
                tax = None
                # Break here - everything under an uncultured does not make sense
                break
            d_silva[otuname][key] = tax
    return d_silva

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--rdp-table', '-r', type=str, dest='rdp_table',
        help='Location of RDP file')
    parser.add_argument('--confidence-threshold', '-c', type=float, dest='confidence_threshold',
        help='This is the minimum confidence required for us to use the classification')
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you want to save the parsed dataset.')
    args = parser.parse_args()

    md2.config.LoggingConfig(level=logging.INFO)

    logging.info('Parsing RDP')
    df = parse_rdp(fname=args.rdp_table, confidence_threshold=args.confidence_threshold)

    for dset in ['healthy', 'uc', 'replicates', 'inoculum']:
        logging.info('Replacing {}'.format(dset))
        study_fname = os.path.join(args.basepath, 'gibson_{dset}_agg.pkl'.format(dset=dset))
        study = md2.Study.load(study_fname)

        study.taxas.generate_consensus_taxonomies(df)
        study_fname = os.path.join(args.basepath, 'gibson_{dset}_agg_taxa.pkl'.format(dset=dset))
        study.save(study_fname)



    
    



    

    
