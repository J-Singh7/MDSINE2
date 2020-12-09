'''Write the contents of the Study object into 3 files that are used to 
pass data into the clv paper: metadata.txt, counts.txt, biomass.txt

'''
import mdsine2 as md2
import pandas as pd
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--dataset', '-d', type=str, dest='dataset',
        help='This is the Gibson dataset we want to do cross validation on')
    parser.add_argument('--basepath', '-o', type=str, dest='basepath',
        help='This is the folder you want to save the documents')
    args = parser.parse_args()
    md2.config.LoggingConfig()

    os.makedirs(args.basepath, exist_ok=True)
    study = md2.Study.load(args.dataset)

    # metadata.txt
    columns = ['sampleID', 'isIncluded', 'subjectID', 'measurementid', 'perturbid']
    data = []
    sampleid = 0
    for subj in study:
        for t in subj.times:
            sampleid += 1
            ppp = 0
            for pidx, perturbation in enumerate(study.perturbations):
                if perturbation.isactive(time=t, subj=subj.name):
                    ppp = pidx+1
                    break

            temp = [sampleid, 1, int(subj.name), t, ppp]
            data.append(temp)

    df = pd.DataFrame(data, columns=columns)
    df.to_csv(os.path.join(args.basepath, 'metadata.txt'), 
        sep='\t', index=False, header=True)

    # counts.txt
    columns = ['#OTU ID'] + [a+1 for a in range(sampleid)]
    data = []
    for taxa in study.taxas:
        temp = [taxa.name]
        sampleid = 1
        for subj in study:
            for t in subj.times:
                temp.append(subj.reads[t][taxa.idx])
        data.append(temp)

    
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(os.path.join(args.basepath, 'counts.txt'), 
        sep='\t', index=False, header=True)

    # biomass.txt
    columns = ['mass1', 'mass2', 'mass3']
    data = []
    for subj in study:
        for t in subj.times:
            data.append(subj.qpcr[t].data)
    
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(os.path.join(args.basepath, 'biomass.txt'), 
        sep='\t', index=False, header=True)

    