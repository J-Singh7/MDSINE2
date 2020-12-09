import mdsine2 as md2
import matplotlib.pyplot as plt
import numpy as np

def _unbias_var_estimate(vals):
    '''Unbiased variance estimate of the values in `vals`
    '''
    vals = np.asarray(vals)
    mean = np.mean(vals)
    a = np.sum((vals - mean)**2)
    return a / (len(vals)-1)

if __name__ == '__main__':
    
    # Alpha
    subjset_healthy = md2.Study.load('../gibson_output/gibson_healthy_agg.pkl')
    subjset_uc = md2.Study.load('../gibson_output/gibson_uc_agg.pkl')
    subjset_inoc = md2.Study.load('../gibson_output/gibson_inoculum_agg.pkl')

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)

    val_healthy = {}
    for subj in subjset_healthy:
        for t in subj.times:
            val = md2.diversity.alpha.normalized_entropy(subj.reads[t])
            if t not in val_healthy:
                val_healthy[t] = []
            val_healthy[t].append(val)
    val_uc = {}
    for subj in subjset_uc:
        for t in subj.times:
            val = md2.diversity.alpha.normalized_entropy(subj.reads[t])
            if t not in val_uc:
                val_uc[t] = []
            val_uc[t].append(val)

    val_inoc_healthy = md2.diversity.alpha.normalized_entropy(
        subjset_inoc['Healthy'].reads[0])
    val_inoc_uc = md2.diversity.alpha.normalized_entropy(
        subjset_inoc['Ulcerative Colitis'].reads[0])

    times = np.sort(list(val_healthy.keys()))
    means_healthy = np.zeros(len(times))
    std_healthy = np.zeros(len(times))
    means_uc = np.zeros(len(times))
    std_uc = np.zeros(len(times))

    for i, t in enumerate(times):
        means_healthy[i] = np.mean(val_healthy[t])
        std_healthy[i] = np.sqrt(_unbias_var_estimate(val_healthy[t]))
        means_uc[i] = np.mean(val_uc[t])
        std_uc[i] = np.sqrt(_unbias_var_estimate(val_uc[t]))

    

        



    





    
