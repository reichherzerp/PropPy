import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file_name = sys.argv[1]
    nr_simulations = int(sys.argv[2])
    t = np.array([])
    print('server merge script with parameters: ', file_name, nr_simulations)
    for i in range(1, nr_simulations):
        df = pd.read_pickle(file_name+str(i)+'.pkl')
        t = np.append(df['d'].values.tolist(), t)
    np.save(file_name+'_merged', t)
    print('saved merged file')

    # plot lin-lin-hist
    plt.figure(figsize=(5,3))
    bins = 30
    plt.hist(t, bins=bins, alpha=0.5)
    plt.title('total observed particles = {:.0e}'.format(len(t)))
    plt.xlabel('D [m]')
    plt.ylabel('# particles')
    plt.tight_layout()
    plt.savefig(file_name+'_hist_lin.png')

    # plot log-log-hist
    plt.figure(figsize=(5,3))
    d = t
    hist, bins = np.histogram(d, bins=bins)
    logbins = np.logspace(np.log10(min(d)),np.log10(max(d)),len(bins))
    plt.hist(d, bins=logbins, alpha=0.5)
    plt.title('total observed particles = {:.0e}'.format(len(t)))
    plt.xlabel('D [m]')
    plt.ylabel('# particles')
    plt.loglog()
    plt.tight_layout()
    plt.savefig(file_name+'_hist_log.png')
    plt.savefig(file_name+'_hist_log.pdf')

    print('created histogram plots')
    
