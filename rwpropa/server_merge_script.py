import time
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw
import sys
import matplotlib.pyplot as plt




if __name__ == '__main__':
    file_name = sys.argv[1]
    nr = int(sys.argv[2])
    t = np.array([])
    print('server merge script with parameters: ', file_name, nr)
    for i in range(1,nr):
        df = pd.read_pickle(file_name+str(i)+'.pkl')
        t = np.append(df['d'].values.tolist(), t)
    np.save('merged', t)