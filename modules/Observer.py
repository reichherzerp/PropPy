import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

class Plotter():
    def __init__(self):
        print('init plotter')
        
    def distribution(self, data):
        print('start plotting')
        plt.figure()
        # Fit a normal distribution to the data:
        mu, std = norm.fit(data)

        # Plot the histogram.
        plt.hist(data, bins=100, density=True, alpha=0.6, color='g')

        # Plot the PDF.
        xmin, xmax = plt.xlim()
        x_values = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x_values, mu, std)
        plt.plot(x_values, p, 'k', linewidth=2)
        title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
        plt.title(title)

        plt.show()