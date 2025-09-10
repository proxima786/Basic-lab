from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

cigale_results=Table.read('results.txt', format='ascii')
cigale_results['ssfr'] = cigale_results['bayes.sfh.sfr10Myrs']/cigale_results['bayes.stellar.m_star']

cigale_results['log_ssfr'] = np.log10(cigale_results['ssfr'])

mask_quiescent = (cigale_results['log_ssfr'] < -10.5)
mask_starforming = (cigale_results['log_ssfr'] > -10.5)

quiescent_table = cigale_results[mask_quiescent]
starforming_table = cigale_results[mask_starforming]

quiescent_table.write('results_quiescent.txt', format='ascii', overwrite=True)
starforming_table.write('results_starforming.txt', format='ascii', overwrite=True)
