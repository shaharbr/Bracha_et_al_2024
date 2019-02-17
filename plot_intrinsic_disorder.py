import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy import stats

# define full list of constructs tested
proteins =['GRA16-HA-MECP2', 'GRA16-HA-TFEB', 'GRA16-HA-SMN1', 'GRA16-HA-ASPA', 'GRA16-HA-ASPAopt', 'GRA16-HA-GALC', 'GRA16-HA-GALCopt','GRA16-HA-GALC-TAT']
# add mean nuclear intensity (nuclear localization) of each protein from separate analysis,
# in the same order as the protein names
nucintens =[80.96, 76.94, 51.49, 4.958, 3.539, 6.551, 2.372, 7.073]
# add them together in a dictionary
proteins = dict((proteins[i], {'nucintens': nucintens[i]}) for i in range(len(proteins)))
# define list of proteins sequences (opt and non-opt encode for the same protein,
# so to avoid repeats opt are not included)
proteinlist = ['GRA16-HA-MECP2', 'GRA16-HA-TFEB', 'GRA16-HA-SMN1', 'GRA16-HA-ASPA', 'GRA16-HA-GALC', 'GRA16-HA-GALC-TAT']

# define figure
f, axes = plt.subplots(3, 2, figsize=(14, 9))

# define empty list which will be iteratively filled with the calculated disorder score for each protein
disorder_df = []
# define disorder threshold
threshold = 0.5

# for each unique protein sequence:
for protein in proteins:
    # import table with calculated IPUred2 and ANCHOR2 scores, calculated by https://iupred2a.elte.hu/
    df = pd.read_csv(protein+ '.csv', header=3)
    # add column with the protein name
    df['protein'] = protein
    # calculate the average IUPRED score and ANCHOR score, over the span of the therapeutic protein
    # (excluding GRA16-HA = from amino acid 519 until the end of the protein)
    IUPRED_AVG = np.average(df.loc[df['# POS']>518]['IUPRED SCORE']) - threshold
    ANCHOR_AVG = np.average(df.loc[df['# POS']>518]['ANCHOR SCORE']) - threshold
    # calculate the average of the two scores into one combined average score = the disorder score
    # (in a new column)
    proteins[protein]['disorder_score'] = (IUPRED_AVG + ANCHOR_AVG) / 2
    disorder_df.append(df)

# concatenate all the protein dataframes which contain the calculated disorder score
df = pd.concat(disorder_df)

# plot the protein disorder profiles of each protein
for protein in proteinlist:
    # the dataframe for that protein
    df_prot = df.loc[df['protein'] == protein]
    # to plot the threshold line, generate list with constatt y = threshold and the length of the protein
    threshold_line = [threshold] * len(df_prot['# POS'])
    # get index of this protein in the list
    i = proteinlist.index(protein)
    # in order to organize the subplots in two columns with 3 subplots each, use this formula
    loc = (i%3,int(i/3))
    # background semi-transparent grid
    axes[loc].grid(alpha=0.5, zorder=1, linestyle=':')
    # plot the threshold line
    axes[loc].plot(df_prot['# POS'], threshold_line, color='black', zorder=2, linewidth=1)
    # plot the IUPRED score of that protein
    axes[loc].plot(df_prot['# POS'], df_prot['IUPRED SCORE'], color='blue', label='IUPred2', zorder=2, linewidth=1)
    # plot the ANCHOR score of that protein
    axes[loc].plot(df_prot['# POS'], df_prot['ANCHOR SCORE'], color='red', label='ANCHOR2', zorder=2, linewidth=1)
    # add text box with the averaged disorder score calculated for that protein
    axes[loc].text(0.99, 1.4, 'Disorder Score: ' + str(round(proteins[protein]['disorder_score'], 2)), fontsize=15,
                   horizontalalignment='right', verticalalignment='top',
                   transform=axes[loc].transAxes, bbox=dict(facecolor='white', alpha=0.5))
    # draw semi-transparent light gray block on top of the common GRA16-HA part of the protein,
    # and also a dashed black line where the fused therapeutic protein starts
    axes[loc].axvspan(0, 518, alpha=0.7, color='#f2f2f2', zorder=2.5)
    axes[loc].axvline(x=518, color='black', linestyle='--', zorder=2.7)
    # set title as the name of the protein
    axes[loc].set_title(protein, fontsize=19, loc='left', y=1.3)
    # formatting
    axes[loc].set_ylabel('Disorder score (AU)', fontsize=16)
    axes[loc].set_xlabel('Position (AA)', fontsize=16)
    axes[loc].legend(bbox_to_anchor=(0.01, 0.32), loc=2, fontsize=10, borderaxespad=0., fancybox=True, framealpha=0.5)
    axes[loc].set_ylim(0, 1)
    axes[loc].yaxis.set_major_locator(ticker.MultipleLocator(0.25))
    axes[loc].set_xlim(left=0)
    axes[loc].tick_params(axis='both', which='major', labelsize=14)
f.subplots_adjust(hspace = 0.9, wspace=0.2, left=0.05, right=0.95)

f.savefig("disorder_profiles.png", bbox_inches="tight")
f.savefig("disorder_profiles.svg", bbox_inches="tight")



# generate linear regression model for the correlation between the intrinsic disorder of the protein and
# the nuclear localization (intensity) of the protein fused to GRA16
x = [proteins[prot]['disorder_score'] for prot in proteins]
y = [proteins[prot]['nucintens'] for prot in proteins]
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

# generate figure showing the correlation between intrinsic disorder score and nuclear localization
f1, ax1 = plt.subplots(1, 1, figsize=(20, 10))
# plot the calculated regression line
ax1.plot(x, intercept + slope*np.array(x), 'gray', label='fitted line')
# scatter plot all the data points
ax1.scatter(x, y, zorder=10, color='#009999')
# add a text box with the calculated R^2 and P value
ax1.text(0.03, 0.97, '$R^2$: '+str(round(r_value, 3)) + '\nP value: '+str(round(p_value, 5)), fontsize=15, horizontalalignment = 'left', verticalalignment = 'top',
                   transform = ax1.transAxes, bbox=dict(facecolor='white', alpha=0.5))
# formatting
ax1.set_title('Disorder Score vs. Host Nuclear Localization', fontsize=25, y=1.1)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_xlabel('Disorder score (AU)', fontsize=18)
ax1.set_ylabel('Host Nuclear Fluorescence \nIntensity (AFU)', fontsize=18)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
f1.subplots_adjust(bottom=0.2, top=0.8, left=0.3)

f1.savefig("disorder_vs_nuclear_localization.png", bbox_inches="tight")
f1.savefig("disorder_vs_nuclear_localization.svg", bbox_inches="tight")


plt.show()
