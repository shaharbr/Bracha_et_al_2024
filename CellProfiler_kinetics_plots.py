import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker



# load dataframe after cleaning and organization by "HFF_cellprofiler_output"
df_HFF = pd.read_csv('statdf_HFF_perimage_all_25.11.csv', usecols=range(0,10))
df_HFF_inf = pd.read_csv('statdf_HFF_perimage_inf_25.11.csv', usecols=range(0,10))

df_luhmes = pd.read_csv('statdf_luhmes_perimage_all_25.11.csv', usecols=range(0,10))
df_luhmes_inf = pd.read_csv('statdf_luhmes_perimage_inf_25.11.csv', usecols=range(0,10))

# define color palette (should be colorblind-friendly)
greens = ['#63B7B7', '#007C80', '#014D4B']
purples = ['#D3A7FF', '#9148E5', '#59228F']
oranges = ['#FFA500', '#D38B07', '#A7710E']
# arrange colors in list according to how they will be called for in the loop
colors = [oranges, greens, purples]

# define function for plotting a requested variable over all strains and conditions.
# received the dataframe, the name of the column that should be plotted, the name of the column that should be
# used as errorbars, whether yaxis should be shared, the title of the figure, the labled of the xaxis,
# whether there should be an ylim, and the requested xticks interval.
def timeplot(df, ycolumn, stdcolumn, sharey, title, ylabel, ylim, xticks, colorlist=colors):
    """"
    Plotting a requested variable over all strains and conditions.
    Receives:
    df: dataframe
    ycolumn: the name of the column that should be plotted
    stdcolumn: the name of the column that should be used as errorbars
    sharey: whether yaxis should be shared
    title: the title of the figure
    ylabel: the label of the y-axis
    ylim: The ylim (default None)
    xticks: the requested xticks interval
    colorlist: using the default color list "colors" defined above.
    """
    # sort the dataframes
    mois = sorted(df['moi'].unique())
    strains = sorted(df['strain'].unique())
    # get number of strains, which will dictate the number of subplots
    numstrains = len(strains)
    # make figure with number of subplots = number of strains
    if numstrains == 1:
        f, axes = plt.subplots(1, numstrains, sharey=sharey, sharex=True, figsize=(7, 6))
    else:
        f, axes = plt.subplots(1, numstrains, sharey=sharey, sharex=True, figsize=(9, 6))
    plt.rcParams['legend.handlelength'] = 0
    # in order to get around matplotlib's deafult behaviour of turning a list of "one subplot" 
    # into a non-indexable variable, I put it inside a list again so that 
    # I could call it by an index in the loop below
    if numstrains == 1 :
        axes = [axes]
    for i in range(numstrains):
        for m in range(len(mois)):
            df_strain_moi = df.loc[(df['strain'] == strains[i]) & (df['moi'] == mois[m])]
            # plot a line plot for each moi, for each strain subplot. plot the variable defined in the input.
            df_strain_moi.plot(x='timepoint', y=ycolumn, yerr=stdcolumn, capsize=3, capthick=1, linewidth=2,
                               marker='o', ax=axes[i], label=str(mois[m]), color=colorlist[i][m])
            # xaxis always time. set only for the middle graph
            if i == int(numstrains / 2):
                axes[i].set_xlabel('Hours post infection', fontsize=20)
            else:
                axes[i].set_xlabel('')
        # ylable title by input. set only for the left-most graph.
        if i == 0:
            axes[i].set_ylabel(ylabel, fontsize=20)
        # title of subplot is the strain
        axes[i].set_title('GRA16-' + strains[i], fontsize=20, y=1.1)
        axes[i].tick_params(axis="both", which="both", bottom="off", top="off",
                            labelbottom="on", left="off", right="off", labelleft="on", labelsize=20)
        # xticks interval by input
        axes[i].xaxis.set_major_locator(ticker.MultipleLocator(xticks))
        # remove the errorbars from the legends
        handles, labels = axes[i].get_legend_handles_labels()
        handles.reverse()
        labels.reverse()
        handles = [h[0] for h in handles]
        # make legends of moi label, with the title "MOI"
        legend = axes[i].legend(handles, labels, title="MOI", loc='upper left', fancybox=True, framealpha=0.5,
                                fontsize=16)
        plt.setp(legend.get_title(), fontsize='16')
        # remove black box around subplots
        axes[i].spines["top"].set_visible(False)
        axes[i].spines["bottom"].set_visible(False)
        axes[i].spines["right"].set_visible(False)
        axes[i].spines["left"].set_visible(False)
        # add semi-transparent yaxis grid in the background
        axes[i].yaxis.grid(True, alpha=0.4, linestyle='--')
    # title of the figure by input
    f.suptitle(title, fontsize=26, y=0.95, fontweight='bold')
    # adjust spaces in the figure
    f.subplots_adjust(top=0.7, left=0.15, bottom=0.15, wspace=0.3)
    # if ylim given in input, define it by the input
    if ylim:
        plt.ylim(ylim)
    # save figures as PNG and SVG
    f.savefig(str(numstrains)+title+".png", bbox_inches="tight")
    f.savefig(str(numstrains)+title + ".svg", bbox_inches="tight")
    return f




f_HFF_nucintens = timeplot(df_HFF_inf, ycolumn='mean_intens', stdcolumn='std_intens', sharey=False,
                           title='Host nuclear localization', ylabel='Normalized mean\n' + r'$\alpha$HA fluorescence (AFU)',
                           ylim=False, xticks=8, colorlist=colors)
f_HFF_infection = timeplot(df_HFF, ycolumn='mean_percentinf', stdcolumn='std_percentinf', sharey=True,
                           title='Infection rate', ylabel='Infected cells (%)',
                           ylim=(-10, 100), xticks=8, colorlist=colors)
f_HFF_posnuc = timeplot(df_HFF, ycolumn='mean_percentposnuc', stdcolumn='std_percentposnuc', sharey=True,
                        title='Protein-delivered nuclei', ylabel='Protein-delivered Nnclei (%)',
                        ylim=(-10, 100), xticks=8, colorlist=colors)
f_luhmes_infection = timeplot(df_luhmes, ycolumn='mean_percentinf', stdcolumn='std_percentinf', sharey=True,
                              title='Infection rate', ylabel='Infected cells (%)',
                              ylim=(-10, 100), xticks=6, colorlist=colors)


plt.show()
