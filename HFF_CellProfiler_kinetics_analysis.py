import pandas as pd
import os


def label_conditions(df, plate_columns_map, plate_rows_map):
    """Assigning the conditions used in each well, according to the plate maps.
    Receives dataframe, plate_columns_map, plate_rows_map.
    Returns the dataframe."""

    # add strain and moi info
    df['strain'] = df['column'].apply(lambda value: plate_columns_map[value])
    df['moi'] = df['row'].apply(lambda value: plate_rows_map[value])
    # change the strain of all moi=0 wells to null (so that they are grouped together)
    df.ix[df.moi == 0, 'strain'] = 'null'
    # condition is defined as the strain + moi
    df['condition'] = df.strain.str.cat(df.moi.astype(str))
    return df


def calc_image_infection(image_grouped):
    """Calculate infection rate in each image."""
    image_grouped['image_infection'] = image_grouped['infection'].mean()
    return image_grouped


def calc_image_stat(image_grouped, param):
    """
    calculating the image background or mean fluorescence.
    Background fluorescence is defined as as the 10th lower percentile of the uninfected nuclei
    in the image (cells that shouldn't have received any labeled protein because they are not infected).
    The 10th lower percentile was selected for defining the image background by testing different values for
    their ability to accurately represent the image background while maintaining the distribution of cell values.
    Mean fluorescence defined as mean of all cells in image.
    The calculated value will be added in a new column "image_background" or "image_mean".
    """
    if param == 'background':
        image_grouped['image_'+param] = image_grouped.loc[image_grouped['count'] == 0]['nuc_intens'].quantile(0.1)
    if param == 'mean':
        image_grouped['image_'+param] = image_grouped['nuc_intens'].mean()
    return image_grouped


def remove_outliers_image(condition_grouped):
    """
    Remove outlier images which have mean fluorescence >5 std away from the mean of images from the same conditon.
    Function receives the data grouped by the comparison group (condition + timepoint).
    """
    # calculate image mean for each image (already grouped by the comparison group, and grouped by well+field
    # for calculating each image mean)
    condition_grouped = condition_grouped.groupby(['well', 'field']).apply(calc_image_stat, param='mean')
    # calculate the mean and std of mean image fluorescece for all images in the comparison group
    condition_mean, condition_std = condition_grouped['image_mean'].mean() , condition_grouped['image_mean'].std()
    # remove images with mean fluorescence that is >5 std away from the mean of that condition
    condition_grouped = condition_grouped.loc[(condition_mean - condition_grouped['image_mean']).abs() < (5 * condition_std)]
    return condition_grouped


def make_t0(df):
    """
    Make "timepoint 0" data for each condition by copying the group of uninfected wells (moi=0 / strain=null)
    for each moi and strain.
    Receive and return dataframe.
    """
    # copy all cells of uninfected wells to a new dataframe t0_df
    t0_df = df.loc[df['strain'] == 'null'].copy(deep=True)
    # define them as timepoint 0
    t0_df['timepoint'] = 0
    # remove uninfected wells from the original dataframe
    df = df[df['strain'] != 'null']
    # extract list of mois and strains (=conditions) used in the experiment
    mois = sorted(df['moi'].unique())
    strains = df['strain'].unique()
    # copy the "uninfected wells" dataframe for each moi and strain, so it acts as the "timepoint 0" of each condition
    for strain in strains:
        for moi in mois:
            t0_df['moi'] = moi
            t0_df['strain'] = strain
            t0_df['condition'] = str(strain)+str(moi)
            # add it to the big dataframe
            df = df.append(t0_df)
    return df


def per_group_stats(df, group):
    """
    Get per condition stats tables for all cells, only infected cells, 
    and only protein-delivered (nucleus-positive) cells.
    """
    if group == 'condition':
        # condition = unique by moi+timepoint+strain
        statdf_group = df.groupby(['moi', 'timepoint', 'strain'])['nuc_intens'].describe()
        statdf_group.rename(columns={'mean': 'mean_intens'}, inplace=True)
        statdf_group.rename(columns={'std': 'std_intens'}, inplace=True)
        statdf_group.drop(['min', '25%', '50%', '75%', 'max'], axis=1, inplace=True)
        statdf_group2 = df.groupby(['moi', 'timepoint', 'strain'])['pos_nuc'].describe()
        statdf_group2.rename(columns={'mean': 'mean_percentposnuc'}, inplace=True)
        statdf_group2.rename(columns={'std': 'std_percentposnuc'}, inplace=True)
        statdf_group3 = df.groupby(['moi', 'timepoint', 'strain'])['infection'].describe()
        statdf_group3.rename(columns={'mean': 'mean_percentinf'}, inplace=True)
        statdf_group3.rename(columns={'std': 'std_percentinf'}, inplace=True)
        statdf_group = statdf_group.join(statdf_group2[['mean_percentposnuc', 'std_percentposnuc']]).join(
            statdf_group3[['mean_percentinf', 'std_percentinf']])
    elif group == 'well':
        # condition = unique by moi+timepoint+strain+well
        statdf_group = df.groupby(['moi', 'timepoint', 'strain', 'well'])['nuc_intens'].describe()
        statdf_group.rename(columns={'mean': 'mean_intens'}, inplace=True)
        statdf_group.rename(columns={'std': 'std_intens'}, inplace=True)
        statdf_group.drop(['min', '25%', '50%', '75%', 'max'], axis=1, inplace=True)
        statdf_group2 = df.groupby(['moi', 'timepoint', 'strain', 'well'])['pos_nuc'].describe()
        statdf_group2.rename(columns={'mean': 'mean_percentposnuc'}, inplace=True)
        statdf_group2.rename(columns={'std': 'std_percentposnuc'}, inplace=True)
        statdf_group3 = df.groupby(['moi', 'timepoint', 'strain', 'well'])['infection'].describe()
        statdf_group3.rename(columns={'mean': 'mean_percentinf'}, inplace=True)
        statdf_group3.rename(columns={'std': 'std_percentinf'}, inplace=True)
        statdf_group = statdf_group.join(statdf_group2[['mean_percentposnuc', 'std_percentposnuc']]).join(
            statdf_group3[['mean_percentinf', 'std_percentinf']])
    elif group == 'well_conditions':
        # condition = unique by moi+timepoint+strain, but this option receives the dataframe with values averaged per well
        statdf_group = df.groupby(['moi', 'timepoint', 'strain'])['mean_intens'].describe()
        statdf_group.rename(columns={'mean': 'mean_intens'}, inplace=True)
        statdf_group.rename(columns={'std': 'std_intens'}, inplace=True)
        statdf_group.drop(['min', '25%', '50%', '75%', 'max'], axis=1, inplace=True)
        statdf_group2 = df.groupby(['moi', 'timepoint', 'strain'])['mean_percentposnuc'].describe()
        statdf_group2.rename(columns={'mean': 'mean_percentposnuc'}, inplace=True)
        statdf_group2.rename(columns={'std': 'std_percentposnuc'}, inplace=True)
        statdf_group3 = df.groupby(['moi', 'timepoint', 'strain'])['mean_percentinf'].describe()
        statdf_group3.rename(columns={'mean': 'mean_percentinf'}, inplace=True)
        statdf_group3.rename(columns={'std': 'std_percentinf'}, inplace=True)
        statdf_group = statdf_group.join(statdf_group2[['mean_percentposnuc', 'std_percentposnuc']]).join(
            statdf_group3[['mean_percentinf', 'std_percentinf']])
    return statdf_group


# ### import data, label and arrange ###

folder_join = os.path.join
working_dir = r"hff_cellprofiler_outputs"  # working folder

# import csv files, merge sheets into one big dataframe, 
# re-arrange and change column names so they are easier to work with
# NOTE: because I ran the cellprofiler in batches, 
# the data I'm using must be merged from multiple cellprofiler output files
filenames = os.listdir(working_dir)
print(filenames)
# reading and merging
df = pd.concat([pd.read_csv(folder_join(working_dir, file), header=1, usecols=[0, 1, 5, 6, 7, 36, 54]) for file in filenames], ignore_index=True)
# column 0 : ImageNumber (in the batch this image was processed in)
df = df.rename(index=str, columns={"ObjectNumber": "nuc_num",  # column 1 : serial number of the nucleus (= serial number of the cell)
                                   "Metadata_field": "field",  # column 5 : field of view of the image in respect to the well it was recorded from (5 fields per well)
                                   "Metadata_timepoint": "timepoint",  # column 6 : timepoint (hours after infection)
                                   "Metadata_well": "well",  # column 7 : well in the 384-well plate
                                   "Intensity_MeanIntensity_Green": "nuc_intens",  # column 36 : nuc_intens = mean fluorescence intensity in the nucleus, in the green channel (anti-HA stain)
                                   "Children_parasites_Count": "count",  # column 54 : count = number of parasites identified in this cell.
                                   # The CellProfiler analysis associates parasites with human cell nuclei by proximity, meaning that
                                   # its assumed that each identified parasite is inside the same cell as the closest nucleus to the parasite.
                                   })
# extract row and column from 'well'
df['row'], df['column'] = df['well'].str.split(' - ').str
# turn timepoint and count to numeric
df['timepoint'] = pd.to_numeric(df['timepoint'], errors='coerce')
df['count'] = pd.to_numeric(df['count'], errors='coerce')
# 24 hours was mis-labled as 21 hours (because of my handwriting)
df['timepoint'] = df['timepoint'].replace(21, 24)
# A was removed as for some reason there weren't any parasites there (probably human error when manually
# infecting the plates)
df = df[df.row != 'A']

# label wells by strain and moi (according to location in the plate):

# plate arrangement
plate_columns_map = {'01': 'HA', '02': 'HA', '03': 'HA', '04': 'HA', '05': 'HA', '06': 'HA', '07': 'HA', '08': 'HA',
                         '09': 'MECP', '10': 'MECP', '11': 'MECP', '12': 'MECP', '13': 'MECP', '14': 'MECP', '15' :'MECP', '16': 'MECP',
                         '17': 'TFEB', '18': 'TFEB', '19': 'TFEB', '20': 'TFEB', '21': 'TFEB', '22': 'TFEB', '23': 'TFEB', '24': 'TFEB', }
plate_rows_map = {'A': 3, 'B': 3, 'C': 3, 'D': 3, 'E': 3,
                      'F': 1, 'G': 1, 'H': 1, 'I': 1, 'J': 1,
                      'K': 0.33, 'L': 0.33, 'M': 0.33, 'N': 0.33, 'O': 0.33,
                      'P': 0}

df = label_conditions(df, plate_columns_map, plate_rows_map)

# define infected cells = cells infected with at least one parasite
df.ix[df['count'] > 0, 'infection'] = 100
df['infection'].fillna(value=0, inplace=True)

# ### clean data from artifacts, normalize fluorescence, and define timepoint 0 ###

# record number of cells in the initial dataset (before cleaning from artefacts and outliers)
numcells = df.shape[0]
print(str(numcells) + ' total cells in dataset')

# remove images of wells which were inoculated with parasites (moi>0 / strain not null),
# but in which no parasites were identified (assumed to be a problem in the image analysis, resulting from
#  parasite staining being too weak in some images)
df = df.groupby(['condition', 'timepoint', 'well', 'field']).apply(calc_image_infection)
df = pd.concat([df.loc[df['moi'] == 0], df.loc[(df['moi'] != 0) & (df['image_infection'] > 0)]])

# report how many cells were removed from the dataset in this step
noparasite_cells = numcells - df.shape[0]
print(str(noparasite_cells)+' cells ('+str(noparasite_cells/numcells*100) +
      '%) removed due to images of inoculated wells with no identified parasites - insufficient immunostaining')

#  calc per image background. Images can be defined by unique combination of condition+timepoint+well+field.
df = df.groupby(['condition', 'timepoint', 'well', 'field']).apply(calc_image_stat, param='background')

# normalize nuclear background fluorescence value by subtracting the image background,
# and multiply by 10000 for ease of plotting (anyway units are arbitrary)
df['nuc_intens'] = (df['nuc_intens'] - df['image_background']) * 10000

# remove images with outlier background fluorescence (in comparison to all images in that condition)
df = df.groupby(['condition', 'timepoint']).apply(remove_outliers_image)
# report how many cells were removed from the dataset in this step
outliercells = numcells - noparasite_cells - df.shape[0]
print(str(outliercells)+' cells ('+str(outliercells/numcells*100) + '%) removed due to wells with outlier mean fluorescence - fluorescence artefact')

# define nuclei that received the labeled protein (positive nuclei) as nuclei with fluorescence above the
# 0.99 quantile of uninfected cells (cells from moi=0 / strain=null). The 0.99 quantile of uninfected cell
# nuclei was selected as a strict definition for "uninfected nuclei fluorescence" (backgroud/nuclei negative
# for the protein). 0.99 was selected instead of "max" in order to avoid defining this value by the outliers.
uninfectednucintens = df.loc[df['strain'] == 'null']['nuc_intens'].quantile(0.99)
df.ix[df['nuc_intens'] > uninfectednucintens, 'pos_nuc'] = 100
df['pos_nuc'].fillna(value=0, inplace=True)

# make "timepoint 0" data for each condition by copying the group of uninfected wells (moi=0 / strain=null)
# for each moi and strain
df = make_t0(df)

df.reset_index(drop=True, inplace=True)

# ### saving results and generating statistical summary tables ###

# save processed df as csv file
df.to_csv('cellprofiler_HFF_sorted_25.11.csv')

# for the statistical dataframes, remove timepoints>24, because after that vacuoles have became too big,
# making image segmentation too inaccurate
df = df.loc[df['timepoint'] < 25]

# define df of only single-infected and only positive nuclei cells
df_pos = df.loc[df['pos_nuc'] == 100]
df_inf = df.loc[df['count'] == 1]

# for the "only single-infected" and "only positive nuclei" dataframes, add the timepoint 0 from the "all cells"
# dataframe(because there weren't really infected or positive-nuclei cells in timepoint 0)
df_pos = pd.concat([df.loc[df['timepoint'] == 0], df_pos.loc[df_pos['timepoint'] != 0]])
df_inf = pd.concat([df.loc[df['timepoint'] == 0], df_inf.loc[df_inf['timepoint'] != 0]])


# per condition stats, by data from all cells in the dataset
statdf_HFF_condition_all = per_group_stats(df, 'condition')
# per condition stats, by data from only infected cells in the dataset
statdf_HFF_condition_inf = per_group_stats(df_inf, 'condition')
# per condition stats, by data from only positive-nuclei cells in the dataset
statdf_HFF_condition_pos = per_group_stats(df_pos, 'condition')

# save stat tables as csv files
statdf_HFF_condition_all.to_csv('statdf_HFF_condition_all_25.11.csv')
statdf_HFF_condition_inf.to_csv('statdf_HFF_condition_inf_25.11.csv')
statdf_HFF_condition_pos.to_csv('statdf_HFF_condition_pos_25.11.csv')


# per well stats, by data from all cells in the dataset
df_HFF_perimage_all = per_group_stats(df, 'well')
# per well stats, by data from only infected cells in the dataset
df_HFF_perimage_inf = per_group_stats(df_inf, 'well')
# per well stats, by data from only positive-nuclei cells in the dataset
df_HFF_perimage_pos = per_group_stats(df_pos, 'well')

# save stat tables as csv files
df_HFF_perimage_all.to_csv('df_HFF_perimage_all_25.11.csv')
df_HFF_perimage_inf.to_csv('df_HFF_perimage_inf_25.11.csv')
df_HFF_perimage_pos.to_csv('df_HFF_perimage_pos_25.11.csv')

# per condition stats, by data of per well averages (for all cells in the well)
statdf_HFF_perimage_all = per_group_stats(df_HFF_perimage_all, 'well_conditions')
# per condition stats, by data of per well averages (for infected cells in the well)
statdf_HFF_perimage_inf = per_group_stats(df_HFF_perimage_inf, 'well_conditions')
# per condition stats, by data of per well averages (for positive-nuclei cells in the well)
statdf_HFF_perimage_pos = per_group_stats(df_HFF_perimage_pos, 'well_conditions')

# save stat tables as csv files
statdf_HFF_perimage_all.to_csv('statdf_HFF_perimage_all_25.11.csv')
statdf_HFF_perimage_inf.to_csv('statdf_HFF_perimage_inf_25.11.csv')
statdf_HFF_perimage_pos.to_csv('statdf_HFF_perimage_pos_25.11.csv')

# End of script