'''
garcia,gil
research script 4 of 7
created: 2/16/2020
last updated: 6/4/2020

purpose: from cscview_results_to_rass.py, we outputted a list of ETS-RASS sources
that matched to CSC sources (once or more than once). We take these lists and input
them into the CSCView again, this time using the ETS coords. CSCview outputs the
CSC counterparts to that list. Here, it is imported and joined with the ETS-RASS list.
Again, we create two seperate list based on how many CSC sources matched to the ETS source
(one or more than one). We compute fluxes and flux ratios and filter to find sources that are
considered highly variable.

'''






import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

csc_results = pd.read_csv('../cscview_results_files/csc_matches_to_ETS_pos.csv',sep='\t',header=11)

csc_results.columns = ['usrid','ETS_CSC_separation','ETS_CSC_name','ETS_CSC_ra','ETS_CSC_dec','ETS_CSC_flux_aper_s','ETS_CSC_flux_aper_lolim_s','ETS_CSC_flux_aper_hilim_s','ETS_CSC_flux_aper_m','ETS_CSC_flux_aper_hilim_m','ETS_CSC_flux_aper_lolim_m']


print()
print()
print('----------------------------------')
print('CSC matches using the single matched ETS positions:')
print('----------------------------------')
print()
print()


def angle_seperation(ra1,dec1,ra2,dec2):
    rad= (np.pi/180)
    angle_1 = np.sin(dec1*rad)*np.sin(dec2*rad)
    angle_2 = np.cos(dec1*rad)*np.cos(dec2*rad)
    ra_diff = (ra1)-(ra2)
    angle_3 = np.cos(ra_diff*rad)
    angle = np.arccos(angle_1+(angle_2*angle_3))
    angle_degree=angle*(1/rad)
    angle_arcsec = angle_degree * 3600
    return angle_arcsec



def coords_hr_to_deg(dataframe,HR_ra_column,HR_dec_column,DEG_ra_column,DEG_dec_column):
    ra_deg_lst = []
    dec_deg_lst = []
    for i in range(dataframe.shape[0]):
        ra = dataframe[HR_ra_column].iloc[i]
        dec = dataframe[HR_dec_column].iloc[i]
        coord = SkyCoord(ra+''+dec,unit=(u.hourangle,u.deg))
        ra_deg_lst += [coord.ra.degree]
        dec_deg_lst += [coord.dec.degree]
    dataframe[DEG_ra_column] = ra_deg_lst
    dataframe[DEG_dec_column] = dec_deg_lst
    return dataframe

csc_results = coords_hr_to_deg(csc_results,'ETS_CSC_ra','ETS_CSC_dec','ETS_CSC_ra_deg','ETS_CSC_dec_deg')

ets_matches = pd.read_csv('../ets_rass_csc_matches/RASS_single_matches.csv')
#inner joining csc on ets_rass based on usrid columns
matches = pd.merge(ets_matches,csc_results,how='inner',on='usrid')


print('# of initial matches:\t', matches.shape)

def within_position_error_circ(dataframe,sources_seperation,position_error,mulitiplier=1):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        #position circle
        pos_err = mulitiplier * dataframe[position_error].iloc[i]
        #seperation of ETS and RASS-matched CSC source
        ets_rass_matched_csc_source_seperation = dataframe['ets_csc_matched_rass_source_sep(arcsec)'].iloc[i] * 1.1
        acceptance_radius = pos_err
        if ets_rass_matched_csc_source_seperation > pos_err:
            acceptance_radius = ets_rass_matched_csc_source_seperation
        #seperation defined by cscview of new matches to ETS source
        seperation = dataframe[sources_seperation].iloc[i]
        if acceptance_radius >= seperation:
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    dataframe1 = dataframe[bool_index_lst]
    return dataframe1

matches = within_position_error_circ(matches,'ETS_CSC_separation','ets_offset',mulitiplier=1)


#matches.to_csv('TEMP.csv',index=False)



def no_data(dataframe):
    bool_lst = []
    for i in range(dataframe.shape[0]):
        flux_m = str(dataframe['ETS_CSC_flux_aper_m'].iloc[i])
        flux_s = str(dataframe['ETS_CSC_flux_aper_s'].iloc[i])
        if flux_m.strip()=='' or flux_s.strip()=='':
            bool_lst+=[False]
        elif float(flux_m)==0 or float(flux_s)==0:
            bool_lst+=[False]
        else:
            bool_lst+=[True]
    dataframe1 = dataframe[bool_lst]
    return dataframe1

matches = no_data(matches)
print('# of matches after eliminating those outside error circle & no data:\t', matches.shape)

#to compute 0.5-2.0 keV band fluxes, we add the CSC med band + CSC soft band
def csc_flux(dataframe,med_band,soft_band,new_col_name):
    flux = []
    flux_csc_err_lst = []
    for i in range(dataframe.shape[0]):
        #calculating the flux
        flux_m_flt = float(dataframe[med_band].iloc[i])
        flux_s_flt = float(dataframe[soft_band].iloc[i])
        flux += [flux_s_flt + flux_m_flt]
        #calculating the flux error
        flux_m_err_flt = float(dataframe['ETS_CSC_flux_aper_hilim_m'].iloc[i]) - float(dataframe[med_band].iloc[i])
        flux_s_err_flt = float(dataframe['ETS_CSC_flux_aper_hilim_s'].iloc[i]) - float(dataframe[soft_band].iloc[i])
        # error propogation: err(x+y) = sqrt(dx^2 + dy^2)
        flux_csc_err = np.sqrt(flux_m_err_flt**2 + flux_s_err_flt**2)
        flux_csc_err_lst += [flux_csc_err]
    dataframe[new_col_name] = flux
    dataframe['ETS_CSC_flux_err'] = flux_csc_err_lst
    return dataframe

matches = csc_flux(matches,'ETS_CSC_flux_aper_m','ETS_CSC_flux_aper_s','ETS_CSC_flux')

#we want to consider the single matches and the mulitple matches seperately
#here, we remove and ets-rass sources that matched to more than 1 csc source
#we output a file with only sources that have a single csc source matched to
def single_matched_sources(df):
    #creating a list of unique ids in our table
    usrid_lst = df['usrid'].tolist()
    d = {}
    for i in usrid_lst: d[i] = i in d
    unique = [k for k in d if not d[k]]  # unordered, loop over the dictionary
    #intializing null set, to be filled w boolean
    bool_index = []
    #if the id is unique, it will be kept in the df. else, it is eliminated
    for i in range(df.shape[0]):
        usrid_num = df['usrid'].iloc[i]
        if usrid_num in unique:
            bool_index += [True]
        else:
            bool_index += [False]
    updated_df = df[bool_index]
    return updated_df

single_matches = single_matched_sources(matches)
print('# of single matches', single_matches.shape)

#for single matches only:
#computing flux ratios for the given flux columns
#returns a df of all the flux ratios and another of variable flux ratios
def flux_ratio(dataframe,flux1_column,flux2_column,flux_column_name,variable_flux_restriction=7):
    flux_ratio = (dataframe[flux1_column]) / (dataframe[flux2_column])
    dataframe[flux_column_name] = flux_ratio
    bool_restriction = (dataframe[flux_column_name] > variable_flux_restriction) | (dataframe[flux_column_name] < (1/variable_flux_restriction))
    variable_dataframe = dataframe[bool_restriction]
    return dataframe, variable_dataframe


def flux_ratio_calc(dataframe):
    flux_ratio1 = dataframe['flux_ets'] / dataframe['flux_csc']
    flux_ratio3 = dataframe['flux_rass'] / dataframe['flux_csc']

    flux_ratio1_err = flux_ratio1 * np.sqrt(  (dataframe['flux_ets_err'] / dataframe['flux_ets'])**2 + (dataframe['ETS_CSC_flux_err'] / dataframe['flux_csc'])**2      )
    flux_ratio3_err = flux_ratio3 * np.sqrt(  (dataframe['flux_rass_err'] / dataframe['flux_rass'])**2 + (dataframe['ETS_CSC_flux_err'] / dataframe['flux_csc'])**2      )

    dataframe['ETS_CSC_flux_ratio'] = flux_ratio1
    dataframe['ETS_CSC_flux_ratio_err'] = flux_ratio1_err
    dataframe['2RXS_CSC_flux_ratio'] = flux_ratio3
    dataframe['2RXS_CSC_flux_ratio_err'] = flux_ratio3_err
    return dataframe

single_matches = flux_ratio_calc(single_matches)

def variability(dataframe,variability_threshold):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        ets_2rxs_ratio = dataframe['flux_ratio'].iloc[i]
        ets_csc_ratio = dataframe['ETS_CSC_flux_ratio'].iloc[i]
        twoRXS_csc_ratio = dataframe['2RXS_CSC_flux_ratio'].iloc[i]

        if ets_2rxs_ratio > variability_threshold or ets_2rxs_ratio < (1.0/variability_threshold) or ets_csc_ratio > variability_threshold or ets_csc_ratio < (1.0/variability_threshold) or twoRXS_csc_ratio > variability_threshold or twoRXS_csc_ratio < (1.0/variability_threshold):
        #if ets_csc_ratio > variability_threshold or ets_csc_ratio < (1.0/variability_threshold):
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    variable_single_matches = dataframe[bool_index_lst]
    return variable_single_matches





variable_single_matches = variability(single_matches,7)
print('# of variable single matches:\t',variable_single_matches.shape)


#single_matches, variable_single_matches = flux_ratio(single_matches,'flux_ets','ETS_CSC_flux','ETS_CSC_flux_ratio')
print('# of variable single matches:\t',variable_single_matches.shape)



heeaders_wanted = [  'c1', 'normalized_offset', 'ets_rass_offset', 'flux_ratio', 'f_err',\
       'std_dev_diff', 'flux_ets', 'flux_ets_err', 'flux_rass',\
       'flux_rass_err', 'c11', 'c12', 'c13', 'c14', 'c15', 'c16', 'c17',\
       'ets_ra', 'ets_dec', 'ets_snr', 'c21', 'c22', 'c23', 'c24', 'c25',\
       'c26', 'c27', 'rass_name1', 'rass_name2', 'rass_ra', 'rass_dec',\
       'rass_pos_err', 'rass_l', 'rass_b', 'rass_ML', 'rass_counts',\
       'rass_counts_err', 'rass_count_rate', 'rass_count_rate_err',\
       'rass_exp_time', 'rass_hardness_ratio', 'rass_HR2', 'rass_HR3',\
       'rass_HR4', 'usrid', 'ets_offset', 'separation', 'name', 'ra', 'dec',\
       'flux_aper_m', 'flux_aper_lolim_m', 'flux_aper_hilim_m', 'flux_aper_s',\
       'flux_aper_lolim_s', 'flux_aper_hilim_s', 'csc_ra_deg', 'csc_dec_deg',\
        'flux_csc',\
       'ETS_CSC_flux_err', 'ETS_CSC_separation',\
        'ETS_CSC_flux_ratio', 'ETS_CSC_flux_ratio_err',\
       '2RXS_CSC_flux_ratio', '2RXS_CSC_flux_ratio_err']

names_wanted = [  'c1', 'normalized_offset', 'ets_rass_offset', 'flux_ratio', 'f_err',\
       'std_dev_diff', 'flux_ets', 'flux_ets_err', 'flux_rass',\
       'flux_rass_err', 'c11', 'c12', 'c13', 'c14', 'c15', 'c16', 'c17',\
       'ets_ra', 'ets_dec', 'ets_snr', 'c21', 'c22', 'c23', 'c24', 'c25',\
       'c26', 'c27', 'rass_name1', 'rass_name2', 'rass_ra', 'rass_dec',\
       'rass_pos_err', 'rass_l', 'rass_b', 'rass_ML', 'rass_counts',\
       'rass_counts_err', 'rass_count_rate', 'rass_count_rate_err',\
       'rass_exp_time', 'rass_hardness_ratio', 'rass_HR2', 'rass_HR3',\
       'rass_HR4', 'usrid', 'ets_offset', 'rass_csc_separation', 'csc_name', 'csc_ra', 'csc_dec',\
       'csc_flux_aper_m', 'csc_flux_aper_lolim_m', 'csc_flux_aper_hilim_m', 'csc_flux_aper_s',\
       'csc_flux_aper_lolim_s', 'csc_flux_aper_hilim_s', 'csc_ra_deg', 'csc_dec_deg',\
        'csc_flux',\
       'csc_flux_err', 'ets_csc_separation',\
        'ets_csc_flux_ratio', 'ets_csc_flux_ratio_err',\
       '2rxs_csc_flux_ratio', '2rxs_csc_flux_ratio_err']

print(single_matches.keys())


single_matches.to_csv('../ets_rass_csc_matches/ETS_CSC_single_matches.csv',columns = heeaders_wanted,index=False)
variable_single_matches.to_csv('../ets_rass_csc_matches/ETS_CSC_variable_single_matches.csv',columns=heeaders_wanted,header=names_wanted,index=False)



#now, we want to consider all sources tha DO match to more than 1 CSC source
#we get rid of all sources that match to only 1 csc sources
#we return a df
def remove_singles(df):
    #creating a list of unique ids in our table
    usrid_lst = df['usrid'].tolist()
    d = {}
    for i in usrid_lst: d[i] = i in d
    unique = [k for k in d if not d[k]]  # unordered, loop over the dictionary
    #intializing null set, to be filled w boolean
    bool_index = []
    #if the id is unique, it will ommited from the df. else, it is kept
    for i in range(df.shape[0]):
        usrid_num = df['usrid'].iloc[i]
        if usrid_num in unique:
            bool_index += [False]
        else:
            bool_index += [True]
    updated_df = df[bool_index]
    return updated_df


multiple_matches = remove_singles(matches)

def lst_of_distinct_values(dataframe):#returns a list of remaining_distinct_values
    remaining_distinct_values = []
    for i in range(dataframe.shape[0]):
        id = dataframe['usrid'].iloc[i]
        if id not in remaining_distinct_values:
            remaining_distinct_values += [id]
    return remaining_distinct_values



remaining_distinct_values = lst_of_distinct_values(multiple_matches)

print('# of ets-rass sources that match to many CSC sources:',len(remaining_distinct_values))
print('there are ',len(remaining_distinct_values),' ets-rass sources matching to',matches.shape[0]-single_matches.shape[0],' csc sources.')






def splitting_data(dataframe):
    columns_lst = list(dataframe.columns)
    #added_cols = 'ETS_CSC_largest_sep'
    #columns_lst.append(added_cols)
    empty_df_w_col_names = pd.DataFrame(columns=columns_lst) #will be our df w split data
    remaining_distinct_values = lst_of_distinct_values(dataframe)
    #indexing dataframe frontwards and backwards by each unique value:
    remaining_usrid_lst = dataframe['usrid'].tolist()
    reversed_remaining_usrid_lst = remaining_usrid_lst[::-1]
    #appending to the empty dataframe, adding a seperator in between
    for element in remaining_distinct_values:
        front_index = remaining_usrid_lst.index(element)
        back_index = len(remaining_usrid_lst) - reversed_remaining_usrid_lst.index(element)
        temp_df = dataframe.iloc[front_index:back_index,:]
        max_dist_lst =[]
        if temp_df.shape[0] < 3:
        #for i in range(temp_df.shape[0]):
        #    max_dist_lst += [temp_df['ets_csc_matched_rass_source_sep(arcsec)'].max()]
        #temp_df['ETS_CSC_largest_sep'] = max_dist_lst
            empty_df_w_col_names = pd.concat([empty_df_w_col_names,temp_df])
            seperator_row = (temp_df.shape[1])*['-']
            df_lines = pd.DataFrame([seperator_row],columns=columns_lst)
            empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
    return empty_df_w_col_names



split_multiple_matches = splitting_data(multiple_matches)



multiple_matches.to_csv('../ets_rass_csc_matches/ETS_CSC_multiple_matches.csv',index=False)


split_multiple_matches.to_csv('../ets_rass_csc_matches/ETS_CSC_multiple_matches_split.csv',index=False)





######
# now we deal with the multiple matches
#####


print()
print()
print('----------------------------------')
print('CSC matches using the multi matched ETS positions:')
print('----------------------------------')
print()
print()


csc_results_multi = pd.read_csv('../cscview_results_files/csc_matches_to_ETS_multiple.csv',sep='\t',header=11)
csc_results_multi.columns =  ['usrid','ETS_CSC_separation','ETS_CSC_name','ETS_CSC_ra','ETS_CSC_dec','ETS_CSC_flux_aper_s','ETS_CSC_flux_aper_lolim_s','ETS_CSC_flux_aper_hilim_s','ETS_CSC_flux_aper_m','ETS_CSC_flux_aper_lolim_m','ETS_CSC_flux_aper_hilim_m']

#adding the coordinates in degrees
csc_results_multi = coords_hr_to_deg(csc_results_multi,'ETS_CSC_ra','ETS_CSC_dec','ETS_CSC_ra_deg','ETS_CSC_dec_deg')
#reading in all the ETS sources that matched to many CSC sources and then combining
ets_matches_multi = pd.read_csv('../ets_rass_csc_matches/RASS_multiple_matches.csv')
matches_multi = pd.merge(ets_matches_multi,csc_results_multi,how='inner',on='usrid')

print('# of initial matches:\t', matches_multi.shape)


def within_position_error_circ_multi(dataframe,sources_seperation,position_error,mulitiplier=1):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        #position circle
        pos_err = mulitiplier * dataframe[position_error].iloc[i]
        #seperation of ETS and RASS-matched CSC source
        ets_rass_matched_csc_source_seperation = dataframe['ETS_CSC_largest_sep'].iloc[i] * 1.1
        acceptance_radius = pos_err
        if ets_rass_matched_csc_source_seperation > pos_err:
            acceptance_radius = ets_rass_matched_csc_source_seperation
        #seperation defined by cscview of new matches to ETS source
        seperation = dataframe[sources_seperation].iloc[i]
        if acceptance_radius >= seperation:
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    dataframe1 = dataframe[bool_index_lst]
    return dataframe1

matches_multi = within_position_error_circ_multi(matches_multi,'ETS_CSC_separation','ets_offset',mulitiplier=1)
matches_multi = no_data(matches_multi)

print('# of matches after eliminating those outside error circle & no data:\t', matches_multi.shape)

#computing the fluxes
matches_multi = csc_flux(matches_multi,'ETS_CSC_flux_aper_m','ETS_CSC_flux_aper_s','ETS_CSC_flux')


distinct_values_multi_matches = lst_of_distinct_values(matches_multi)








matches_multi_split = splitting_data(matches_multi)



print('there are ', len(distinct_values_multi_matches), 'ets-rass sources matching to',matches_multi.shape,'csc sources')
#print(matches_multi_split)

matches_multi.to_csv('../ets_rass_csc_matches/ETS_CSC_multi_multi_matches.csv',index=False)
matches_multi_split.to_csv('../ets_rass_csc_matches/ETS_CSC_multi_multi_matches_split.csv',index=False)
