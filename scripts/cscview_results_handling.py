import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

#importing results from csc_view using CSC 2.0
csc_results = pd.read_csv('../cscview_results_files/csc_matches_to_ets_rass_sources.csv',sep='\t',header=11)
print('total amount of csc matches:',csc_results.shape)

#converting csc coordinates from hour to deg
def coords_hr_to_deg(dataframe):
    ra_deg_lst = []
    dec_deg_lst = []
    for i in range(dataframe.shape[0]):
        ra = dataframe['ra'].iloc[i]
        dec = dataframe['dec'].iloc[i]
        coord = SkyCoord(ra+''+dec,unit=(u.hourangle,u.deg))
        ra_deg_lst += [coord.ra.degree]
        dec_deg_lst += [coord.dec.degree]
    dataframe['csc_ra_deg'] = ra_deg_lst
    dataframe['csc_dec_deg'] = dec_deg_lst
    return dataframe

csc_results = coords_hr_to_deg(csc_results)


#importing ets-rass sources
ets_rass = pd.read_csv('../ets_rass_matches_files/match_theta_copy.csv')

#left joining csc on ets_rass based on usrid columns
matches = pd.merge(ets_rass,csc_results,how='inner',on='usrid')
print('# of matches after initial merge:',matches.shape)

matches.to_csv('TEMP_matches.csv',index=False)

#keeping significant matches: want seperation of sources to be within RASS_err circle
#RASS_pos_err < seperation of ets-rass source w CSC source
def within_position_error_circ(dataframe):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        rass_pos_err = 2* dataframe['rass_pos_err'].iloc[i]
        seperation = dataframe['separation'].iloc[i]
        if rass_pos_err > seperation:
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    dataframe1 = dataframe[bool_index_lst]
    return dataframe1

matches = within_position_error_circ(matches)

print('shape after eliminating sources outside of error circle: ',matches.shape)

def no_data(dataframe):
    bool_lst = []
    for i in range(dataframe.shape[0]):
        flux_m = str(dataframe['flux_aper_m'].iloc[i])
        flux_s = str(dataframe['flux_aper_s'].iloc[i])
        if flux_m.strip()=='' or flux_s.strip()=='':
            bool_lst+=[False]
        elif float(flux_m)==0 or float(flux_s)==0:
            bool_lst+=[False]
        else:
            bool_lst+=[True]
    dataframe1 = dataframe[bool_lst]
    return dataframe1

matches = no_data(matches)
print('# of matches after eliminating empty data:',matches.shape)


def csc_flux(dataframe,med_band,soft_band):
    flux_ratio = []
    for i in range(dataframe.shape[0]):
        flux_m_flt = float(dataframe[med_band].iloc[i])
        flux_s_flt = float(dataframe[soft_band].iloc[i])
        flux_ratio += [flux_s_flt + flux_m_flt]
    dataframe['flux_csc'] = flux_ratio
    return dataframe

matches = csc_flux(matches,'flux_aper_m','flux_aper_s')


matches.to_csv('../ets_rass_csc_matches/all_matches.csv',index=False)




'''
#########
#making a dataframe for only the single matched sources
#########
'''

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

def big_flux_ratio(dataframe,flux1_column,flux2_column):
    flux_ratio = (dataframe[flux1_column]) / (dataframe[flux2_column])
    dataframe['rass_csc_flux_ratio'] = flux_ratio
    bool_restriction = (dataframe['rass_csc_flux_ratio'] > 10) | (dataframe['rass_csc_flux_ratio'] < 0.1)
    variable_dataframe = dataframe[bool_restriction]
    return dataframe, variable_dataframe

single_matches,variable_single_matches = big_flux_ratio(single_matches,'flux_rass','flux_csc')

print('# of single matches:',single_matches.shape)
print('# of variable single matches:',variable_single_matches.shape)

single_matches.to_csv('../ets_rass_csc_matches/single_matches.csv',index=False)
variable_single_matches.to_csv('../ets_rass_csc_matches/variable_single_matches.csv',index=False)

'''
############
# making a dataframe for multiple matches only
############
'''

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


#splitting multi match file by source matches:
def splitting_data(dataframe):
    columns_lst = list(dataframe.columns)
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
        empty_df_w_col_names = pd.concat([empty_df_w_col_names,temp_df])
        df_lines = pd.DataFrame([['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-']],columns=columns_lst)
        empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
        '''
        for i in range(temp_df.shape[0]):
            empty_df_w_col_names = pd.concat([empty_df_w_col_names,temp_df])
            df_lines = pd.DataFrame([['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-']],columns=columns_lst)
            empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
        '''
    return empty_df_w_col_names

split_multiple_matches = splitting_data(multiple_matches)
#print(split_multiple_matches)
split_multiple_matches.to_csv('../ets_rass_csc_matches/split_multiple_matches.csv',index=False)





print()
print('Three .csv files exported.')
print('done.')
