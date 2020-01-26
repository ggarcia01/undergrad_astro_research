import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

'''
defining fxns
'''

#converting csc coordinates from hour to deg
#pass dataframe and the column names for the ra and dec columns, must be in hour units
#then, choose the names of the converted ra and dec columns that are in degrees now
#returns the same dataframe imported but w 2 new columns: ra and dec in degrees
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

#some csc sources have no flux data.
#we remove any csc sources who have a 0 or an empty flux cell
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

#to compute 0.5-2.0 keV band fluxes, we add the CSC med band + CSC soft band
def csc_flux(dataframe,med_band,soft_band):
    flux_ratio = []
    for i in range(dataframe.shape[0]):
        flux_m_flt = float(dataframe[med_band].iloc[i])
        flux_s_flt = float(dataframe[soft_band].iloc[i])
        flux_ratio += [flux_s_flt + flux_m_flt]
    dataframe['flux_csc'] = flux_ratio
    return dataframe

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

#for single matches only:
#computing flux ratios for the given flux columns
#returns a df of all the flux ratios and another of variable flux ratios
def flux_ratio(dataframe,flux1_column,flux2_column,flux_column_name,variable_flux_restriction=10):
    flux_ratio = (dataframe[flux1_column]) / (dataframe[flux2_column])
    dataframe[flux_column_name] = flux_ratio
    bool_restriction = (dataframe[flux_column_name] > variable_flux_restriction) | (dataframe[flux_column_name] < (1/variable_flux_restriction))
    variable_dataframe = dataframe[bool_restriction]
    return dataframe, variable_dataframe

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

#we want to know what ets-rass sources remain.
#we return a list of them
def lst_of_distinct_values(dataframe):#returns a list of remaining_distinct_values
    remaining_distinct_values = []
    for i in range(dataframe.shape[0]):
        id = dataframe['usrid'].iloc[i]
        if id not in remaining_distinct_values:
            remaining_distinct_values += [id]
    return remaining_distinct_values

#splitting multi match file by source matches:
#returns a df
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
        seperator_row = dataframe.shape[1]*['-']
        df_lines = pd.DataFrame([seperator_row],columns=columns_lst)
        empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
    return empty_df_w_col_names


'''
#handling ets-rass and CSC data
#in this instance, we are using the RASS positions + RASS poss errors to find CSC sources
'''

print()
print()
print('----------------------------------')
print('CSC matches using RASS positions:')
print('----------------------------------')
print()
print()

#importing results from csc_view using CSC 2.0
csc_results = pd.read_csv('../cscview_results_files/csc_matches_to_ets_rass_sources.csv',sep='\t',header=11)
print('total amount of csc matches:',csc_results.shape)
#importing ets-rass sources
ets_rass = pd.read_csv('../ets_rass_matches_files/match_theta_copy.csv')
#converting csc coords to degrees
csc_results = coords_hr_to_deg(csc_results,'ra','dec','csc_ra_deg','csc_dec_deg')
#inner joining csc on ets_rass based on usrid columns
matches = pd.merge(ets_rass,csc_results,how='inner',on='usrid')
print('# of matches after initial merge:',matches.shape)
#exporting all matches
matches.to_csv('../ets_rass_csc_matches/RASS_all_matches.csv',index=False)
#removing all csc matches that are farther than 2-times the pos_err_circ of the RASS source
matches = within_position_error_circ(matches)
print('shape after eliminating sources outside of error circle: ',matches.shape)
#removing csc sources with no flux data. cant compute fluxes and flux ratios otherwise
matches = no_data(matches)
print('# of matches after eliminating empty data:',matches.shape)
#computing csc fluxes and adding them to the dataframe
matches = csc_flux(matches,'flux_aper_m','flux_aper_s')
#making a dataframe for only the single matched sources
single_matches = single_matched_sources(matches)
#computing flux ratios, finding all variable sources
single_matches,variable_single_matches = flux_ratio(single_matches,'flux_rass','flux_csc','rass_csc_flux_ratio')
#returning stats
print('# of single matches:',single_matches.shape)
print('# of variable single matches:',variable_single_matches.shape)
#exporting all single matches and all variable sources
single_matches.to_csv('../ets_rass_csc_matches/RASS_single_matches.csv',index=False)
variable_single_matches.to_csv('../ets_rass_csc_matches/RASS_variable_single_matches.csv',index=False)
# making a dataframe for multiple matches only
multiple_matches = remove_singles(matches)
#list of ets-rass sources w mulitple matches
remaining_distinct_values = lst_of_distinct_values(multiple_matches)
#reporing stats
print('# of ets-rass sources that match to many CSC sources:',len(remaining_distinct_values))
print('there are ',len(remaining_distinct_values),' ets-rass sources matching to',matches.shape[0]-single_matches.shape[0],' csc sources.')
#split our data for easier readability
split_multiple_matches = splitting_data(multiple_matches)
#export multi matches file
split_multiple_matches.to_csv('../ets_rass_csc_matches/RASS_multiple_matches.csv',index=False)




print()
print()
print('----------------------------------')
print('CSC matches using ETS positions:')
print('----------------------------------')
print()
print()
'''
#handling ets-rass and CSC data
#in this instance, we are using the ETS positions + ETS poss errors to find CSC sources
'''
csc_matches_to_ets = pd.read_csv('../cscview_results_files/csc_matches_to_ETS_pos.csv',sep='\t',header=12)
print('total # of csc matches:',csc_matches_to_ets.shape)
#converting coords to degrees
csc_matches_to_ets = coords_hr_to_deg(csc_matches_to_ets,'ra','dec','csc_ra_deg','csc_dec_deg')
#inner joining csc on ets_rass based on usrid column
matches1 = pd.merge(ets_rass,csc_matches_to_ets,how='inner',on='usrid')
print('# of matches after initial merge:',matches1.shape)
#exporting all matches
matches1.to_csv('../ets_rass_csc_matches/ETS_all_matches.csv',index=False)
#removing csc sources with no flux data. cant compute fluxes and flux ratios otherwise
matches1 = no_data(matches1)
print('# of matches after eliminating empty data:',matches1.shape)
#computing csc fluxes and adding them to the dataframe
matches1 = csc_flux(matches1,'flux_aper_m','flux_aper_s')
#exporting remaining matches with csc flux calculations
matches1.to_csv('../ets_rass_csc_matches/ETS_all_matches_w_fluxes.csv',index=False)
#making a dataframe for only the single matched sources
single_matches1 = single_matched_sources(matches1)
#computing flux ratios, finding all variable sources
single_matches1,variable_single_matches1 = flux_ratio(single_matches1,'flux_ets','flux_csc','ets_csc_flux_ratio')
#returning stats
print('# of single matches:',single_matches1.shape)
print('# of variable single matches:',variable_single_matches1.shape)
#exporting all single matches and all variable sources
single_matches1.to_csv('../ets_rass_csc_matches/ETS_single_matches.csv',index=False)
variable_single_matches1.to_csv('../ets_rass_csc_matches/ETS_variable_single_matches.csv',index=False)
# making a dataframe for multiple matches only
multiple_matches1 = remove_singles(matches1)
#list of ets-rass sources w mulitple matches
remaining_distinct_values1 = lst_of_distinct_values(multiple_matches1)
#reporing stats
print('# of ets-rass sources that match to many CSC sources:',len(remaining_distinct_values1))
print('there are ',len(remaining_distinct_values1),' ets-rass sources matching to',matches1.shape[0]-single_matches1.shape[0],' csc sources.')
#split our data for easier readability
split_multiple_matches1 = splitting_data(multiple_matches1)
#export multi matches file
split_multiple_matches1.to_csv('../ets_rass_csc_matches/ETS_multiple_matches.csv',index=False)

print()
print('done.')
