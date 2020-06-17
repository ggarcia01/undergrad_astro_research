'''
garcia,gil
research script 3 of 7
created: 1/20/2020
last updated: 6/4/2020

purpose: we have inputted the ETS-RASS catalog into the CSCView using the RASS coords, which outputted
a .csv file contained all the CSC matches to the ETS-RASS catalog. Here, we input the
ETS-RASS catalog and the CSC matches in order to join them. We create a list of all ETS-RASS
sources that match to 1 CSC source. We create a seperate list of all ETS-RASS sources that
matched to more than 1 CSC source.

in the proccess of doing so, we compute and add to the DF:
-convert coordinates from hours to degrees
-compute CSC flux and flux errors
-compute angular seperation between the CSC source and the ETS source.

'''





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
        rass_pos_err = 2.5* dataframe['rass_pos_err'].iloc[i]
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
#we also compute the error associated with this flux observation
def csc_flux(dataframe,med_band,soft_band,new_col_name):
    flux = []
    flux_csc_err_lst = []
    for i in range(dataframe.shape[0]):
        #calculating the flux
        flux_m_flt = float(dataframe[med_band].iloc[i])
        flux_s_flt = float(dataframe[soft_band].iloc[i])
        flux += [flux_s_flt + flux_m_flt]
        #calculating the flux error
        flux_m_err_flt = float(dataframe['flux_aper_hilim_m'].iloc[i]) - float(dataframe[med_band].iloc[i])
        flux_s_err_flt = float(dataframe['flux_aper_hilim_s'].iloc[i]) - float(dataframe[soft_band].iloc[i])
        # error propogation: err(x+y) = sqrt(dx^2 + dy^2)
        flux_csc_err = np.sqrt(flux_m_err_flt**2 + flux_s_err_flt**2)
        flux_csc_err_lst += [flux_csc_err]
    dataframe[new_col_name] = flux
    dataframe['ETS_CSC_flux_err'] = flux_csc_err_lst
    return dataframe

#we want to consider the single matches and the mulitple matches seperately
#here, we remove any ets-rass sources that matched to more than 1 csc source
#this fxn outputs a file with sources that only have a single csc source matched to
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
#this fxn adds a row of "-" in between each set of multi matches
#returns a more readiable df
def splitting_data(dataframe):
    columns_lst = list(dataframe.columns)
    added_cols = 'ETS_CSC_largest_sep'
    columns_lst.append(added_cols)
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
        for i in range(temp_df.shape[0]):
            max_dist_lst += [temp_df['ets_csc_matched_rass_source_sep(arcsec)'].max()]
        temp_df['ETS_CSC_largest_sep'] = max_dist_lst
        empty_df_w_col_names = pd.concat([empty_df_w_col_names,temp_df])
        seperator_row = (temp_df.shape[1])*['-']
        df_lines = pd.DataFrame([seperator_row],columns=columns_lst)
        empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
    return empty_df_w_col_names

#returns angle separation in arcseconds
def angle_seperation(dataframe,ra1_col,dec1_col,ra2_col,dec2_col,new_col_name):
    seperations = []
    for i in range(dataframe.shape[0]):
        ra1,dec1 = dataframe[ra1_col].iloc[i],dataframe[dec1_col].iloc[i]
        ra2,dec2 = dataframe[ra2_col].iloc[i],dataframe[dec2_col].iloc[i]
        rad= (np.pi/180)
        angle_1 = np.sin(dec1*rad)*np.sin(dec2*rad)
        angle_2 = np.cos(dec1*rad)*np.cos(dec2*rad)
        ra_diff = (ra1)-(ra2)
        angle_3 = np.cos(ra_diff*rad)
        angle = np.arccos(angle_1+(angle_2*angle_3))
        angle_degree=angle*(1/rad)
        angle_arcsec = angle_degree * 3600
        seperations += [angle_arcsec]
        #print(  dataframe['usrid'].iloc[i]  ,angle_arcsec)
    dataframe[new_col_name] = seperations
    return dataframe

#this fxn finds the unique values in a column and/or list
def unique_values(dataframe,col_name):
    id_lst = dataframe[col_name].tolist()
    unique_lst = []
    bool_lst = []
    for element in id_lst:
        if element not in unique_lst:
            unique_lst += [element]
            bool_lst += [True]
        else:
            bool_lst += [False]
    dataframe1 = dataframe[bool_lst]
    return dataframe1



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
#removing all csc matches that are farther than 2-times the pos_err_circ of the RASS source
matches = within_position_error_circ(matches)
print('shape after eliminating sources outside of error circle: ',matches.shape)
#removing csc sources with no flux data. cant compute fluxes and flux ratios otherwise
matches = no_data(matches)
print('# of matches after eliminating empty data:',matches.shape)
#angular separation between ETS pos and CSC matched RASS source:
matches = angle_seperation(matches,'ets_ra','ets_dec','csc_ra_deg','csc_dec_deg','ets_csc_matched_rass_source_sep(arcsec)')
#computing csc fluxes and adding them to the dataframe
matches.to_csv('../ets_rass_csc_matches/RASS_all_matches.csv',index=False)
matches = csc_flux(matches,'flux_aper_m','flux_aper_s','flux_csc')
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
multiple_matches_split = splitting_data(multiple_matches)
multiple_matches_split.to_csv('../ets_rass_csc_matches/RASS_CSC_multi_matches_split.csv',index=False)
#list of ets-rass sources w mulitple matches
remaining_distinct_values = lst_of_distinct_values(multiple_matches)
#reporing stats
print('# of ets-rass sources that match to many CSC sources:',len(remaining_distinct_values))
print('there are ',len(remaining_distinct_values),' ets-rass sources matching to',matches.shape[0]-single_matches.shape[0],' csc sources.')
#split our data for easier readability
split_multiple_matches = splitting_data(multiple_matches)
single_split_multiple_matches = unique_values(split_multiple_matches,'usrid')





#export multi matches file
single_split_multiple_matches.to_csv('../ets_rass_csc_matches/RASS_multiple_matches.csv',index=False)


print()
print('done.')
