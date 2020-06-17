'''
garcia,gil
research script 5 of 7

created: 3/5/2020
last updated: 4/14/2020

purpose: from the list of ETS matches the had no RASS counterpart, we look for
CSC counterparts in the CSCview. Here, we merge the two files and compute flux ratios
to find high variable sources.

'''



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord


csc_results = pd.read_csv('../cscview_results_files/csc_matches_to_ets_no_matches.csv',sep='\t',header=11)
csc_results.columns = ['usrid','separation','name','ra','dec','flux_aper_m','flux_aper_lowlim_m','flux_aper_hilim_m','flux_aper_s','flux_aper_lolim_s','flux_aper_hilim_s']




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

csc_results = coords_hr_to_deg(csc_results,'ra','dec','ra_deg','dec_deg')

ets_no_match = pd.read_csv('../ets_rass_matches_files/no_match_theta.csv')
#inner joining csc on ets_rass based on usrid columns
matches = pd.merge(ets_no_match,csc_results,how='inner',on='usrid')

print('# of initial matches:\t', matches.shape)


def within_position_error_circ(dataframe,sources_seperation,position_error,mulitiplier=1):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        #position circle
        pos_err = mulitiplier * dataframe[position_error].iloc[i]
        #seperation defined by cscview of new matches to ETS source
        seperation = dataframe[sources_seperation].iloc[i]
        if pos_err >= seperation:
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    dataframe1 = dataframe[bool_index_lst]
    return dataframe1


matches = within_position_error_circ(matches,'separation','ets_offset',mulitiplier=1)


print('# of matches within 1sigma of ETS circle:\t',matches.shape)


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
        flux_m_err_flt = float(dataframe['flux_aper_hilim_m'].iloc[i]) - float(dataframe[med_band].iloc[i])
        flux_s_err_flt = float(dataframe['flux_aper_hilim_s'].iloc[i]) - float(dataframe[soft_band].iloc[i])
        # error propogation: err(x+y) = sqrt(dx^2 + dy^2)
        flux_csc_err = np.sqrt(flux_m_err_flt**2 + flux_s_err_flt**2)
        flux_csc_err_lst += [flux_csc_err]
    dataframe[new_col_name] = flux
    dataframe['ETS_CSC_flux_err'] = flux_csc_err_lst
    return dataframe

matches = csc_flux(matches,'flux_aper_m','flux_aper_s','csc_flux')


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



def flux_ratio(dataframe,flux1_column,flux2_column,flux_column_name,variable_flux_restriction=7):
    flux_ratio = (dataframe[flux1_column]) / (dataframe[flux2_column])
    dataframe[flux_column_name] = flux_ratio
    bool_restriction = (dataframe[flux_column_name] > variable_flux_restriction) | (dataframe[flux_column_name] < (1/variable_flux_restriction))
    variable_dataframe = dataframe[bool_restriction]
    return dataframe, variable_dataframe


def flux_ratio_calc(dataframe):
    flux_ratio1 = dataframe['ets_flux'] / dataframe['csc_flux']
    #flux_ratio3 = dataframe['flux_rass'] / dataframe['csc_flux']

    flux_ratio1_err = flux_ratio1 * np.sqrt(  (dataframe['ets_flux_err'] / dataframe['ets_flux'])**2 + (dataframe['ETS_CSC_flux_err'] / dataframe['csc_flux'])**2      )
    #flux_ratio3_err = flux_ratio3 * np.sqrt(  (dataframe['flux_rass_err'] / dataframe['flux_rass'])**2 + (dataframe['csc_flux'] / dataframe['ETS_CSC_flux_err'])**2      )

    dataframe['ETS_CSC_flux_ratio'] = flux_ratio1
    dataframe['ETS_CSC_flux_ratio_err'] = flux_ratio1_err
    #dataframe['2RXS_CSC_flux_ratio'] = flux_ratio3
    #dataframe['2RXS_CSC_flux_ratio_err'] = flux_ratio3_err
    return dataframe

single_matches = flux_ratio_calc(single_matches)

def variability(dataframe,variability_threshold):
    bool_index_lst=[]
    for i in range(dataframe.shape[0]):
        #ets_2rxs_ratio = dataframe['flux_ratio'].iloc[i]
        ets_csc_ratio = dataframe['ETS_CSC_flux_ratio'].iloc[i]
        #twoRXS_csc_ratio = dataframe['2RXS_CSC_flux_ratio'].iloc[i]

#    if ets_2rxs_ratio > variability_threshold or ets_2rxs_ratio < (1.0/variability_threshold) or ets_csc_ratio > variability_threshold or ets_csc_ratio < (1.0/variability_threshold) or twoRXS_csc_ratio > variability_threshold or twoRXS_csc_ratio < (1.0/variability_threshold):
        if ets_csc_ratio > variability_threshold or ets_csc_ratio < (1.0/variability_threshold):
            bool_index_lst += [True]
        else:
            bool_index_lst += [False]
    variable_single_matches = dataframe[bool_index_lst]
    return variable_single_matches


print(single_matches.keys())

columns_wanted = ['c1', 'c2', 'c3', 'c4', 'ets_flux', 'ets_flux_err', 'ets_ra', 'ets_dec',\
       'ets_snr', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15', 'c16', 'usrid',\
       'ets_offset', 'separation', 'name', 'ra', 'dec', 'flux_aper_m',\
       'flux_aper_lowlim_m', 'flux_aper_hilim_m', 'flux_aper_s',\
       'flux_aper_lolim_s', 'flux_aper_hilim_s', 'ra_deg', 'dec_deg',\
       'csc_flux', 'ETS_CSC_flux_err', 'ETS_CSC_flux_ratio',\
       'ETS_CSC_flux_ratio_err']

names_wanted = ['c1', 'c2', 'c3', 'c4', 'ets_flux', 'ets_flux_err', 'ets_ra', 'ets_dec',\
       'ets_snr', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15', 'c16', 'usrid',\
       'ets_offset', 'separation', 'csc_name', 'csc_ra', 'csc_dec', 'csc_flux_aper_m',\
       'csc_flux_aper_lowlim_m', 'csc_flux_aper_hilim_m', 'csc_flux_aper_s',\
       'csc_flux_aper_lolim_s', 'csc_flux_aper_hilim_s', 'csc_ra_deg', 'csc_dec_deg',\
       'csc_flux', 'csc_flux_err', 'ets_csc_flux_ratio',\
       'ets_csc_flux_ratio_err']




variable_single_matches = variability(single_matches,7)
print('# of variable single matches:\t',variable_single_matches.shape)

variable_single_matches.to_csv('../ets_rass_csc_matches/no_match_theta_variable_single_matches.csv',index=False)
single_matches.to_csv('../ets_rass_csc_matches/no_match_theta_single_matches.csv',index=False)
