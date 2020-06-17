'''
garcia,gil
research script 7 of 7

created: 4/13/2020
last updated: 4/13/2020

purpose: for each set of matches in the multi match list, we plot the position of the ETS, RASS, and CSC
source. Because they are multi matches, there will be 1 ETS source, 1 RASS source
and more than 1 CSC source (varies depending on the match set).

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
rcParams['figure.figsize'] = 10, 10


def lst_of_distinct_values(dataframe):#returns a list of remaining_distinct_values
    remaining_distinct_values = []
    for i in range(dataframe.shape[0]):
        id = dataframe['usrid'].iloc[i]
        if id not in remaining_distinct_values:
            remaining_distinct_values += [id]
    return remaining_distinct_values


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
        count=0
        for i in range(temp_df.shape[0]):
            if count == 0:
                plt.scatter(temp_df['ets_ra'].iloc[i],temp_df['ets_dec'].iloc[i],color='k',s=15)
                plt.scatter(temp_df['rass_ra'].iloc[i],temp_df['rass_dec'].iloc[i],color='k',s=15)
                plt.scatter(temp_df['ETS_CSC_ra_deg'].iloc[i],temp_df['ETS_CSC_dec_deg'].iloc[i],label='CSC2.0 sources',s=15,color='k',marker='*')
                plt.text( x=(temp_df['ets_ra'].iloc[i]),y= (temp_df['ets_dec'].iloc[i])+0.0005,s='ETS'    )
                plt.text(  x=(temp_df['rass_ra'].iloc[i]),y=( temp_df['rass_dec'].iloc[i])+0.0005,s='2RXS'   )
                ets_offset = (temp_df['ets_offset'].iloc[i]) / 3600
                ets_csc_separation = (temp_df['ets_csc_matched_rass_source_sep(arcsec)'].iloc[i] * 1.1) /3600
                if ets_offset > ets_csc_separation:
                    ets_radius = ets_offset
                else:
                    ets_radius = ets_csc_separation
                rass_offset = (temp_df['rass_pos_err'].iloc[i] * 2.5) / 3600
                print('rass offset in arcsec',rass_offset * 3600)
                print('ets_radius in arcsec',ets_radius * 3600)
                circle1 = plt.Circle( (temp_df['ets_ra'].iloc[i], temp_df['ets_dec'].iloc[i]),ets_radius ,color='salmon',alpha=0.4)
                circle2 = plt.Circle( (temp_df['rass_ra'].iloc[i], temp_df['rass_dec'].iloc[i]),rass_offset  ,color='skyblue',alpha=0.4)
                plt.gcf().gca().add_artist(circle1)
                plt.gcf().gca().add_artist(circle2)
                count = count +1

            x_pt = temp_df['ETS_CSC_ra_deg'].iloc[i]
            y_pt = temp_df['ETS_CSC_dec_deg'].iloc[i]
            multiplier=2
            plt.xlim( temp_df['ets_ra'].iloc[i]-(multiplier*ets_offset), temp_df['ets_ra'].iloc[i]+(multiplier*ets_offset)     )
            plt.ylim( temp_df['ets_dec'].iloc[i]-(multiplier*ets_offset), temp_df['ets_dec'].iloc[i]+(multiplier*ets_offset)     )

            plt.scatter(x_pt,y_pt,color='k',s=15,marker='*')
            #plt.text(   x=temp_df['ETS_CSC_ra_deg'].iloc[i],y=temp_df['ETS_CSC_dec_deg'].iloc[i]  +0.0005 ,s='CSC2.0 source'  )
        #plt.legend()
        plt.xlabel('R.A.')
        plt.ylabel('Decl.')
        plt.legend()
        plt.ticklabel_format(useOffset=False)
        plt.show()
        empty_df_w_col_names = pd.concat([empty_df_w_col_names,temp_df])
        seperator_row = (temp_df.shape[1])*['-']
        df_lines = pd.DataFrame([seperator_row],columns=columns_lst)
        empty_df_w_col_names = empty_df_w_col_names.append(df_lines)
    return empty_df_w_col_names




#multi_matches = pd.read_csv('../ets_rass_csc_matches/ETS_CSC_multiple_matches.csv')
multi_matches =  pd.read_csv('../ets_rass_csc_matches/ETS_CSC_multi_multi_matches.csv')

splitting_data(multi_matches)
