'''
garcia, gil
research script 2 of 7

created: 3/5/2020
last updated 3/5/2020

purpose: adds a usrid and a offset value to the "no match sources" (ETS sources that did
not have a RASS counterpart)
'''



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('../ets_rass_matches_files/no_match_theta.csv')

#creating a usrid column (setting up an inner join)
usrid_lst =[]
for i in range(df.shape[0]):
    usrid_lst += [i+1]
df['usrid'] = usrid_lst


def ets_offsets_from_snr(dataframe,snr_column):
    offset = []
    for i in range(dataframe.shape[0]):
        snr = dataframe[snr_column].iloc[i]
        if snr <=3:
            offset += [73]
        elif snr > 3 and snr <=4:
            offset +=[53]
        elif snr >4 and snr <=6:
            offset +=[41]
        elif snr > 6:
            offset +=[30]
    dataframe['ets_offset'] = offset
    return dataframe

df = ets_offsets_from_snr(df,'ets_snr')


df.to_csv('../ets_rass_matches_files/no_match_theta.csv',index=False)
