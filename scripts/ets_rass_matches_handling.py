import pandas as pd

df = pd.read_csv('../ets_rass_matches_files/match_theta_copy.csv')
size_before = df.shape


#creating a usrid column (setting up an inner join)
usrid_lst =[]
for i in range(df.shape[0]):
    usrid_lst += [i+1]
df['usrid'] = usrid_lst

#normalized offset restriction: only want sources w norm_offset < 1.4
restriction1 = df['normalized_offset'] < 1.4
df = df[restriction1]
size_after = df.shape

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


#overwrite the data file:
df.to_csv('../ets_rass_matches_files/match_theta_copy.csv',index=False)

print('done.')
