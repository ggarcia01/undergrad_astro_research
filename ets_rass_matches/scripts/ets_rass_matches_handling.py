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
#overwrite the data file:
df.to_csv('../ets_rass_matches_files/match_theta_copy.csv',index=False)

print('done.')
