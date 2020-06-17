'''
garcia,gil
research script 6 of 7

created: 2/27/2020
last updated: 4/13/2020


purpose: for all the ETS-RASS sources that matched to a single CSC source, we plot
their position along with their position error circle. Each of these plots contain a
single ETS, RASS, and CSC source since they are singly matched.
'''



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

single_matches=pd.read_csv('../ets_rass_csc_matches/ETS_CSC_single_matches.csv')



for i in range(35,39): #range(single_matches.shape[0]):
    usrid = single_matches['usrid'].iloc[i]
    print()
    print(usrid)
    print()
    #if single_matches['name'].iloc[i] != single_matches['ETS_CSC_name'].iloc[i]:
    ra = [single_matches['ets_ra'].iloc[i], single_matches['rass_ra'].iloc[i], single_matches['csc_ra_deg'].iloc[i], single_matches['ETS_CSC_ra_deg'].iloc[i] ]
    dec = [single_matches['ets_dec'].iloc[i], single_matches['rass_dec'].iloc[i], single_matches['csc_dec_deg'].iloc[i], single_matches['ETS_CSC_dec_deg'].iloc[i] ]

    print('ets_ra,dec:',single_matches['ets_ra'].iloc[i],single_matches['ets_dec'].iloc[i] )
    print('rass ra,dec:', single_matches['rass_ra'].iloc[i], single_matches['rass_dec'].iloc[i] )
    print('rass matched csc:', single_matches['csc_ra_deg'].iloc[i], single_matches['csc_dec_deg'].iloc[i]  )
    print('csc matched csc:', single_matches['ETS_CSC_ra_deg'].iloc[i],single_matches['ETS_CSC_dec_deg'].iloc[i])


    ets_offset = (single_matches['ets_offset'].iloc[i]) / 3600
    ets_csc_separation = (single_matches['ets_csc_matched_rass_source_sep(arcsec)'].iloc[i] * 1.1) /3600
    if ets_offset > ets_csc_separation:
        ets_radius = ets_offset
    else:
        ets_radius = ets_csc_separation
    rass_offset = (single_matches['rass_pos_err'].iloc[i] * 2.5) / 3600
    print('rass offset in arcsec',rass_offset * 3600)
    print('ets_radius in arcsec',ets_radius * 3600)
    circle1 = plt.Circle( (single_matches['ets_ra'].iloc[i], single_matches['ets_dec'].iloc[i]),ets_radius ,color='salmon',alpha=0.4)
    circle2 = plt.Circle( (single_matches['rass_ra'].iloc[i], single_matches['rass_dec'].iloc[i]),rass_offset  ,color='skyblue',alpha=0.4)
    plt.gcf().gca().add_artist(circle1)
    plt.gcf().gca().add_artist(circle2)
    multiplier=2
    #plt.scatter( single_matches['ets_ra'].iloc[i], single_matches['ets_dec'].iloc[i],color ='black',s=15,label='ETS source'  )
    #plt.scatter( single_matches['rass_ra'].iloc[i], single_matches['rass_dec'].iloc[i],color ='black',s=15,label='RASS source'  )
    #plt.scatter(  single_matches['ETS_CSC_ra_deg'].iloc[i], single_matches['ETS_CSC_dec_deg'].iloc[i],color='black',s=15,label='CSC2.0 source'    )
    plt.xlim( single_matches['ets_ra'].iloc[i]-(multiplier*ets_offset), single_matches['ets_ra'].iloc[i]+(multiplier*ets_offset)     )
    plt.ylim( single_matches['ets_dec'].iloc[i]-(multiplier*ets_offset), single_matches['ets_dec'].iloc[i]+(multiplier*ets_offset)     )
    plt.scatter(ra,dec,color='k',s=15)
    plt.text( x=(single_matches['ets_ra'].iloc[i]),y= (single_matches['ets_dec'].iloc[i])+0.001,s='ETS'    )
    plt.text(  x=(single_matches['rass_ra'].iloc[i])-0.001,y=( single_matches['rass_dec'].iloc[i])+0.001,s='2RXS'   )
    plt.text(   x= single_matches['csc_ra_deg'].iloc[i]+0.001,y=single_matches['csc_dec_deg'].iloc[i]+0.001,s='CSC2.0 source'   )
    #plt.text(   x=single_matches['ETS_CSC_ra_deg'].iloc[i],y=single_matches['ETS_CSC_dec_deg'].iloc[i]  +0.002 ,s='ets matched csc'                  )
    #plt.grid()
    plt.xlabel('R.A.')
    plt.ylabel('Decl.')
    #plt.title('usrid: '+str(usrid))
    #plt.legend()
    #plt.legend()
    plt.ticklabel_format(useOffset=False)
    plt.show()
