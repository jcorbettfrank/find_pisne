
# coding: utf-8

# In[55]:


import requests
import os
import pandas as pd
import numpy as np
from io import StringIO
import errno
import pdb
import string

# In[56]:

OverArchDirs = ['Data','Data_BVRI']

SN_types = ['SLSN','1a','IIn','IIp','IIb']
for OverArchDir in OverArchDirs:
    for SN_type in SN_types:
        catalog = 'astrocats'
        fn = 'The Open Supernova Catalog '+ SN_type+ '.csv'
        data_dir = os.path.join(OverArchDir,SN_type)
        Phot_lim = 20 #if doesn't have Phot lim photometry datapoints then we won't look at it
        if OverArchDir =='Data':
            Filters = ['g','r','i','z'] #if it doesnt have photometry in these filters we won't look at it
        else:
            Filters = ['B','V','R','I']
        try:
            os.makedirs(data_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            

        # In[57]:


        api_str_base = 'https://api.'+catalog+'.space/'
        df = pd.read_csv(fn)
        sort_df = df.sort_values(by=['Phot.'],axis=0,ascending=False)
        exp_df = sort_df.drop(['Radio','X-ray'],axis = 1, inplace=False)
        exp_df['Fit?'] = np.zeros(len(exp_df))
        exp_df = exp_df[exp_df['Phot.'] > Phot_lim]
        print(exp_df)
        #need to delete entries that do not have redshift

        # In[58]:


        #add column to dataframe containing absolute magnitudes - assumes apparent magnitudes are in AB system
        def appAB2Abs(df,DL):
            #DL in Mpc
            app = df['magnitude'].values
            Abs = app-5*np.log10(DL*10**6)+5
            df['AbsMag'] = pd.Series(Abs,index=df.index)
            return(df)


        # In[60]:


        arr_holder = np.zeros((1,len(exp_df.values[0])))
        goodSNs = 0
        badSNs = 0
        for i in range(len(exp_df)):
            rowDf = pd.DataFrame(exp_df.iloc[[i]])
            SN_name = rowDf['Name'].values[0]
            
            #get luminosity distance
            dLum = requests.get(api_str_base+SN_name+'/lumdist+redshift?complete&format=csv')
            api_response = dLum.text
            #check if data is available
            if api_response[2:9] == 'message':
                #remove from SNinfodf
                badSNs+=1
            else:
                dLum_df = pd.read_csv(StringIO(api_response))
                lumdist = dLum_df['lumdist'].values[0]
                redshift = dLum_df['redshift'].values[0]
                try:
                    redshift = redshift[0] #think the first is most believed
                except IndexError:
                    pass
                redshift = float(redshift)
                #make sure redshift makes sense
                if redshift >= 0:
                    rowDf['z'] = redshift
                    #get SN photometry
                    data = requests.get(api_str_base+SN_name+'/photometry/time+magnitude+e_magnitude+upperlimit+band+zeropoint+system+kcorrected+scorrected+mcorrected?format=csv')
                    api_response2 = data.text
                    if api_response2[2:9] == 'message':
                        #remove from SNinfodf
                        badSNs+=1

                    else:
                        SNdf = pd.read_csv(StringIO(api_response2))
                        #add Absolute Magnitudes
                        newDf = appAB2Abs(SNdf,lumdist)

                        #make sure the filters we want are in this file
                        bands = newDf['band'].values
                        lower_bands = []
                        for b in bands:
                            try:
                                lower_bands.append(b)
                            except AttributeError:
                                pass
                        FilterTroll = 0
                        for filt in Filters:
                            if filt in lower_bands or filt+"'" in lower_bands:
                                FilterTroll += 1
                        if FilterTroll != len(Filters): #not all the filters are in this file
                            badSNs+=1
                        else:
                            try:
                                invalidChars = set(string.punctuation.replace(":",""))
                                if any(char in invalidChars for char in SN_name):
                                    SN_name = SN_name.replace(":","-")
                                    exp_df.at[i,'Name'] = SN_name
                                os.makedirs(os.path.join(data_dir,SN_name)) #try to make a directory, will return OSerror if already exist
                            except OSError as e:
                                if e.errno != errno.EEXIST:
                                    raise
                            #we made it
                            rowDf.to_csv(os.path.join(data_dir,SN_name,'info.csv'))
                            newDf.to_csv(os.path.join(data_dir,SN_name,'Photometry.csv'))
                            arr_holder = np.concatenate((arr_holder,np.array([exp_df.iloc[i].values])),axis=0)
        arr_holder = arr_holder[1:,:]#get rid of first row
        master_holderDF = pd.DataFrame(arr_holder,columns = exp_df.columns.values)
        print(master_holderDF)
        master_holderDF.to_csv(os.path.join(data_dir,'OverallSNinfo.csv'))


        # In[55]:


        master_holderDF.to_csv(os.path.join(data_dir,'OverallSNinfo.csv'))


        # In[49]:

        '''

        SN_name = 'SN2015bn'
        dLum = requests.get(api_str_base+SN_name+'/lumdist+redshift?complete&format=csv')
        print(dLum.text)
        data = requests.get(api_str_base+SN_name+'/photometry/time+magnitude+e_magnitude+upperlimit+band+zeropoint+system+kcorrected+scorrected+mcorrected?format=csv')
        api_response2 = data.text
        SNdf = pd.read_csv(StringIO(api_response2))
        display(SNdf)
        '''

        # In[ ]:


        "'"

