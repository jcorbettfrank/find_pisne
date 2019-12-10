#This fits the observed light curves 



# coding: utf-8

# In[13]:


import requests
import os
import pandas as pd
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import errno
import pdb
from scipy.optimize import curve_fit
import pystan
import pickle


# In[14]:


#Seems like alpha is always -1, which is weird 

SNtype = 'IIn'
overArchDatDir = 'Data'
InfoCsv = 'OverallSNinfo.csv' #requires an info csv File to pull names from
modelName = 'StanModel'
zLim = .05
filters = ['g','r','i','z']
scale = 7
unc = .2 #did a very large uncertainty 
MaxNumIt = 20
seasLimBound = .3 #if normalized diff is below than will not be considered a new season
SeasonLimDays = 100 #there has to be at least one gap of SeasonLimDays in obs to break into seasons
ApriorLowerunc = .3 #exp(priorLowerunc) will be the prior uncertainties in alpha and Mp if the std is less than bc unc can't be neg
LpriorLowerunc = .3 #alpha and Mp


# In[15]:


overDir = os.path.join(overArchDatDir,SNtype)
OverSNinfoDF = pd.read_csv(os.path.join(overDir,InfoCsv))
Names = OverSNinfoDF['Name'].values
def getBadInds(filters,FiltID):
    #takes a list of filters as input
    all_inds = np.arange(len(FiltID),dtype=int)
    good_inds = np.array([],dtype=int)
    for filt in filters:
        filt_idxs = np.where(np.logical_or(FiltID == filt.upper(),FiltID == filt.lower()))[0].astype(int)
        good_inds = np.append(good_inds,filt_idxs)
    bad_inds = np.delete(all_inds,good_inds)
    return bad_inds

def Filters2Num(filters,FilterID):
    count = 1
    for filt in filters:
        idx = np.where(np.logical_or(FilterID == filt.lower(),FilterID==filt.upper()))
        FilterID[idx] = count
        count+=1
    return FilterID

def getIdxSep(filters,FilterID):
    #FilterID should be sorted first
    IdxSep = [[] for i in range(len(filters))]
    for i in range(len(filters)):
        filt = filters[i]
        filtIdxs = np.where(np.logical_or(FilterID == filt.lower(),FilterID==filt.upper()))
        lower = np.min(filtIdxs)
        upper = np.max(filtIdxs)
        upper+=1 #for slicing
        IdxSep[i] = [lower,upper]
    return IdxSep

def getZNames(Names,zlim,maxCount):
    #returns SN names & redshifs with redshift under zlim
    #max number of returned names = maxCount
    ZNames = []
    Zreds = []
    count = 0
    for n in Names:
        Datdir = os.path.join(overDir,n)
        info_fn = os.path.join(Datdir,'info.csv')
        SN_info = pd.read_csv(info_fn)
        redshifts = SN_info['z'].values[0]
        try:
            redshift = redshifts[0] #think the first is most believed
        except IndexError:
            redshift = redshifts
        if float(redshift) < zLim:
            ZNames.append(n)
            Zreds.append(float(redshift))
        count+=1
        if count >= maxCount:
            break
    return ZNames,Zreds


# In[16]:


#do it for 1 SN
numFilt = len(filters)
ZNames,Zreds = getZNames(Names,zLim,MaxNumIt)

#need some way of keeping track of names that were fitted (ideally good fits) so I can
# 1) not fit them again
# 2) load them for clustering later
ZNamesUnfitted = []
ZNamesUnfittedDrow = []
for ZName in ZNames:
    dfRowBool = OverSNinfoDF['Name'] == ZName
    dfRow = OverSNinfoDF[dfRowBool]
    dfRowIdx = np.where(dfRowBool==True)[0][0]
    #if dfRow['Fit?'].values[0] == 0:
    #hasn't been fit yet
    ZNamesUnfitted.append(ZName)
    ZNamesUnfittedDrow.append(dfRowIdx)

print(ZNamesUnfitted)


# In[ ]:


for SNidx,Name in enumerate(ZNamesUnfitted):
    redshift = Zreds[SNidx]
    Datdir = os.path.join(overDir,Name)
    info_fn = os.path.join(Datdir,'info.csv')
    phot_fn = os.path.join(Datdir,'Photometry.csv')
    preProcess_fn = os.path.join(Datdir,'preProcess.png')
    FitImg_fn = os.path.join(Datdir,'FitImg.png')
    phot_df = pd.read_csv(phot_fn)
    sortF_df = phot_df.sort_values(by=['band'])
    Mags = sortF_df['AbsMag'].values
    T = sortF_df['time'].values
    FilterID = sortF_df['band'].values
    
    #get rid of datapoints corresponding to a filter not in 'filters'
    bad_inds = getBadInds(filters,FilterID)

    FilterID = np.delete(FilterID,bad_inds)
    FilterIDN = np.copy(FilterID)
    FilterIDN = Filters2Num(filters,FilterIDN)
    FilterIDN=FilterIDN.astype(int)
    
    T = np.delete(T,bad_inds)
    Mags = np.delete(Mags,bad_inds)
    sig_Mags = unc*np.random.random_sample(size=np.shape(Mags))
    l = 10**((Mags/(-2.5))-scale)
    sig_l = 10**((sig_Mags/-2.5)-scale)
    IdxSep = getIdxSep(filters,FilterID)
    fig2 = plt.figure(figsize=(15,15))
    fig = plt.figure(figsize=(15,15))

    ax1 = fig.add_subplot(1,1,1)
    colorFilts = ['c','g','y','orange']
    ax1.set_xlabel('MJD',fontsize = 20)
    ax1.set_ylabel('Absolute Magnitude',fontsize = 20)
    ax1.set_title(Name,fontsize=20,pad = 20)
    ax1.yaxis.set_label_coords(-.08,.5)
    ax1.xaxis.set_label_coords(.5,-.08,)
    ax1.set_xticks([])
    ax1.set_yticks([])

    ax3 = fig2.add_subplot(1,1,1)
    colorFilts = ['c','g','y','orange']
    ax3.set_xlabel('MJD',fontsize = 20)
    ax3.set_ylabel('Absolute Magnitude',fontsize = 20)
    ax3.set_title(Name,fontsize=20,pad = 20)
    ax3.yaxis.set_label_coords(-.08,.5)
    ax3.xaxis.set_label_coords(.5,-.08,)
    ax3.set_xticks([])
    ax3.set_yticks([])

    count = 1
    alpha_guesses = []
    L_maxs = []
    T0s = []
    l_stan = np.array([])
    t_stan = np.array([])
    sig_l_stan = np.array([])
    FID_stan = np.array([],dtype=int)
    duringseasons=[]
    FilterTracker=np.array(list(range(numFilt)))
    for j in IdxSep:
        lower = j[0]
        upper = j[1]
        FID = FilterIDN[lower:upper]
        L = l[lower:upper]
        t = T[lower:upper]
        sig_L = sig_l[lower:upper]
        #now I could sort by date
        orderedI = t.argsort()
        t = t[orderedI]
        L = L[orderedI]
        FID = FID[orderedI]
        sig_L = sig_L[orderedI]
        
        ax4 = fig2.add_subplot(len(filters),1,count)
        ax4.set_xlim([np.min(T),np.max(T)])
        ax4.set_ylim([np.min(l),np.max(l)])
        ax4.set_title(filters[count-1])
        ax4.scatter(t,L, c=colorFilts[count-1])

        #break up into seasons
        diff = np.diff(t)
        second_diff = np.diff(diff)
        if np.max(diff) > SeasonLimDays:
            orderedD = second_diff.argsort()
            norm_second_diff = second_diff/np.max(second_diff)
            maxIdx = int(np.where(norm_second_diff==np.max(norm_second_diff))[0])

            seasIdx = np.where(norm_second_diff>seasLimBound)[0] #so last index will be before jump
            seasIdx +=2 #now it should not be inclusive, keep at +2
            #now we need to enter boundaries so we can index properly
            seasIdx = np.insert(seasIdx,len(seasIdx),len(t)-1)
            seasIdx = np.insert(seasIdx,0,0)
            #but there's a chance that the first or last index is already a bound
            #so take the unique vals just to make sure
            seasIdx = np.unique(seasIdx)
            #at this point we have len(seasIdx)-1 seasons, whos start indices go up to len(seasIdx-1)
        else:
            seasIdx = [0, len(t)-1]
            
        num_seasons = len(seasIdx)-1
        seasAvgs = np.zeros(num_seasons)
        for Sidx in range(num_seasons):
            LSeas = L[seasIdx[Sidx]:seasIdx[Sidx+1]]
            seasAvgs[Sidx] = np.mean(LSeas)
        ExplosionIdx = np.where(seasAvgs == np.max(seasAvgs))[0][0]
        #How it stands now, is this will cut out seasons that only have 1 data point
        '''
        print("seasIdx",seasIdx)
        print("seasAvgs",seasAvgs)
        print(ExplosionIdx)
        print(seasIdx,num_seasons)
        '''
        seasPlot = np.zeros(len(t))
        seasPlot[seasIdx] = np.max(l)/2
        ax4.plot(t,seasPlot)
        
        #Now have to say which seasons to use
        seas2use = [ExplosionIdx]

        if num_seasons == 2:
            if ExplosionIdx == 0:
                #first season holds explosion
                seas2use.append(ExplosionIdx+1)
            elif ExplosionIdx == num_seasons:
                seas2use.append(ExpSeas-1)
        elif num_seasons > 2:
            if ExplosionIdx == 0:
                #first season holds explosion
                seas2use.append(ExplosionIdx+1)
            elif ExplosionIdx == num_seasons:
                seas2use.append(ExplosionIdx-1)
            else:
                seas2use.append(ExplosionIdx+1) #plus two bc last index is upper bound python indexing
                seas2use.append(ExplosionIdx-1)
        #now have to add the end of the last season for python indexing
        if max(seas2use) != len(seasIdx)-1:
            seas2use.append(max(seas2use)+1)
        seasIdx = np.array(seasIdx) #for list indexing
        seas2use = np.array(seas2use)
        
        seasIdx2use = seasIdx[seas2use]
        #print(seasIdx,seas2use,seasIdx2use)
        lowerS = np.min(seasIdx2use)
        upperS = np.max(seasIdx2use)
    
    


        tCrop = t[lowerS:upperS]
        LCrop = L[lowerS:upperS]
        FIDCrop = FID[lowerS:upperS]
        sigLCrop = sig_L[lowerS:upperS]
        FilterTracker[count-1] =len(tCrop)
        l_stan = np.append(l_stan,LCrop)
        t_stan = np.append(t_stan,tCrop)
        FID_stan = np.append(FID_stan,FIDCrop)
        sig_l_stan = np.append(sig_l_stan,sigLCrop)
        ax2 = fig.add_subplot(len(filters),1,count)
        ax2.set_xlim([np.min(t),np.max(t)])
        ax2.set_ylim([np.min(L),np.max(L)])
        ax2.set_title(filters[count-1])
        ax2.scatter(tCrop,LCrop, c=colorFilts[count-1])

        #next thing is to determine if during season
        Lfirst = L[seasIdx[ExplosionIdx]]
        LExpSeas = L[seasIdx[ExplosionIdx]:seasIdx[ExplosionIdx+1]]
        if Lfirst == np.max(LExpSeas):
            #then could've exploded during gap and we wouldnt see
            duringseason = 0
            #use last time point of prev season to do alpha fit
            Tfit_st = seasIdx[ExplosionIdx]-1 #an index
        else:
            #there's some Lmax later in time
            duringseason = 1
            if ExplosionIdx != 0:
                Tfit_st = seasIdx[ExplosionIdx] #still the same thing..fit starting in explosion season depends on how "dim" the preceding points are
            else:
                Tfit_st = seasIdx[ExplosionIdx] #use the first day in explostion fit
        Tfit_st = seasIdx[ExplosionIdx] #always going to use explosion seas..if not duringseason taking the L
        duringseasons.append(duringseason)
        Mp = np.max(LExpSeas)
        Ip = int(np.where(L==Mp)[0][0])
        l_fit = L[Tfit_st:Ip+1]
        Tfit = t[Tfit_st:Ip+1]
        L_maxs.append(Mp)
        if len(Tfit) > 0:
            Tp = Tfit[-1]
            TFIT = Tfit-np.min(Tfit)
            Tp -= np.min(Tfit)
            #TFIT = Tfit
            ######### THINK T AND TP NEED TO SUBTRACT FROM ZERO
            def alphaRise(t,alpha):
                return Mp*(t/Tp)**np.exp(alpha)
            pVals,pcov = curve_fit(alphaRise,TFIT,l_fit,p0 = -1)
            ax2.plot(Tfit,alphaRise(TFIT,pVals[0]))
            alpha_guesses.append(pVals[0])
            T0s.append(Tfit[0])
            '''
            print(alpha_guesses)
            print(Tfit)
            print(t[seasIdx[ExplosionIdx]])
            print(t[Ip])
            '''
        count+=1
        plt.savefig(preProcess_fn,bbox_inces = "tight")
  

    N_obs = len(l_stan)
    L_max = np.mean(L_maxs)
    L_max_unc = np.std(L_maxs)
    alpha_avg = np.mean(alpha_guesses)
    alpha_avg_unc = np.std(alpha_guesses)
    #unc cant be negative, why would it be negative in the first place
    
    if alpha_avg_unc < np.exp(ApriorLowerunc):
        alpha_avg_unc = np.exp(ApriorLowerunc)
    if L_max_unc < np.exp(LpriorLowerunc):
        L_max_unc = np.exp(LpriorLowerunc)
    #alpha_avg = np.exp(-1)
    #alpha_avg_unc = np.exp(.5)
    if alpha_avg > 0:
        #we need a negative alpha
        alpha_avg = -1
    #alpha_avg_unc = .5 #TOTALLY IGNORING FIT, JUST BREAK UP INTO SEASONS IS WHAT'S USEFUL
    t0_mean = np.mean(T0s)
    fluxscale = 10**scale
    duringseason = np.any(duringseasons) #if during_season is diff for diff filters then will screw up
    if duringseason == False:
        duringseason = 0
        print(Name, 'Did Not Capture Rise.')
        OverSNinfoDF.at[ZNamesUnfittedDrow[SNidx],'Fit?'] = -1
    else:
        duringseason = 1

        Kcor_N = int(np.max(T)-np.min(T))
        Kcor = np.zeros((numFilt,Kcor_N))
        SN_dat = {'N_obs': N_obs,
                    'N_filt': numFilt,
                    't': t_stan,
                    'fL': l_stan,
                    'dfL': sig_l_stan,
                    'z': redshift,
                    't0_mean': t0_mean,
                    'J': FID_stan,
                    'Kcor_N': Kcor_N,
                    'Kcor': Kcor,
                    'duringseason': duringseason,
                    'fluxscale': fluxscale
                    }
        SN = """
        data {
        int<lower=0> N_obs;
        int<lower=0> N_filt;
        vector[N_obs] t;
        vector[N_obs] fL;
        vector[N_obs] dfL;
        real z;
        real t0_mean;
        int<lower=1,upper=N_filt> J[N_obs];
        int<lower=0> Kcor_N;
        real Kcor[N_filt,Kcor_N];
        real<lower=0> fluxscale;
        real<lower=0,upper=1> duringseason;
        }
        transformed data {
        vector[N_filt] prior_tp;
        vector[N_filt] prior_sig_tp;
        vector[N_filt] prior_lbeta1;
        vector[N_filt] prior_sig_lbeta1;
        vector[N_filt] prior_lbeta2;
        vector[N_filt] prior_sig_lbeta2;
        prior_tp[1] = log(20); 
        prior_tp[2] = log(20);
        prior_tp[3] = log(20);
        prior_tp[4] = log(20);
        //prior_tp[5] = log(30);
        prior_sig_tp[1] = 0.3;
        prior_sig_tp[2] = 0.3;
        prior_sig_tp[3] = 0.3;
        prior_sig_tp[4] = 0.3;
        //prior_sig_tp[5] = 0.3;
        //prior_lbeta1[1] = -2.1;
        //prior_lbeta1[2] = -2.3;
        prior_lbeta1[1] = -3.8;
        prior_lbeta1[2] = -3.8;
        prior_lbeta1[3] = -3.8;
        prior_lbeta1[4] = -3.8;
        //prior_lbeta1[5] = -4.0;
        //prior_sig_lbeta1[1] = 0.6;
        //prior_sig_lbeta1[2] = 0.8;
        prior_sig_lbeta1[1] = 1.5;
        prior_sig_lbeta1[2] = 1.5;
        prior_sig_lbeta1[3] = 1.5;
        prior_sig_lbeta1[4] = 1.5;
        //prior_sig_lbeta1[5] = 2;
        prior_lbeta2[1] = -3.4;
        prior_lbeta2[2] = -4.9;
        prior_lbeta2[3] = -4.9;
        prior_lbeta2[4] = -4.9;
        //prior_lbeta2[5] = -4.9;
        prior_sig_lbeta2[1] = 1.5;
        prior_sig_lbeta2[2] = 1.5;
        prior_sig_lbeta2[3] = 1.5;
        prior_sig_lbeta2[4] = 1.5;
        //prior_sig_lbeta2[5] = 1.5;
        }
        parameters {
        real pt0;
        vector<lower=0>[N_filt] t1;
        vector<lower=0>[N_filt] t2;
        vector<lower=0>[N_filt] td;
        vector<lower=0>[N_filt] tp;
        vector<upper=0>[N_filt] lalpha;
        vector<upper=0>[N_filt] lbeta1;
        vector<upper=0>[N_filt] lbeta2;
        vector<upper=0>[N_filt] lbetadN;
        vector<upper=0>[N_filt] lbetadC;
        vector<lower=0>[N_filt] Mp;
        vector[N_filt] Yb;
        vector<lower=0>[N_filt] V;
        }

        transformed parameters {
        vector[N_obs] mm;
        vector[N_obs] dm;
        vector<lower=0>[N_filt] M1;
        vector<lower=0>[N_filt] M2;
        vector<lower=0>[N_filt] Md;
        M1 = Mp ./ exp(exp(lbeta1) .* tp); 
        M2 = Mp .* exp(-exp(lbeta2) .* t2);
        Md = M2 .* exp(-exp(lbetadN) .* td);
        for (n in 1:N_obs) {
        real N_SNc;
        int Kc_up;
        int Kc_down;
        real t_exp;
        int j;
        int k;
        real mm_1;
        real mm_2;
        real mm_3;
        real mm_4;
        real mm_5;
        real mm_6;
        j = J[n];
        t_exp = (t[n] - (t0_mean + pt0)) / (1 + z);
        if (t_exp<0) {
        mm_1 = Yb[j];
        } else {
        mm_1 = 0; }
        if ((t_exp>=0) && (t_exp < t1[j])) {
        mm_2 = Yb[j] + M1[j] * pow(t_exp / t1[j] , exp(lalpha[j]));
        } else {
        mm_2 = 0;
        }
        if ((t_exp >= t1[j]) && (t_exp < t1[j] + tp[j])) {
        mm_3 = Yb[j] + M1[j] * exp(exp(lbeta1[j]) * (t_exp - t1[j]));
        } else {
        mm_3 = 0;
        }
        if ((t_exp >= t1[j] + tp[j]) && (t_exp < t1[j] + tp[j] + t2[j])) {
        mm_4 = Yb[j] + Mp[j] * exp(-exp(lbeta2[j]) * (t_exp - t1[j] - tp[j]));
        } else {
        mm_4 = 0;
        }
        if ((t_exp >= t1[j] + tp[j] + t2[j]) && (t_exp < t1[j] + tp[j] + t2[j] + td[j])) {
        mm_5 = Yb[j] + M2[j] * exp( - exp(lbetadN[j]) * (t_exp - t1[j] - tp[j] - t2[j]));
        } else {
        mm_5 = 0;
        }
        if (t_exp >= t1[j] + tp[j] + t2[j] + td[j]) {
        mm_6 = Yb[j] + Md[j] * exp(-exp(lbetadC[j]) * (t_exp - t1[j] - tp[j] - t2[j] - td[j]));
        } else {
        mm_6 = 0;
        }
        dm[n] = sqrt(pow(dfL[n],2) + pow(V[j],2));
        if (t_exp<0) {
        N_SNc = 0;
        } else if (t_exp < Kcor_N - 2){
        Kc_down = 0;
        while ((Kc_down+1) < t_exp) {
        Kc_down = Kc_down + 1;
        }
        Kc_up = Kc_down+1;
        N_SNc = Kcor[j,Kc_down+1] + (t_exp -
        floor(t_exp)) * (Kcor[j,Kc_up+1] - Kcor[j,Kc_down+1]);
        } else {
        N_SNc = Kcor[j,Kcor_N];
        }
        mm[n]=(mm_1+mm_2+mm_3+mm_4+mm_5+mm_6)
        / (pow(10, N_SNc/( -2.5)));
        }
        }

        model {
        if (duringseason == 1) {
        pt0 ~ skew_normal(-1, 1, -0.5);
        } else {
        pt0 ~ skew_normal(-30, 20, -1);
        }
        t1 ~ lognormal(log(1), 0.3);
        tp ~ lognormal(prior_tp, prior_sig_tp);
        t2 ~ lognormal(log(100), 0.3);
        td ~ lognormal(log(10), 0.5);
        lalpha ~ normal("""+str(alpha_avg)+""","""+str(alpha_avg_unc)+""");
        lbeta1 ~ normal(prior_lbeta1, prior_sig_lbeta1);
        lbeta2 ~ normal(prior_lbeta2, prior_sig_lbeta2);
        lbetadN ~ normal(-3, .5);
        lbetadC ~ normal(-5, 1);
        Mp ~ lognormal(log("""+str(L_max)+"""), log("""+str(L_max_unc)+"""));
        Yb ~ normal(0, 0.3);
        V ~ cauchy(0, 0.01);
        fL ~ normal(mm,dm);
        }
        """
        try:
            sm = pickle.load(open(os.path.join(Datdir,modelName)))
        except FileNotFoundError: 
            sm = pystan.StanModel(model_code = SN)
            with open(os.path.join(Datdir,modelName)+'.pkl','wb') as f:
                pickle.dump(sm,f) 

        fit = sm.sampling(data = SN_dat, iter = 250, chains = 4,n_jobs=1,algorithm = "NUTS")

        fit_data = fit.extract()
        fit_mags = fit_data['mm']
        fit_sig_mags = fit_data['dm']
        avgfit_mags = np.median(fit_mags,0)
        avgfit_sig_mags = np.median(fit_sig_mags,0)
        st_idx = 0
        count = 1
        fig = plt.figure(figsize=(15,15))
        ax1 = fig.add_subplot(1,1,1)
        colorFilts = ['c','g','y','orange']
        ax1.set_xlabel('MJD',fontsize = 20)
        ax1.set_ylabel('Absolute Magnitude',fontsize = 20)
        ax1.yaxis.set_label_coords(-.08,.5)
        ax1.xaxis.set_label_coords(.5,-.08,)
        ax1.set_title(Name,fontsize=20,pad = 20)
        ax1.set_xticks([])
        ax1.set_yticks([])
        lower = 0
        upper = 0
        for j in IdxSep:
            upper += FilterTracker[count-1]
            t = t_stan[lower:upper]
            L_exp = l_stan[lower:upper]
            L_exp_sig = sig_l_stan[lower:upper]

            L_fit = avgfit_mags[lower:upper]
            L_fit_sig = avgfit_sig_mags[lower:upper]
            lower=upper
            ax2 = fig.add_subplot(len(filters),1,count)
            ax2.set_xlim([np.min(t_stan),np.max(t_stan)])
            ax2.set_ylim([np.min(l_stan),np.max(l_stan)])
            ax2.set_title(filters[count-1])
            ax2.scatter(t,L_exp, c=colorFilts[count-1],label='Data')
            ax2.errorbar(t,L_fit,L_fit_sig,label='Fit')
            ax2.legend()
            count+=1
        plt.savefig(FitImg_fn,bbox_inces = "tight")

        outDf = pd.DataFrame() #creates a new dataframe that's empty
        #filter order follows order in filter list [g, r, i , z]
        for key,val in fit_data.items():
            dim = np.shape(val)
            if len(dim) == 1:
                outDf[key] = val
            elif dim[1] == numFilt:
                for j in range(numFilt):
                    filt = filters[j]
                    outDf[key+filt.upper()] = val[:,j]
        outDf.to_csv(os.path.join(Datdir,'StanFit.csv'))
        #Successful Fit (not necessarily a good fit)
        #now need to mark down in SNinfoDf and save
        #could eventually take standard deviation of fit into account to determine good fit
        OverSNinfoDF.at[ZNamesUnfittedDrow[SNidx],'Fit?'] = 1
        OverSNinfoDF.to_csv(os.path.join(overDir,InfoCsv))

