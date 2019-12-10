# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:33:31 2018

@author: Jordan
"""
#This just fits the theoretical light curves 


import pystan
import pickle
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[3]:
np.random.seed(7)
#input_files = ['He80','He100','He130']
#input_files = ['He80','He100','He130','R150','R175','R200','R225','R250','B200','B250']
#input_files = ['B200','B250']
#input_files = ['B250','R150','R175','R200']
#input_files = ['He80']
input_files = ['B200']
filters = ['g','r','i','z']
cadences = [1, 1, 1, 1]
MJD_zero = 2458205
scale = 7
numFilt = len(filters)

for fn in input_files:
    fn += '_BVRIMag'
    col_names = ['Days','Bmag','Vmag','Rmag','Imag']
    df = pd.read_csv(fn+'.txt', header = None,names = col_names,dtype = np.float64)
    L = np.array([])
    Sig_L = np.array([])
    FilterID = np.array([],dtype = 'int')
    T = np.array([])
    data = df.values
    data = data[20:,:] #cut off shock breakout
    epoch = data[:,0]
    MJD = epoch + MJD_zero
    Filter_Tracker = np.array(list(range(numFilt)))
    upper_unc = .05 
    lower_unc = .03 #in mags
    L_max = 0
    L_maxs = []
    T0s = []
    alpha_guesses = []
    for i in range(numFilt):
        mag = data[:,i+1]  
        sig_mag = (upper_unc-lower_unc)*np.random.random_sample(size=np.shape(mag))+lower_unc        
        l=10**((mag/(-2.5))-scale)
        sig_l = 10**((sig_mag/(-2.5))-scale)
        #now let's desample
        L_max += np.max(l)
        cad = cadences[i]
        if i < 3:
            good_inds = np.where(np.logical_or(np.logical_or(epoch % cad == 0,(epoch+.5)% cad == 0) ,(epoch-.5)% cad == 0))
        else:
            good_inds = np.where(epoch % cad <= 1)
        L_filt = l[good_inds]
        T_filt = MJD[good_inds]
        L = np.append(L,L_filt)
        Filter_Tracker[i] = np.shape(good_inds)[1]
        Sig_L = np.append(Sig_L,sig_l[good_inds])
        filterId = [int(i+1)for j in range(len(l[good_inds]))]
        FilterID = np.append(FilterID,filterId)
        T = np.append(T,T_filt)
        
        Tfit_st = 0
        Mp = np.max(L_filt)
        Ip = int(np.where(L_filt==Mp)[0][0])
        l_fit = L_filt[Tfit_st:Ip+1]
        Tfit = T_filt[Tfit_st:Ip+1]
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
            
            alpha_guesses.append(pVals[0])
            T0s.append(Tfit[0])
        
    
    L_max = np.max(L_maxs)
    L_max_unc = 3*np.std(L_maxs)
    FilterID = FilterID.astype(np.int)
    N_obs = len(L)
    z = .1
    t0_mean = MJD_zero +20
    fluxscale = 10**scale
    duringseason = 1
    Kcor_N = int(epoch[-1] - epoch[0])
    Kcor = np.zeros((numFilt,Kcor_N))
    SN_dat = {'N_obs': N_obs,
            'N_filt': numFilt,
            't': T,
            'fL': L,
            'dfL': Sig_L,
            'z': z,
            't0_mean': t0_mean,
            'J': FilterID,
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
    lalpha ~ normal(-1, .5);
    lbeta1 ~ normal(prior_lbeta1, prior_sig_lbeta1);
    lbeta2 ~ normal(prior_lbeta2, prior_sig_lbeta2);
    lbetadN ~ normal(-3, .5);
    lbetadC ~ normal(-5, 1);
    Mp ~ lognormal(log("""+str(L_max)+"""),log("""+str(L_max_unc)+"""));
    Yb ~ normal(0, 0.3);
    V ~ cauchy(0, 0.01);
    fL ~ normal(mm,dm);
    }
    """
    
    sm = pystan.StanModel(model_code = SN)
    with open('model4'+fn+'.pkl','wb') as f:
        pickle.dump(sm,f)
    
    fit = sm.sampling(data = SN_dat, iter = 500, chains = 4,n_jobs=1,algorithm = "NUTS")
    with open('fit_data'+fn+'.txt','w') as f:
        print(fit,file=f)
    
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
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title(fn)
    for i in Filter_Tracker:
        L_exp = L[st_idx:st_idx+i]
        Sig_L_exp = Sig_L[st_idx:st_idx+i]
        L_fit = avgfit_mags[st_idx:st_idx+i]
        Sig_L_fit = avgfit_sig_mags[st_idx:st_idx+i]
        #convert to magnitude
        M_exp = -2.5*np.log10((10**scale)*L_exp)
        Sig_M_exp = -2.5*np.log10((10**scale)*Sig_L_exp)
        M_fit = -2.5*np.log10((10**scale)*L_fit)
        Sig_M_fit = -2.5*np.log10((10**scale)*Sig_L_fit)
        date = T[st_idx:st_idx+i]
        st_idx = st_idx+i
        ax2 = fig.add_subplot(numFilt,1,count)
        ax2.errorbar(date,L_exp,yerr = Sig_L_exp, fmt = 'o',label = 'data',color = colorFilts[count-1])
        ax2.plot(date,L_fit,label = 'fit',color = colorFilts[count-1])
        ax2.set_title(filters[count-1],fontsize = 20)
        ax2.legend()
        #ax2.set_ylim(-18,-14)
        #ax2.invert_yaxis()
        count = count+1
    plt.savefig('stan_fit'+fn+'.png',bbox_inces = "tight")
    
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
    outDf.to_csv('StanFit'+fn+'.csv')
        