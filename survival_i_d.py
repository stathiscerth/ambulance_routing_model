# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 12:33:44 2023

@author: user"""

import numpy as np

def survival_prob(tmax,case):
    t=np.array(range(tmax))
    if case==1:
        bi0=0.09
        bi1=17
        bi2=1.01
        bd0=0.57
        bd1=61
        bd2=2.03
        surv_pr_i=bi0/(np.power((t/bi1),bi2)+1)
        surv_pr_d=bd0/(np.power((t/bd1),bd2)+1)
    elif case==2:
        bi0=0.24
        bi1=47
        bi2=1.3
        bd0=0.76
        bd1=138
        bd2=2.17
        surv_pr_i=bi0/(np.power((t/bi1),bi2)+1)
        surv_pr_d=bd0/(np.power((t/bd1),bd2)+1)
    elif case==3:
        bi0=0.56
        bi1=91
        bi2=1.58
        bd0=0.81
        bd1=160
        bd2=2.41
        surv_pr_i=bi0/(np.power((t/bi1),bi2)+1)
        surv_pr_d=bd0/(np.power((t/bd1),bd2)+1)
    tl=t.tolist()    
    surv_pr_il=surv_pr_i.tolist() 
    surv_pr_dl=surv_pr_d.tolist() 
    surv_pr_all=[tl,surv_pr_il,surv_pr_dl]
    return surv_pr_all

ax=survival_prob(600,3)            

