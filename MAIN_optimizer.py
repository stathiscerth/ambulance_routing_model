# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 13:21:48 2023

@author: user
"""
import survival_i_d as sid
import numpy as np
import pyomo
import pyomo.environ as pyo
import sys
import math
import time
import distance as pad
import os
# trying pyomo glpk
solvername = 'glpk'

# Get the directory of the current script (if the script is in the project root)
current_dir = os.path.dirname(__file__)  # __file__ is the current script's path

# Build the paths relative to the current script's directory
solverpath_folder = os.path.join(current_dir, 'glpkzip', 'glpk-4.65', 'w64')
solverpath_exe = os.path.join(solverpath_folder, 'glpsol.exe')

# does not need to be directly on c drive
sys.path.append(solverpath_folder)

#robust test... create a scenario based on the duration assumptions based schedule
#def scenario_schedule(m_c,h_c,a_c,initial_sch):
      
    
#    return final_scenario

def survival_score(sur_pro,amb_sche):
    s=0
    si=len(amb_sche)
    for i in range(si):
       # print(si)
        tc=amb_sche[i][0]
        tm_ar=amb_sche[i][3]
        s=s+sur_pro[tc][tm_ar]
    return s

#DATA, manually inserted!
#max evacuation time
max_time=240
#max time should be enough for all the patients to be transported. 
#it dependes on the number of ambulance, the distance between hospitals and mci, and the number of patients
#if not large enough, an Error like KeyError: "Index '(0, 11, 1)' is not valid for indexed component 'f'"  will appear to the screen

#survival probabilities
print('which of the 3 case 1-pessimistic, 2-moderate, 3-optimistic')

tmax = max_time+60
a_input = input()
case = int(a_input)
survival_probability_per_class = sid.survival_prob(tmax,case)
survival_prob_triage_halfhour = [survival_probability_per_class[1],survival_probability_per_class[2]]
#mci sites coordinates

mci1 = [40.63265, 22.94120]
mci2 = [40.59920, 22.94969]
mci3 = [40.56536, 22.98186]
mci_coo = [mci1, mci2, mci3]
mci_number = len(mci_coo)
#patients per triage at each mci, immediate and delayed at mci1, and immediate and delayed at mci2    
list_patient_triage_class = [[15,5],[12,12],[10,10]]    

#coordinates for hospitals of the scenario
h1=[40.57758,22.97037]
h2=[40.63050,23.04415]
h3=[40.63736,22.95888]
h4=[40.62965,22.96056]
h5=[40.63394,22.95572]
h6=[40.61283,22.96277]
h7=[40.67564,22.96160]
h_coo=[h1,h2,h3,h4,h5,h6,h7]
hospital_number = len(h_coo)
#distances from mcis
hospital_dist = [[0]*mci_number for _ in range(hospital_number)]
max_dist = 0
for mci in range(mci_number):
    for h in  range(hospital_number):
        ds=pad.distance_between(h_coo[h][1],h_coo[h][0],mci_coo[mci][1],mci_coo[mci][0])
        hospital_dist[h][mci] = ds
        if ds> max_dist:
            max_dist =ds
#hospital capacities  ....MUST be more in total than the total number of victims          
capacities = [7,16,12,12,6,14,18]
    
# initial coordinates of scenarios ambulance        
a1=[40.57758,22.97037]
a2=[40.63050,23.04415]
a3=[40.63736,22.95888]
a4=[40.62965,22.96056]
a5=[40.63394,22.95572]
a6=[40.61283,22.96277]
a7=[40.67564,22.96160]
a8=[40.57758,22.97037]
a9=[40.63050,23.04415]
a10=[40.63736,22.95888]
a11=[40.62965,22.96056]
a_coo=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11]        
ambulance_number = len(a_coo)
list_arrival_times = [[0]*mci_number for _ in range(ambulance_number)]
for i in range(0, ambulance_number):
    for mci in range(0, mci_number):
        ds=pad.distance_between(a_coo[i][1],a_coo[i][0],mci_coo[mci][1],mci_coo[mci][0])
        list_arrival_times[i][mci] =ds 
        
        
print('model STARTSSSS')
start1 = time.process_time()
start = time.process_time()
# your code here    


# model
model3 = pyo.ConcreteModel()
model3.S = pyo.RangeSet(0, 1)
model3.I = pyo.RangeSet(0, mci_number-1)
model3.Tb = pyo.RangeSet(-(max_dist+2), max_time)
model3.Ta = pyo.RangeSet(0, max_time+max_dist+1)
model3.To = pyo.RangeSet(0, max_time)
model3.H = pyo.RangeSet(0, hospital_number-1)
model3.A = pyo.RangeSet(0, ambulance_number-1)
model3.x = pyo.Var(model3.S, model3.I, model3.Tb, model3.H,
                   domain=pyo.NonNegativeIntegers, bounds=(0, ambulance_number))
model3.a = pyo.Var(model3.I, model3.Ta, model3.H,
                   domain=pyo.NonNegativeIntegers, bounds=(0, ambulance_number))
model3.f = pyo.Var(model3.I, model3.To, model3.A,
                   domain=pyo.NonNegativeIntegers, bounds=(0, 1))
survivals = sum(sum(sum(sum(model3.x[s, i, t, h]*survival_prob_triage_halfhour[s][t+hospital_dist[h][i]] for s in model3.S) for i in model3.I) for t in model3.To) for h in model3.H)
model3.OBJ = pyo.Objective(expr=survivals, sense=pyo.maximize)
model3.cons = pyo.ConstraintList()
# for j in model.K:
#  model.cons.add(sum(model.x[i,j] for i in model.I)<=pa_per_h)
# trick of negative values...all zero
for t in range(-(max_dist+2), 0):
    model3.cons.add(sum(sum(sum(
        model3.x[s, i, t, h] for s in model3.S)for i in model3.I)for h in model3.H) == 0)
# constraints  of patiens at the field
for s in model3.S:
    for i in model3.I:
        model3.cons.add(sum(sum(model3.x[s, i, t, h] for t in model3.To)
                       for h in model3.H) <= list_patient_triage_class[i][s])
# CAPACITY CONSTRAINTS
for h in model3.H:
        model3.cons.add(sum(sum(sum(model3.x[s, i, t, h] for t in model3.To)
                        for i in model3.I)for s in model3.S) <= capacities[h])        


#ambulance arrivals constraint
for a in range(0,ambulance_number):
    model3.cons.add(sum(model3.f[i,list_arrival_times[a][i],a] for i in model3.I) ==1)
for a in range(0,ambulance_number):
    for t in model3.To:
        for i in model3.I:
            if list_arrival_times[a][i]!=t:
                model3.cons.add(model3.f[i,t,a] ==0)
# ambulance availability constraints
model3.aa=pyo.Expression(model3.To,model3.I)
for i in model3.I:
    for t in model3.To:
        model3.aa[t,i]=sum(model3.a[i,t,h] for h in model3.H)+sum(model3.f[i,t,a] for a in model3.A)
for t in model3.To:
    for i in model3.I:
        model3.cons.add(sum(sum(model3.x[s, i, t, h] for s in model3.S)
                        for h in model3.H) <= model3.aa[t,i])
for t in model3.To:
    for h in model3.H:
        model3.cons.add(sum(model3.a[i,t+hospital_dist[h][i],h] for i in model3.I)<=sum(sum(model3.x[s, i, t-hospital_dist[h][i], h] for s in model3.S)for i in model3.I))
for t in range(max_time+1,max_time+1+max_dist+1):
    model3.cons.add(sum(sum(model3.a[ i, t, h] for h in model3.H)for i in model3.I)==0)  
for i in model3.I:
    for h in model3.H:
        for t in range(0,hospital_dist[h][i]):
            model3.cons.add(model3.a[i,t,h]==0)

list_for_plot_red=[[0,0]]   
list_for_plot_yellow=[[0,0]] 
# solve
opt = pyo.SolverFactory('glpk', executable=solverpath_exe)
# You can add a time limit by using the following command line
opt.options['tmlim'] = 120
solution = opt.solve(model3)
model3.OBJ.display()
REAL_OBJ=pyo.value(model3.OBJ)
schedule_process=[]
for t in model3.To:
    for i in model3.I:
        for h in model3.H:
            for s in model3.S:
                vl = pyo.value(model3.x[s,i,t,h])
                if (vl > 0):
                    print('minute '+str(t))
                    print('had '+str(vl)+' patients from place '+str(i)+' and class ' +
                          str(s)+' at hospital '+str(h))
                    if s==0:
                        new_ele=[t+hospital_dist[h][i],vl]
                        list_for_plot_red.append(new_ele)
                    elif s==1:
                        new_ele=[t+hospital_dist[h][i],vl]
                        list_for_plot_yellow.append(new_ele)
                    for jj in range(int(vl)):
                        schedule_process.append([s,i,t,h])
 
list_for_plot_red.sort()
list_for_plot_yellow.sort()                       
print(time.process_time() - start)  
# Code resumes execution after the delay

start = time.process_time()
mission=[]
amb_to_hosp=[]
first_a=list(range(ambulance_number))
for ii in schedule_process:
    fa=0
    for a in first_a:
       if pyo.value(model3.f[ii[1],ii[2],a])>0:
           mission.append(a)
           tt=ii[2]+hospital_dist[ii[3]][ii[1]]
           amb_to_hosp.append([a,ii[3],tt])
           first_a.remove(a)
           fa=1
           break
    if fa==0:
        mission.append(-1)
amb_to_hosp.sort(key=lambda x: x[0]) 

a_list=[[[0 for _ in range(0,hospital_number)] for _ in range(0,max_time+1 )] for _ in range(0,mci_number)]
for i in range(mci_number):
    for t in range(max_time+1):
        for h in range(hospital_number):
            a_list[i][t][h]=pyo.value(model3.a[i,t,h])  

for ii in range(len(mission)):
    jj=mission[ii]
    if jj==-1:
        current_mission=schedule_process[ii]
        place=current_mission[1]
        mtime=current_mission[2]
        destination=current_mission[3]
        for hh in model3.H:
            are=False
            if a_list[place][mtime][hh]>0:
                hosp_time=mtime-hospital_dist[hh][place]
                found_ambulance = -1
                are=True
                for row in amb_to_hosp:
                    if row[1] == hh and row[2] == hosp_time:
                        found_ambulance = row[0]
                        row[1]=destination
                        row[2]=mtime+hospital_dist[destination][place]
                        break
                if found_ambulance>-1:
                    mission[ii]=found_ambulance
                    a_list[place][mtime][hh]=a_list[place][mtime][hh]-1
                    break
                else:
                    print('mistake...oups times do not match')

            if hh==hospital_number-1 and not are:
                 print('mistake...oups... no available ambulance', ii)
            
for row, value in zip(schedule_process, mission):
    row.append(value)        


for am in model3.A:
    filtered_list = [row for row in schedule_process if row[-1] == am]
    print('ambulance '+str(am)+' will do these tasks:')    
    for row in filtered_list:
        row.insert(3,hospital_dist[row[3]][row[1]]+row[2])
        print('go to mci place '+str(row[1])+' at minute ' +str(row[2])+' and take patient of triage '+str(row[0])+' to hospital '+str(row[4]))
ambulance_schedule_opt=schedule_process.copy()
ambulance_schedule_opt.sort(key = lambda i: (i[5],i[2]))
score=survival_score(survival_prob_triage_halfhour,ambulance_schedule_opt)
print('optimal schedule expected_survivals:',REAL_OBJ)
#ambulance_schedule.... 0 class, 1 mci, 2 time of departure, 3 hospital, 4 ambulance... sorted by ambulance and by time
print("running time : ",time.process_time() - start1) 
print('MODEL END')



                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
print('model red first with weights')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
# model
model3 = pyo.ConcreteModel()
model3.S = pyo.RangeSet(0, 1)
model3.I = pyo.RangeSet(0, mci_number-1)
model3.Tb = pyo.RangeSet(-(max_dist+2), max_time)
model3.Ta = pyo.RangeSet(0, max_time+max_dist+1)
model3.To = pyo.RangeSet(0, max_time)
model3.H = pyo.RangeSet(0, hospital_number-1)
model3.A = pyo.RangeSet(0, ambulance_number-1)
model3.x = pyo.Var(model3.S, model3.I, model3.Tb, model3.H,
                   domain=pyo.NonNegativeIntegers, bounds=(0, ambulance_number))
model3.a = pyo.Var(model3.I, model3.Ta, model3.H,
                   domain=pyo.NonNegativeIntegers, bounds=(0, ambulance_number))
model3.f = pyo.Var(model3.I, model3.To, model3.A,
                   domain=pyo.NonNegativeIntegers, bounds=(0, 1))
survivals = sum(sum(sum(1000*(model3.x[0, i, t, h]*survival_prob_triage_halfhour[0][t+hospital_dist[h][i]])+0.01*(model3.x[1, i, t, h]*survival_prob_triage_halfhour[1][t+hospital_dist[h][i]])  for i in model3.I) for t in model3.To) for h in model3.H)
model3.OBJ = pyo.Objective(expr=survivals, sense=pyo.maximize)
model3.real_obj=pyo.Var(domain=pyo.NonNegativeReals)
model3.cons = pyo.ConstraintList()
model3.cons.add(model3.real_obj==sum(sum(sum(sum(model3.x[s, i, t, h]*survival_prob_triage_halfhour[s][t+hospital_dist[h][i]] for s in model3.S) for i in model3.I) for t in model3.To) for h in model3.H))
# for j in model.K:
#  model.cons.add(sum(model.x[i,j] for i in model.I)<=pa_per_h)
# trick of negative values...all zero
for t in range(-(max_dist+2), 0):
    model3.cons.add(sum(sum(sum(
        model3.x[s, i, t, h] for s in model3.S)for i in model3.I)for h in model3.H) == 0)
# constraints  of patiens at the field
for s in model3.S:
    for i in model3.I:
        model3.cons.add(sum(sum(model3.x[s, i, t, h] for t in model3.To)
                        for h in model3.H) <= list_patient_triage_class[i][s])
# CAPACITY CONSTRAINTS
#for h in model3.H:
    #for s in model3.S:
       # model3.cons.add(sum(sum(model3.x[s, i, t, h] for t in model3.To)
                      #  for i in model3.I) <= capacities[h][s])
        
for h in model3.H:
         model3.cons.add(sum(sum(sum(model3.x[s, i, t, h] for t in model3.To)
                         for i in model3.I)for s in model3.S) <= capacities[h])        

       
#ambulance arrivals constraint
for a in range(0,ambulance_number):
    model3.cons.add(sum(model3.f[i,list_arrival_times[a][i],a] for i in model3.I) ==1)
for a in range(0,ambulance_number):
    for t in model3.To:
        for i in model3.I:
            if list_arrival_times[a][i]!=t:
                model3.cons.add(model3.f[i,t,a] ==0)
# ambulance availability constraints
model3.aa=pyo.Expression(model3.To,model3.I)
for i in model3.I:
    for t in model3.To:
        model3.aa[t,i]=sum(model3.a[i,t,h] for h in model3.H)+sum(model3.f[i,t,a] for a in model3.A)
for t in model3.To:
    for i in model3.I:
        model3.cons.add(sum(sum(model3.x[s, i, t, h] for s in model3.S)
                        for h in model3.H) <= model3.aa[t,i])
for t in model3.To:
    for h in model3.H:
        model3.cons.add(sum(model3.a[i,t+hospital_dist[h][i],h] for i in model3.I)<=sum(sum(model3.x[s, i, t-hospital_dist[h][i], h] for s in model3.S)for i in model3.I))
for t in range(max_time+1,max_time+1+max_dist+1):
    model3.cons.add(sum(sum(model3.a[ i, t, h] for h in model3.H)for i in model3.I)==0)  
for i in model3.I:
    for h in model3.H:
        for t in range(0,hospital_dist[h][i]):
            model3.cons.add(model3.a[i,t,h]==0)
    
# solve
opt = pyo.SolverFactory('glpk', executable=solverpath_exe)
# You can add a time limit by using the following command line
opt.options['tmlim'] = 120
solution = opt.solve(model3)
schedule_process=[]
model3.OBJ.display()
list_for_plot_red_st=[[0,0]]   
list_for_plot_yellow_st=[[0,0]]
for t in model3.To:
    for i in model3.I:
        for h in model3.H:
            for s in model3.S:
                vl = pyo.value(model3.x[s,i,t,h])
                if (vl > 0):
                    print('minute '+str(t))
                    print('had '+str(vl)+' patients from place '+str(i)+' and class ' +
                          str(s)+' at hospital '+str(h))
                    if s==0:
                        new_ele=[t+hospital_dist[h][i],vl]
                        list_for_plot_red_st.append(new_ele)
                    elif s==1:
                        new_ele=[t+hospital_dist[h][i],vl]
                        list_for_plot_yellow_st.append(new_ele)
                    for jj in range(int(vl)):
                        schedule_process.append([s,i,t,h])    
                    
print('real objective value ',pyo.value(model3.real_obj))
REAL_OBJ_st=pyo.value(model3.real_obj)
print(time.process_time() - start)

mission=[]
amb_to_hosp=[]
first_a=list(range(ambulance_number))
for ii in schedule_process:
    fa=0
    for a in first_a:
       if pyo.value(model3.f[ii[1],ii[2],a])>0:
           mission.append(a)
           tt=ii[2]+hospital_dist[ii[3]][ii[1]]
           amb_to_hosp.append([a,ii[3],tt])
           first_a.remove(a)
           fa=1
           break
    if fa==0:
        mission.append(-1)
amb_to_hosp.sort(key=lambda x: x[0]) 

a_list=[[[0 for _ in range(0,hospital_number)] for _ in range(0,max_time+1 )] for _ in range(0,mci_number)]
for i in range(mci_number):
    for t in range(max_time+1):
        for h in range(hospital_number):
            a_list[i][t][h]=pyo.value(model3.a[i,t,h])  

for ii in range(len(mission)):
    jj=mission[ii]
    if jj==-1:
        current_mission=schedule_process[ii]
        place=current_mission[1]
        mtime=current_mission[2]
        destination=current_mission[3]
        for hh in model3.H:
            are=False
            if a_list[place][mtime][hh]>0:
                hosp_time=mtime-hospital_dist[hh][place]
                found_ambulance = -1
                are=True
                for row in amb_to_hosp:
                    if row[1] == hh and row[2] == hosp_time:
                        found_ambulance = row[0]
                        row[1]=destination
                        row[2]=mtime+hospital_dist[destination][place]
                        break
                if found_ambulance>-1:
                    mission[ii]=found_ambulance
                    a_list[place][mtime][hh]=a_list[place][mtime][hh]-1
                    break
                else:
                    print('mistake...oups times do not match')

            if hh==hospital_number-1 and not are:
                 print('mistake...oups... no available ambulance', ii)
            
for row, value in zip(schedule_process, mission):
    row.append(value)        

print("running time : ", time.process_time() - start) 

for am in model3.A:
    filtered_list = [row for row in schedule_process if row[-1] == am]
    print('ambulance '+str(am)+' will do these tasks:')    
    for row in filtered_list:
        row.insert(3,hospital_dist[row[3]][row[1]]+row[2])
        print('go to mci place '+str(row[1])+' at minute ' +str(row[2])+' and take patient of triage '+str(row[0])+' to hospital '+str(row[4]))
ambulance_schedule_opt_st=schedule_process.copy()
ambulance_schedule_opt_st.sort(key = lambda i: (i[5],i[2]))
score_st=survival_score(survival_prob_triage_halfhour,ambulance_schedule_opt_st)
print('immediate first, near hospital first. expected_survivals: ',REAL_OBJ_st)


