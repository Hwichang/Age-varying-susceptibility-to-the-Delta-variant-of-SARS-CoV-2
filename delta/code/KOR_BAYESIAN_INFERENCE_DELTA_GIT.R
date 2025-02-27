########################################################################################################
################################ - Korea Bayesian Inference - #########################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(rGammaGamma)
library(stats)
library(STAR)
library(dplyr)
library(matrixStats)
library(lubridate)

rm(list=ls())
gc()
load("contact_school.RData") 
load("contact_work.RData") 
load("contact_others.RData") 
load("contact_all.RData") 
skage.groups_new = as.vector(read.csv('korea_population.csv',header=T)$x)
corona_new = read.csv('corona_daily_delta.csv')
skage = read.csv('skage.csv',header=T)
skage = skage[1:101,5]
skage = as.numeric(sapply(skage, function(x) gsub(',','',x)))
skage_0_2 = sum(skage[1:3])
skage_3_4 = sum(skage[4:5])
skage_5_6 = sum(skage[6:7])
skage_7_9 = sum(skage[8:10])
skage_10_12 = sum(skage[11:13])
skage_13_14 = sum(skage[14:15])
skage_15 = skage[16]
skage_16_18 = sum(skage[17:19])
skage_19 = skage[20]


###################School closure############################
school = as.matrix(read.csv('school.csv'))


###################Contact matrix#############################
contact_matrix = function(S){
  
  school_contact = as.matrix(contact_school$KOR)
  school_contact[1,] = (1- school[S,1]/100) * school_contact[1,]
  school_contact[2,] = (skage_5_6/(skage_5_6+skage_7_9))*(1-school[S,1]/100) * school_contact[2,] +  (skage_7_9/(skage_5_6+skage_7_9))*(1-school[S,2]/100) * school_contact[2,]
  school_contact[3,] = (skage_10_12/(skage_10_12+skage_13_14))*(1-school[S,2]/100) * school_contact[3,] +  (skage_13_14/(skage_10_12+skage_13_14))*(1-school[S,3]/100) * school_contact[3,]
  school_contact[4,] = (skage_15/(skage_15+skage_16_18+skage_19))*(1-school[S,3]/100) * school_contact[4,] +  (skage_16_18/(skage_15+skage_16_18+skage_19))*(1-school[S,4]/100) * school_contact[4,]
  
  res = as.matrix(contact_all$KOR) - school_contact - (0.07)*as.matrix(contact_work$KOR) - (0.02)*as.matrix(contact_others$KOR) 
  
  return(res)
}


###################Vaccine calender #######################
AZ_1_res = as.matrix(read.csv('AZ_1.csv'))
AZ_2_res = as.matrix(read.csv('AZ_2.csv'))
PF_1_res = as.matrix(read.csv('PF_1.csv'))
PF_2_res = as.matrix(read.csv('PF_2.csv'))
M_1_res = as.matrix(read.csv('M_1.csv'))
M_2_res = as.matrix(read.csv('M_2.csv'))
AZ_PF_res = as.matrix(read.csv('AZ_PF.csv'))
JJ_res = as.matrix(read.csv('JJ.csv'))


##########################efficacy_lower###################
# AZ_1_eff = 0.35
# AZ_2_eff = 0.62
# PF_1_eff = 0.50
# PF_2_eff = 0.77
# M_1_eff = 0.64
# M_2_eff = 0.84
# JJ_eff = 0.67


############################efficacy########################
AZ_1_eff = 0.46
AZ_2_eff = 0.67
PF_1_eff = 0.57
PF_2_eff = 0.8
M_1_eff = 0.75
M_2_eff = 0.8
JJ_eff = 0.69


##########################efficacy_upper#####################
# AZ_1_eff = 0.55
# AZ_2_eff = 0.71
# PF_1_eff = 0.63
# PF_2_eff = 0.83
# M_1_eff = 0.83
# M_2_eff = 0.89
# JJ_eff = 0.71


#######################InCUBATION PERIOD######################
Incu_param1 = 4.544
Incu_param2 = 1/0.709


######################TRANSMISSION ONSET######################
tran_dist_mu = -4
tran_param1 = 5.2662158
tran_param2 = 1/0.8709042 


##################Infection to Recover########################
I_R_param1 = 4
I_R_param2 = 4/5


##################Infection to Quarantine######################
C_param = 1.7


#####################Symptom to Quarantine#####################
symp_q_dist = read.csv('symp_q_dist.csv')$x


#####################Set seed#################################
set.seed(0814)


####################initial value of W#######################
#seiq_matrix = read.csv('initial_seiq_korea_delta_asym_0.04.csv')
seiq_matrix = read.csv('initial_seiq_korea_delta.csv')
#seiq_matrix = read.csv('initial_seiq_korea_delta_asym_0.4.csv')


##################initial susceptible########################
confirmed_temp = c(6763,10900,22836,20984,23077,27671,22622,10175,5851)

skage.groups_temp = skage.groups_new
skage.groups_temp[16] = 1587676
confirmed = rep(0,16)
for(i in 1:8){
  confirmed[2*i-1] = round(confirmed_temp[i]*skage.groups_temp[2*i-1]/(skage.groups_temp[2*i-1]+skage.groups_temp[2*i]))
  confirmed[2*i] = round(confirmed_temp[i]*skage.groups_temp[2*i]/(skage.groups_temp[2*i-1]+skage.groups_temp[2*i]))
}
confirmed[16] = confirmed[16]+confirmed_temp[9]


vaccine_sus = matrix(0, nrow=300,ncol=16)
for(t in 30:300){
  vaccine_sus[t,] = round(AZ_1_eff*colSums(AZ_1_res[(1:(t-21)),]) + (AZ_2_eff-AZ_1_eff)*colSums(AZ_2_res[(1:(t-14)),])+
                            PF_1_eff*colSums(PF_1_res[(1:(t-21)),]) + (PF_2_eff-PF_1_eff)*colSums(PF_2_res[(1:(t-14)),])+
                            M_1_eff*colSums(M_1_res[(1:(t-21)),]) + (M_2_eff-M_1_eff)*colSums(M_2_res[(1:(t-14)),])+
                            JJ_eff*colSums(JJ_res[(1:(t-14)),]) + (PF_2_eff- AZ_1_eff)*colSums(AZ_PF_res[(1:(t-14)),]))
}

suscept = matrix(0,nrow = 300 ,ncol=16)

for(t in 1:nrow(suscept)){
  suscept[t,] = sapply(1:16,function(x) {skage.groups_new[x] - vaccine_sus[t,x]- nujuck[x] - length(which((seiq_matrix$E_date<(t))&(seiq_matrix$age/5+1==x)))})
}

suscept_new = suscept


###############initial Exposed###############################
exposed = matrix(0,nrow = 300,ncol=16)
for(t in 1:nrow(exposed)){
  exposed[t,] = sapply(1:16,function(x) {length(which((seiq_matrix$E_date==(t))&(seiq_matrix$age/5+1==x)))})
}


##############initial Infectious#############################

I_total_sym = matrix(0,nrow = 300,ncol=16)
for(t in 1:nrow(I_total_sym)){
  I_total_sym[t,] = sapply(1:16,function(x){length(which((seiq_matrix$I_date<=(t))&(seiq_matrix$Q_date>(t))&(seiq_matrix$age/5+1==x)&(!is.na(seiq_matrix$Y_date))))})
}

I_total_asym = matrix(0,nrow = 300,ncol=16)
for(t in 1:nrow(I_total_asym)){
  I_total_asym[t,] = sapply(1:16,function(x){length(which((seiq_matrix$I_date<=(t))&(seiq_matrix$Q_date>(t))&(seiq_matrix$age/5+1==x)&(is.na(seiq_matrix$Y_date))))})
}

I_total = I_total_sym + 0.5*I_total_asym
I_total_new = I_total


####################### q_matrix ###############################

# naive estimate
naive_p = function(t){
  res = (1/(skage.groups_new)) * (contact_matrix(t-1)%*%I_total[(t-1),])
  return(res)
}

init_q = matrix(0, nrow=6, ncol=16)
naive_q = (exposed[2:nrow(exposed),]/suscept[1:(nrow(exposed)-1),])/t(sapply(2:nrow(exposed), function(x) {naive_p(x)}))
naive_q = replace(naive_q , is.infinite(naive_q ),NaN)

q_standard = c(1,177,233,253,300)
for( i in 1:4){
  temp = naive_q[q_standard[i]:(q_standard[i+1]-1),]
  temp = na.omit(temp)
  init_q[i,] = colMeans(temp)
}


init_q[2,] = sapply(1:16, function(x) { max(rgamma(1,0.001,rate = 0.001),1e-300)})  

q_matrix = matrix(0, nrow = 299,ncol = 16)
q_matrix[1,] = init_q[1,]
for( i in 1:4){
  for(j in q_standard[i]:(q_standard[i+1]-1)){
    q_matrix[j,] = init_q[i,]
  }
}

############################ p(t)#############################
## p_unit : force of infection at time t
## p_unit_new : force of infection at time t divided by q
## tilde_p_unit : force of infection at time t about newly imputed value
## tilde_p_unit_new : force of infection at time t divided by q newly imputed value

p_unit = function(t){
  res = (q_matrix[(t-1),]/(skage.groups_new)) * (contact_matrix(t-1)%*%I_total[t-1,])
  return(res)
}

p_unit_new = function(t){
  res = (1/(skage.groups_new)) * (contact_matrix(t-1)%*%I_total[t-1,])
  return(res)
}

tilde_p_unit = function(t){
  res = (q_matrix[(t-1),]/(skage.groups_new)) * (contact_matrix(t-1)%*%I_total_new[t-1,])
  return(res)
}

tilde_p_unit_new = function(t){
  res = (1/(skage.groups_new)) * (contact_matrix(t-1)%*%I_total_new[t-1,])
  return(res)
}


################ For caculate P(E=t) ##########################
p = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {p_unit(x)})))
  return(res)
}

p_new = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {p_unit_new(x)})))
  return(res)
}

tilde_p = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {tilde_p_unit(x)})))
  return(res)
}

tilde_p_new = function(s,t){
  res = -colSums(t(sapply(s:t,function(x) {tilde_p_unit_new(x)})))
  return(res)
}



start_date = 177 # 2021-06-27
end_date = 233 # 2020-08-21

start = which(seiq_matrix$Q_date>=start_date)[1]
end = which(seiq_matrix$Q_date<=end_date + 14)[length(which(seiq_matrix$Q_date<=end_date + 14))]
end-start
q_list = matrix(0, nrow=1100,ncol=16)

k=1

for(iter in 1:1100){
  start_time <- Sys.time()
  
  ######################################### - W update - ##########################################################
  for( i in start:end){
    
    Q_i = seiq_matrix$Q_date[i]
    i_age = seiq_matrix$age[i]/5+1
    E_old = seiq_matrix$E_date[i]
    I_old = seiq_matrix$I_date[i]
    alpha = 0
    
    if(!is.na(seiq_matrix$Y_date[i])){
      c = 1
      Y_new = Q_i -  sample(symp_q_dist,1)
      E_new = Y_new -  ceiling(rgamma(1,Incu_param1,Incu_param2))
      I_new = max(Y_new + min(ceiling(tran_dist_mu + rgamma(1,tran_param1,tran_param2)),14),E_new)

    }else{
      c = 0.5
      Y_new = NA
      I_new = Q_i - ceiling(rexp(1,C_param))
      E_new = I_new - max(ceiling(rgamma(1,Incu_param1,Incu_param2) + tran_dist_mu + rgamma(1,tran_param1,tran_param2)),0)
    }
    
    exposed[E_old,i_age] = exposed[E_old,i_age] - 1 #substract i
    
    m_E = min(E_old,E_new)
    M_E = max(E_old,E_new)
    
    m_I = min(I_old,I_new)
    M_I = max(I_old,I_new)
    
    
    if(E_old > E_new){
      suscept_new[E_new:(E_old-1),i_age] = sapply(E_new:(E_old-1), function(x){suscept_new[x,i_age]-1})
    }else if(E_old < E_new){
      suscept_new[E_old:(E_new-1),i_age] = sapply(E_old:(E_new-1), function(x){suscept_new[x,i_age]+1})
    }
    
    if((I_old >= Q_i) & (I_new < Q_i)){
      I_total_new[I_new:(Q_i-1),i_age] = sapply(I_new:(Q_i-1), function(x){I_total_new[x,i_age]+c})
      
    }else if((I_old > I_new) & (I_old < Q_i)){
      I_total_new[I_new:(I_old-1),i_age] = sapply(I_new:(I_old-1), function(x){I_total_new[x,i_age]+c})
      
    }else if((I_old < I_new) & (I_new < Q_i)){
      I_total_new[I_old:(I_new-1),i_age] = sapply(I_old:(I_new-1), function(x){I_total_new[x,i_age]-c})
      
    }else if((I_new >= Q_i) & (I_old < Q_i)){
      I_total_new[I_old:(Q_i-1),i_age] = sapply(I_old:(Q_i-1), function(x){I_total_new[x,i_age]-c})
    }
    
    ##################### W_(-i) MCMC ratio ############################
    age_num_B = suscept_new[M_I,]
    temp_B = (tilde_p(m_I+1,M_I) - p(m_I+1,M_I))*age_num_B
    alpha = alpha + sum(temp_B)
    
    age_num_A = exposed[(m_I+1):M_I,]
    temp_A = t(sapply((m_I+1):M_I , function(y) { log(tilde_p_unit(y)) + tilde_p(m_I,y-1) - log(p_unit(y)) - p(m_I,y-1)}))
    alpha = alpha + sum(temp_A*age_num_A)
    
    
    ######################### W_i MCMC ratio #########################################################################
    alpha = alpha + (log(tilde_p_unit(E_new)) + tilde_p(m_E-1,E_new-1) - log(p_unit(E_old)) - p(m_E-1,E_old-1))[i_age]
    
    
    ########################### accept or reject ###################################################################
    accept = runif(1,0,1)
    if( log(accept)  <= alpha){
      seiq_matrix$I_date[i] = I_new
      seiq_matrix$Y_date[i] = Y_new
      seiq_matrix$E_date[i] = E_new
      I_total = I_total_new
      suscept = suscept_new
      exposed[E_new,i_age] = exposed[E_new,i_age]+1
    }else{
      I_total_new = I_total
      suscept_new = suscept
      exposed[E_old,i_age] = exposed[E_old,i_age]+1
    }
    
  }
  
  
  ########################################### - q update - ##########################################################
  data_num = exposed[start_date:end_date,]
  
  data_vec = colSums(t(sapply(start_date:end_date, function(x)  {-p_new(start_date,x)}))*data_num)
  data_vec = data_vec - p_new(start_date,end_date)* suscept[end_date,]
  
  q_new = sapply(1:16, function(x) rgamma(1 ,  0.001 + colSums(data_num)[x] , rate = 0.001 + data_vec[x]))
  q_list[k,] = q_new
  
  q_matrix[start_date:end_date,] = t(sapply(start_date:end_date , function(x) {q_new}))
  
  print(q_matrix[start_date,])
  
  k = k+1 
  
  end_time <- Sys.time()
  print(end_time - start_time)
}

