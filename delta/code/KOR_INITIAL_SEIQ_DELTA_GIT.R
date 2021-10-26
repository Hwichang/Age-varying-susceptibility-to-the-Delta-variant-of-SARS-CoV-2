########################################################################################################
################################ - Korea Initial SEIQ - ################################################
########################################################################################################
library(optiSolve)
library(quadprog)
library(polynom)
library(logitnorm)
library(rGammaGamma)
library(stats)
library(STAR)
rm(list=ls())
gc()

corona_daily = read.csv('korea_corona_daily_delta.csv')
corona_in = corona_daily[,11]
corona_in = as.numeric(gsub(',','',corona_in))
corona_daily = corona_daily[,1:10]

skage = read.csv('skage.csv',header=T)
skage = skage[1:101,5]
skage = as.numeric(sapply(skage, function(x) gsub(',','',x)))
age_limits_temp = c(10,20,30,40,50,60)

age_limits = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
skage.groups = rep(0,17)
skage.groups_temp = rep(0,6)
skage.groups_temp[1] = sum(skage[1:10])
for( i in 2:6){
  skage.groups_temp[i] = sum(skage[(age_limits_temp[(i-1)]+1):age_limits_temp[i]])
}

for( i in 1:length(corona_in)){
  corona_daily[i,2:7] = corona_daily[i,2:7] - corona_in[i]*skage.groups_temp/sum(skage.groups_temp)
}


skage.groups[1] = sum(skage[1:5])
for( i in 2:16){
  skage.groups[i] = sum(skage[(age_limits[(i-1)]+1):age_limits[i]])
}
skage.groups[17] = sum(skage[81:101])



age_limits_new = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
skage.groups_new = rep(0,16)
skage.groups_new[1] = sum(skage[1:5])
for( i in 2:15){
  skage.groups_new[i] = sum(skage[(age_limits[(i-1)]+1):age_limits[i]])
}
skage.groups_new[16] = sum(skage[76:101])


######################### - Load data - #########################
corona_new = read.csv('corona_daily_delta.csv')


######################### - INCUBATION PERIOD - #########################
Incu_param1 = 4.544
Incu_param2 = 1/0.709


######################### - TRANSMISSION ONSET - #########################
tran_dist_mu = -4
tran_param1 = 5.2662158
tran_param2 = 1/0.8709042 


######################### - Infection to Recover - #########################
I_R_param1 = 4
I_R_param2 = 4/5


######################### - Infection to Quarantine - #########################
C_param = 1.7


######################### - Symptom to Quarantine - #########################
symp_q_dist = read.csv('symp_q_dist.csv')$x


######################### - Set seed - #########################
set.seed(0814)


######################### - asymptomatic or symptomatic - #########################
total_I = sum(colSums(corona_new[,2:17]))
total_I_all = sum(colSums(corona_new[,2:17]))
index = total_I_all - total_I + 1
asym = rbinom(total_I_all,1,0.16)
#asym = rbinom(total_I_all,1,0.04) #asymomatic rate 0.04
#asym = rbinom(total_I_all,1,0.4) #asymomatic rate 0.4

n_sym = length(which(asym[index:total_I_all]==0))
n_asym = length(which(asym[index:total_I_all]==1))
seiq_matrix = data.frame('age'=NA,'Q_date'=NA,'I_date'=NA,'Y_date'=NA,'E_date'=NA)


######################### - initial sampling - #########################
sym_Y_list = rgamma(n_sym,Incu_param1,Incu_param2)
sym_I_list = tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)
sym_D_list = sample(symp_q_dist,n_sym,replace=TRUE)
asym_C_list = rexp(n_asym,C_param)
asym_L_list = rgamma(n_sym,Incu_param1,Incu_param2) + tran_dist_mu + rgamma(n_sym,tran_param1,tran_param2)


######################### - initial value of W - #########################
k=1
y=1
l=1

sum(corona_new[,2:17])

move = as.Date('2021-06-20')-as.Date('2021-01-01')

for( i in 1:nrow(corona_new)){
  for( j in 2:17){
    if(corona_new[i,j]!=0){
      for( h in 1:corona_new[i,j]){
        if( asym[k]==0){ #symptomatic
          seiq_matrix[k,1] = 5*(j-2) #age
          seiq_matrix[k,2] = i+move #Qurantined date
          seiq_matrix[k,4] = i+move - ceiling(sym_D_list[y]) #Symptom onset date
          seiq_matrix[k,5] = i+move - ceiling(sym_D_list[y] + sym_Y_list[y] ) #Exposed date
          seiq_matrix[k,3] = max(as.integer(i+move - ceiling(sym_D_list[y]) + ceiling(sym_I_list[y])) , seiq_matrix[k,5])  #transmission onset date
          k = k+1
          y = y+1
        }
        else{
          seiq_matrix[k,1] = 5*(j-2) #age
          seiq_matrix[k,2] = i+move #Qurantined date
          seiq_matrix[k,5] = i+move - ceiling(asym_C_list[l]) - max(ceiling(asym_L_list[l]),0) #Exposed date
          seiq_matrix[k,3] = i+move - ceiling(asym_C_list[l]) #transmission onset date
          seiq_matrix[k,4] = NA
          k = k+1
          l = l+1
        }
      }
    }
  }
}


