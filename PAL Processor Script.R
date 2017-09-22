########## Weston QC Processor ##########
#                                       #
#                                       #
#   Created by: Daniel Palmer, PhD      #
#                                       #
#                                       #
#########################################

## Establish Dependencies ##
library(reshape)
library(plyr)

############################

## Settings ##

##############

## Initial Variables ##
# ID #
validID = vector()

# Bad Names #
bad.list = list()
bad.list$bad.site1 = c('UWO', " UWO") #1
bad.list$bad.site2 = c('UOG', " UOG") #2

bad.list$bad.app = c('APPPS1', ' APPPS1') #3
bad.list$bad.5x = c('5FAD',' 5FAD') #4
bad.list$bad.3x = c('3XTG', ' 3XTG') #5

bad.wt = c(' w', 'w') #6
bad.tg = c(' t', 't') #7

bad.list$bad.app.wt = c() #8 
bad.list$bad.app.tg = c() #9
bad.list$bad.5x.wt = c() #10
bad.list$bad.5x.tg = c() #11
bad.list$bad.3x.wt = c() #12
bad.list$bad.3x.tg = c() #13

bad.list$bad.female = c() #14
bad.list$bad.male = c() #15

bad.list$bad.4m = c() #16
bad.list$bad.7m = c() #17
bad.list$bad.10m = c() #18

bad.list$bad.hab1 = c() #19
bad.list$bad.hab2 = c() #20
bad.list$bad.IT = c() #21
bad.list$bad.MT = c() #22
bad.list$bad.MI = c() #23
bad.list$bad.punish = c() #24

bad.list$bad.acq1 = c() #25
bad.list$bad.acq2 = c() #26
bad.list$bad.acq3 = c() #27
bad.list$bad.acq4 = c() #28

bad.list$bad.main1 = c() #29
bad.list$bad.main2 = c() #30
bad.list$bad.main3 = c() #31
bad.list$bad.main4 = c() #32

# Good Names #
good.names = vector()
good.names[1] = 'Site1'
good.names[2] = 'Site2'
good.names[3] = 'APP-PS1'
good.names[4] = '5xFAD'
good.names[5] = '3xTG-AD'
good.names[6] = NA
good.names[7] = NA
good.names[8] = 'C57BL6'
good.names[9] = 'APPPS1'
good.names[10] = 'B6SJLF1/J'
good.names[11] = '5xFAD'
good.names[12] = 'B6129SF2/J'
good.names[13] = '3xTG-AD'
good.names[14] = 'Female'
good.names[15] = 'Male'
good.names[16] = '4'
good.names[17] = '7'
good.names[18] = '10'
good.names[19] = NA
good.names[20] = NA
good.names[21] = "Phase I"
good.names[22] = 'Phase II'
good.names[23] = 'Phase III'
good.names[24] = 'Phase IV'
good.names[25] = 'dPAL Acquisition 1'
good.names[26] = 'dPAL Acquisition 2'
good.names[27] = 'dPAL Acquisition 3'
good.names[28] = 'sPAL Acquisition 1'
good.names[29] = 'dPAL1'
good.names[30] = 'dPAL2'
good.names[31] = 'dPAL3'
good.names[32] = 'dPAL4'

# Column Match to Name Check #
lookup.list = list()
lookup.list$TestSite = c(1:2)
lookup.list$Mouse.Strain = c(3:5)
lookup.list$Genotype = c(6:13)
lookup.list$Sex = c(14:15)
lookup.list$Age.Months = c(16:18)
lookup.list$Task = c(19:32)

#######################

## Function List ##
# Search for Col#, Bad Names, Good Names #
Name.Search = function(x, lookup.list=lookup.list, good.names=good.names, bad.names=bad.names){
  col.num = 0
  good.name = ''
  bad.names = c()
  for(a in 1:length(lookup.list)){
    if(isTRUE(x %in% unlist(lookup.list[a]))){
      
    }
  }
}

# QC Function PAL Pretrain #
QC.Pretrain.Function = function(dataset,validID){
  initial.ID.list = as.vector(unique(as.character(dataset[ ,1])))
}

# QC Function PAL Main #
QC.Main.Function

# Fix Naming Rules - Main & Acq #
Naming.MainAcq.Function = function(dataset, lookup.list=lookup.list, good.names=good.names, bad.names=bad.names){
  colnames(dataset)[1:11] = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day','Week')
  # Fix Spacing on relevant Columns #
  col.list.spacefix = c(1,2,3,4,5,6,7,10,11)
  for(a in col.list.spacefix){
    dataset[ ,a] = gsub(" ", "", dataset[ ,a])
  }
  # Fix Factor Names #
  if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) > 2)){
    geno.long = TRUE
    position.vec = c(1:5,8:13,14:32)
  }else if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) = 2)){
    geno.long = FALSE
    position.vec = c(1:5,6:7,14:32)
  }
  for(a in position.vec){
    if((a <= 1) & (a >= 2)){
      col.num = 3
    }else if((a <=3) & (a >= 5)){
      
    }else if((a <=6) & (a >= 13)){
      
    }else if((a <=14) & (a >= 15))
    
  }
  
}


###################

## Collect Initial Dataset ##
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Pretrain.csv')
raw.data.acquisition = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Acquisition.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Main.csv')
