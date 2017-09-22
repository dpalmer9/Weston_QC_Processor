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

bad.list$bad.wt = c(' w', 'w') #6
bad.list$bad.tg = c(' t', 't') #7

bad.list$bad.app.wt = c(' w', 'w') #8 
bad.list$bad.app.tg = c(' t', 't') #9
bad.list$bad.5x.wt = c(' w', 'w') #10
bad.list$bad.5x.tg = c(' t', 't') #11
bad.list$bad.3x.wt = c(' w', 'w') #12
bad.list$bad.3x.tg = c(' t', 't') #13

bad.list$bad.female = c('F') #14
bad.list$bad.male = c('M') #15

bad.list$bad.4m = c(4,'10') #16
bad.list$bad.7m = c(7,'7') #17
bad.list$bad.10m = c(10,'4') #18

bad.list$bad.hab1 = c('Hab1') #19
bad.list$bad.hab2 = c('Hab2') #20
bad.list$bad.IT = c('MousePALInitialTouchTrainingv3') #21
bad.list$bad.MT = c('MousePALMustTouchTrainingv3') #22
bad.list$bad.MI = c('MousePALMustInitiateTrainingv3') #23
bad.list$bad.punish = c('MousePALPunishIncorrectTrainingv3') #24

bad.list$bad.acq1 = c('MousedPALacquisition1v3') #25
bad.list$bad.acq2 = c('MousedPALacquisition2v3') #26
bad.list$bad.acq3 = c('MousedPALacquisition6v3','MousedPALacquisitionSIMPLEv3') #27
bad.list$bad.acq4 = c('MousesPALacquisition1v3') #28

bad.list$bad.main1 = c('MousedPAL1v3') #29
bad.list$bad.main2 = c('MousedPAL2v3') #30
bad.list$bad.main3 = c('MousedPALSIMPLEv3') #31
bad.list$bad.main4 = c('MousesPAL1v3') #32

# Good Names #
good.list = vector()
good.list[1] = 'Site1'
good.list[2] = 'Site2'
good.list[3] = 'APP-PS1'
good.list[4] = '5xFAD'
good.list[5] = '3xTG-AD'
good.list[6] = 'W'
good.list[7] = 'T'
good.list[8] = 'C57BL6'
good.list[9] = 'APPPS1'
good.list[10] = 'B6SJLF1/J'
good.list[11] = '5xFAD'
good.list[12] = 'B6129SF2/J'
good.list[13] = '3xTG-AD'
good.list[14] = 'Female'
good.list[15] = 'Male'
good.list[16] = '4'
good.list[17] = '7'
good.list[18] = '10'
good.list[19] = NA
good.list[20] = NA
good.list[21] = "Phase I"
good.list[22] = 'Phase II'
good.list[23] = 'Phase III'
good.list[24] = 'Phase IV'
good.list[25] = 'dPAL Acquisition 1'
good.list[26] = 'dPAL Acquisition 2'
good.list[27] = 'dPAL Acquisition 3'
good.list[28] = 'sPAL Acquisition 1'
good.list[29] = 'dPAL1'
good.list[30] = 'dPAL2'
good.list[31] = 'dPAL3'
good.list[32] = 'dPAL4'

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

# QC Function PAL Pretrain #
QC.Pretrain.Function = function(dataset){
  new.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
  colnames(new.data) = colnames(dataset)
  initial.ID.list = as.vector(unique(as.character(dataset[ ,1])))
  for(a in 1:length(initial.ID.list)){
    temp.id = initial.ID.list[a]
    temp.data = dataset[which(dataset[ ,1] == temp.id), ]
    temp.data = temp.data[order(dataset$Date), ]
    temp.data.dateunique = as.vector(unique(temp.data[ ,8]))
    if(length(temp.data.dateunique) != length(temp.data[ ,8])){
      temp.data.2 = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
      colnames(temp.data.2) = temp.data
      for(b in 1:length(temp.data.dateunique)){
        data.check = temp.data[which(temp.data[ ,8] == temp.data.dateunique[b]), ]
        if(nrow(data.check) > 1){
          data.check = data.check[1, ]
          temp.data.2 = rbind(temp.data.2,data.check)
        }else{
          temp.data.2 = rbind(temp.data.2, data.check)
        }
      }
      temp.data = temp.data.2
    }
    new.data = rbind(new.data, temp.data)
  }
  return(new.data)
}

# QC Function PAL Acq #
QC.Pretrain.Function = function(dataset){
  new.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
  colnames(new.data) = colnames(dataset)
  initial.ID.list = as.vector(unique(as.character(dataset[ ,2])))
  for(a in 1:length(initial.ID.list)){
    temp.id = initial.ID.list[a]
    temp.data = dataset[which(dataset[ ,2] == temp.id), ]
    temp.data = temp.data[order(dataset$Date), ]
    temp.data.dateunique = as.vector(unique(temp.data[ ,9]))
    if(length(temp.data.dateunique) != length(temp.data[ ,9])){
      temp.data.2 = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
      colnames(temp.data.2) = temp.data
      for(b in 1:length(temp.data.dateunique)){
        data.check = temp.data[which(temp.data[ ,9] == temp.data.dateunique[b]), ]
        if(nrow(data.check) > 1){
          data.check = data.check[1, ]
          temp.data.2 = rbind(temp.data.2,data.check)
        }else{
          temp.data.2 = rbind(temp.data.2, data.check)
        }
      }
      temp.data = temp.data.2
    }
    new.data = rbind(new.data, temp.data)
  }
  return(new.data)
}
# QC Function PAL Main #
#QC.Main.Function

# Fix Naming Rules - Main & Acq #
Naming.MainAcq.Function = function(dataset, good.names=good.list, bad.names=bad.list){
  colnames(dataset)[1:11] = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day','Week')
  col.list.spacefix = c(1,2,3,4,5,6,7,8,10,11)
  for(a in col.list.spacefix){
    dataset[ ,a] = gsub(" ", "", dataset[ ,a])
  }
  if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) > 2)){
    geno.long = TRUE
    position.vec = c(1:5,8:13,14:32)
  }else if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) == 2)){
    geno.long = FALSE
    position.vec = c(1:5,6:7,14:32)
  }
  for(a in position.vec){
    if((a >= 1) & (a <= 2)){
      col.num = 3
    }else if((a >=3) & (a <= 5)){
      col.num = 4
    }else if((a >=6) & (a <= 13)){
      col.num = 5
    }else if((a >=14) & (a <= 15)){
      col.num = 6
    }else if((a >=16) & (a <= 18)){
      col.num = 7
    }else if((a >= 19)){
      col.num = 8
    }
    temp.badnames = as.vector(unlist(bad.names[a]))
    temp.goodname = good.names[a]
    for(b in 1:nrow(dataset)){
      if(isTRUE(is.element(dataset[b,col.num],temp.badnames))){
        dataset[b,col.num] = temp.goodname
      }
      if(isTRUE(dataset[b,col.num] == 'W')){
        if(isTRUE(dataset[b,4] == "APP-PS1")){
          dataset[b,col.num] = good.list[8]
        }else if(isTRUE(dataset[b,4] == "5xFAD")){
          dataset[b,col.num] = good.list[10]
        }else if(isTRUE(dataset[b,4] == "3xTG-AD")){
          dataset[b,col.num] = good.list[12]
        }
      }
      if(isTRUE(dataset[b,col.num] == 'T')){
        if(isTRUE(dataset[b,4] == "APP-PS1")){
          dataset[b,col.num] = good.list[9]
        }else if(isTRUE(dataset[b,4] == "5xFAD")){
          dataset[b,col.num] = good.list[11]
        }else if(isTRUE(dataset[b,4] == "3xTG-AD")){
          dataset[b,col.num] = good.list[13]
        }
      }
    }
  }
  return(dataset)
}

# Fix Naming Rules - Pretrain #
Naming.Pretrain.Function = function(dataset, good.names=good.list, bad.names=bad.list){
  colnames(dataset)[1:10] = c('AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day','Week')
  col.list.spacefix = c(1,2,3,4,5,6,7,9,10)
  for(a in col.list.spacefix){
    dataset[ ,a] = gsub(" ", "", dataset[ ,a])
  }
  if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) > 2)){
    geno.long = TRUE
    position.vec = c(1:5,8:13,14:32)
  }else if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) == 2)){
    geno.long = FALSE
    position.vec = c(1:5,6:7,14:32)
  }
  for(a in position.vec){
    if((a >= 1) & (a <= 2)){
      col.num = 2
    }else if((a >=3) & (a <= 5)){
      col.num = 3
    }else if((a >=6) & (a <= 13)){
      col.num = 4
    }else if((a >=14) & (a <= 15)){
      col.num = 5
    }else if((a >=16) & (a <= 18)){
      col.num = 6
    }else if((a >= 19)){
      col.num = 7
    }
    temp.badnames = as.vector(unlist(bad.names[a]))
    temp.goodname = good.names[a]
    for(b in 1:nrow(dataset)){
      if(isTRUE(is.element(dataset[b,col.num],temp.badnames))){
        dataset[b,col.num] = temp.goodname
      }
      if(isTRUE(dataset[b,col.num] == 'W')){
        if(isTRUE(dataset[b,4] == "APP-PS1")){
          dataset[b,col.num] = good.list[8]
        }else if(isTRUE(dataset[b,4] == "5xFAD")){
          dataset[b,col.num] = good.list[10]
        }else if(isTRUE(dataset[b,4] == "3xTG-AD")){
          dataset[b,col.num] = good.list[12]
        }
      }
      if(isTRUE(dataset[b,col.num] == 'T')){
        if(isTRUE(dataset[b,4] == "APP-PS1")){
          dataset[b,col.num] = good.list[9]
        }else if(isTRUE(dataset[b,4] == "5xFAD")){
          dataset[b,col.num] = good.list[11]
        }else if(isTRUE(dataset[b,4] == "3xTG-AD")){
          dataset[b,col.num] = good.list[13]
        }
      }
    }
  }
  return(dataset)
}

############################################################################################

## Collect Initial Dataset ##
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Pretrain.csv')
raw.data.acquisition = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Acquisition.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PAL\\Weston Data PAL Main.csv')

## Run Namecheck on Data ##
name.data.pretrain = Naming.Pretrain.Function(raw.data.pretrain)
name.data.acquisition = Naming.MainAcq.Function(raw.data.acquisition)
name.data.main = Naming.MainAcq.Function(raw.data.main)

## Run QC Analysis
qc.data.pretrain = QC.Pretrain.Function(name.data.pretrain)

qc.idlist.pretrain = as.vector(unique(as.character(qc.data.pretrain$AnimalID)))
