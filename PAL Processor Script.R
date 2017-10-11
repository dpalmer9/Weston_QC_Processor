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
bad.list$bad.3x = c('3XTG', ' 3XTG', '3xTG','3TG') #5

bad.list$bad.wt = c(' w', 'w','W', ' W') #6
bad.list$bad.tg = c(' t', 't','T', ' T') #7

bad.list$bad.app.wt = c(' w', 'w','W') #8 
bad.list$bad.app.tg = c(' t', 't','T') #9
bad.list$bad.5x.wt = c(' w', 'w','W') #10
bad.list$bad.5x.tg = c(' t', 't','T') #11
bad.list$bad.3x.wt = c(' w', 'w','W') #12
bad.list$bad.3x.tg = c(' t', 't','T') #13

bad.list$bad.female = c('F') #14
bad.list$bad.male = c('M') #15

bad.list$bad.4m = c(4,'4') #16
bad.list$bad.7m = c(7,'7') #17
bad.list$bad.10m = c(10,'10') #18

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
good.list[32] = 'sPAL1'

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
QC.Pretrain.Function = function(dataset, good.names=good.list){
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
  good.10m = as.vector(unlist(good.names[c(26,27,28,30,31,32)]))
  new.data = new.data[which(new.data[ ,6] == "4"), ]
  new.data = new.data[!new.data[ ,7] %in% good.10m, ]
  return(new.data)
}

# QC Function PAL Acq #
QC.Acq.Function = function(dataset, good.names=good.list){
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
  good.4m = as.vector(unlist(good.names[c(21,22,23,24,25,29)]))
  good.10m = as.vector(unlist(good.names[c(26,27,28,30,31,32)]))
  new.data.4m = new.data[which(new.data[ ,7] == "4"), ]
  new.data.10m = new.data[which(new.data[ ,7] == "10"), ]
  new.data.4m = new.data.4m[!new.data.4m[ ,8] %in% good.10m, ]
  new.data.10m = new.data.10m[!new.data.10m[ ,8] %in% good.4m, ]
  new.data = rbind(new.data.4m, new.data.10m)
  new.data = new.data[which((((new.data[ ,15] == 36 ) | (new.data[ ,15] == 30 )) & (new.data[ ,14] < 3600)) | ((new.data[ ,15] < 36 ) & (new.data[ ,14] > 3599))), ]
  return(new.data)
}
# QC Function PAL Main #

QC.Main.Function = function(dataset, good.names=good.list){
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
  good.4m = as.vector(unlist(good.names[c(21,22,23,24,25,29)]))
  good.10m = as.vector(unlist(good.names[c(26,27,28,30,31,32)]))
  new.data.4m = new.data[which(new.data[ ,7] == "4"), ]
  new.data.10m = new.data[which(new.data[ ,7] == "10"), ]
  new.data.4m = new.data.4m[!new.data.4m[ ,8] %in% good.10m, ]
  new.data.10m = new.data.10m[!new.data.10m[ ,8] %in% good.4m, ]
  new.data = rbind(new.data.4m, new.data.10m)
  new.data = new.data[which((((new.data[ ,15] == 36 ) | (new.data[ ,15] == 30 )) & (new.data[ ,14] < 3600)) | ((new.data[ ,15] < 36 ) & (new.data[ ,14] > 3599))), ]
  return(new.data)
}

# Fix Naming Rules - Main & Acq #
Naming.MainAcq.Function = function(dataset, good.names=good.list, bad.names=bad.list){
  colnames(dataset)[1:11] = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day','Week')
  fixed.colnames = colnames(dataset)[13:317]
  #block.start = list(22:24,27:29,32:34,37:39,42:44,47:49)
  block.start = c(22,27,32,37,42,47)
  for(a in 1:length(block.start)){
    block.start[a] = block.start[a] - 12
  }
  #lat.start = list(52:87)
  lat.start = c(52,90,128,166,204,242,280)
  for(a in 1:length(lat.start)){
    lat.start[a] = lat.start[a] - 12
  }
  block.sub = FALSE
  lat.sub = FALSE
  block.num = 0
  lat.num = 0
  for(a in 1:length(fixed.colnames)){
    fixed.colnames[a] = gsub("End.Summary...", "", fixed.colnames[a], ignore.case=FALSE)
    fixed.colnames[a] = gsub("X12", "Twelve", fixed.colnames[a],ignore.case=FALSE)
    fixed.colnames[a] = gsub("..s.", ".", fixed.colnames[a],ignore.case=FALSE)
    if(is.element(a,block.start) == TRUE){
      block.sub = TRUE
      block.num = 1
    }
    if(is.element(a,lat.start) == TRUE){
      lat.sub = TRUE
      lat.num = 1
    }
    if(block.sub == TRUE){
      if((block.num <= 3) & (block.num > 1)){
        fixed.colnames[a] = substr(fixed.colnames[a],1,(nchar(fixed.colnames[a]) - 2))
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(block.num), sep = "")
        block.num = block.num + 1
      }else if(block.num == 1){
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(block.num), sep = ".")
        block.num = block.num + 1
      }else{
        block.sub = FALSE
        block.num = 0
      }
    }
    if(lat.sub == TRUE){
      if((lat.num <= 36) & (lat.num > 1)){
        fixed.colnames[a] = substr(fixed.colnames[a],1,(nchar(fixed.colnames[a]) - 2))
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(lat.num), sep = "")
        lat.num = lat.num + 1
      }else if(lat.num == 1){
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(lat.num), sep = ".")
        lat.num = lat.num + 1
      }else{
        lat.sub = FALSE
        lat.num = 0
      }
    }
  }
  colnames(dataset)[13:317] = fixed.colnames
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
        if(isTRUE(dataset[b,3] == "APP-PS1")){
          dataset[b,col.num] = good.list[8]
        }else if(isTRUE(dataset[b,3] == "5xFAD")){
          dataset[b,col.num] = good.list[10]
        }else if(isTRUE(dataset[b,3] == "3xTG-AD")){
          dataset[b,col.num] = good.list[12]
        }
      }
      if(isTRUE(dataset[b,col.num] == 'T')){
        if(isTRUE(dataset[b,3] == "APP-PS1")){
          dataset[b,col.num] = good.list[9]
        }else if(isTRUE(dataset[b,3] == "5xFAD")){
          dataset[b,col.num] = good.list[11]
        }else if(isTRUE(dataset[b,3] == "3xTG-AD")){
          dataset[b,col.num] = good.list[13]
        }
      }
    }
  }
  return(dataset)
}

# Fix Date #
Datefix.Function <- function(data, colnum){
  a<-list()
  formats = c('%Y%m%d','%m/%d/%Y')
  for(i in 1:length(formats)){
    a[[i]]<- as.Date(as.character(data[ ,colnum]),format=formats[i])
    a[[1]][!is.na(a[[i]])]<-a[[i]][!is.na(a[[i]])]
  }
  data[ ,colnum] = a
  return(data)
}

# Fix Latency Calc #
LatFix.MainAcq.Function = function(dataset,IQD.num){
  mean.colnums = c(126,164,202,240,278,316)
  for(a in 1:length(mean.colnums)){
    lat.raw.cols = c((mean.colnums[a] - 36):(mean.colnums[a] - 1))
    mean.col = mean.colnums[a]
    std.col = mean.colnums[a] + 1
    for(b in 1:nrow(dataset)){
      lat.data = as.vector(as.numeric(dataset[b,lat.raw.cols]))
      lat.iqr = IQR(lat.data,na.rm=TRUE)
      lat.iqr.quartile = quantile(lat.data,na.rm=TRUE)
      lat.iqd.upper = lat.iqr.quartile[4] + (IQD.num * lat.iqr)
      lat.iqd.lower = lat.iqr.quartile[2] - (IQD.num * lat.iqr)
      lat.data = lat.data[lat.data > lat.iqd.lower]
      lat.data = lat.data[lat.data < lat.iqd.upper]
      lat.newmean = mean(lat.data, na.rm=TRUE)
      lat.newsd = sd(lat.data, na.rm=TRUE)
      dataset[b,mean.col] = lat.newmean
      dataset[b,std.col] = lat.newsd
    }
    #mean.avg.data = as.vector(as.numeric(dataset[ ,mean.col]))
    #mean.iqr = IQR(mean.avg.data, na.rm=TRUE)
    #mean.iqd.multiplied = mean.iqr * IQD.num
    #mean.iqr.quartile = quartile(mean.avg.data,na.rm=TRUE)
    #mean.iqr.quartile25 = as.numeric(mean.iqr.quartile[2])
    #mean.iqr.quartile75 = as.numeric(mean.iqr.quartile[4])
    #mean.avg.data[mean.avg.data < (mean.iqr.quartile25 * mean.iqd.multiplied)] = NA
    #mean.avg.data[mean.avg.data > (mean.iqr.quartile75 * mean.iqd.multiplied)] = NA
    #dataset[ ,mean.col] = mean.avg.data
  }
  return(dataset)
}

# Aggregate Function #
Aggregate.Function = function(dataset){
  agg.list = list(dataset$AnimalID, dataset$TestSite, dataset$Mouse.Strain, dataset$Genotype, dataset$Sex, dataset$Age.Months, dataset$Task, dataset$Week)
  agg.data = aggregate(dataset, by= agg.list, FUN=mean, na.rm = TRUE)
  agg.data[ ,10:16] = agg.data[ ,1:7]
  agg.data[ ,19] = agg.data[ ,8]
  agg.data[ ,c(1:9,18)] = NULL
  return(agg.data)
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

## Fix Dates for Files ##
date.data.pretrain = Datefix.Function(name.data.pretrain,8)
date.data.acquisition = Datefix.Function(name.data.acquisition,9)
date.data.main = Datefix.Function(name.data.main,9)

## Run QC Analysis
qc.data.pretrain = QC.Pretrain.Function(date.data.pretrain)

qc.idlist.pretrain = as.vector(unique(as.character(qc.data.pretrain$AnimalID)))

qc.data.acquisition = QC.Acq.Function(date.data.acquisition)

qc.data.main = QC.Main.Function(date.data.main)

## Lat Fix ##
qc.data.acquisition.lat = LatFix.MainAcq.Function(qc.data.acquisition,3)
qc.data.main.lat = LatFix.MainAcq.Function(qc.data.main,3)

## Remove Non-Completers ##

qc.idlist.final = as.vector(unique(as.character(qc.data.main.lat$AnimalID)))
qc.data.acquisition.final = qc.data.acquisition.lat[qc.data.acquisition.lat$AnimalID %in% qc.idlist.final, ]
qc.data.pretrain.final = qc.data.pretrain[qc.data.pretrain$AnimalID %in% qc.idlist.final, ]

## Session Only Acq Data #

qc.data.acq.sessions = qc.data.acquisition.final[ ,1:10]
qc.data.acq.sessions[ ,10] = 1
qc.data.acq.sessions[ ,9] = NULL
qc.data.acq.agg = aggregate(Day ~ Database + AnimalID + TestSite + Mouse.Strain + Genotype + Sex + Age.Months + Task, FUN=sum, na.rm=TRUE, data=qc.data.acq.sessions)


## Aggregate Main File ##

qc.data.main.agg = Aggregate.Function(qc.data.main.lat)


## Save Raw Data Files ##
write.csv(qc.data.pretrain.final, "Weston PAL Pretraining Oct 11 2017 Updated.csv")
write.csv(qc.data.acquisition.final, "Weston PAL Acquisition Oct 11 2017 Updated.csv")
write.csv(qc.data.acq.agg, "Weston PAL Acquisition Aggregated Oct 11 2017 Updated.csv")
write.csv(qc.data.main.lat, "Weston PAL Main Task Oct 11 2017 Updated.csv")
write.csv(qc.data.main.agg, "Weston PAL Main Task Aggregated Oct 11 2017 Updated.csv")
