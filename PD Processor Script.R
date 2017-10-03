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

bad.list$bad.wt = c(' w', 'w') #6
bad.list$bad.tg = c(' t', 't') #7

bad.list$bad.app.wt = c(' w', 'w', 'C57Bl/6j') #8 
bad.list$bad.app.tg = c(' t', 't') #9
bad.list$bad.5x.wt = c(' w', 'w') #10
bad.list$bad.5x.tg = c(' t', 't') #11
bad.list$bad.3x.wt = c(' w', 'w') #12
bad.list$bad.3x.tg = c(' t', 't') #13

bad.list$bad.female = c('F') #14
bad.list$bad.male = c('M') #15

bad.list$bad.4m = c(4,'4') #16
bad.list$bad.7m = c(7,'7') #17
bad.list$bad.10m = c(10,'10') #18

bad.list$bad.hab1 = c('Hab1') #19
bad.list$bad.hab2 = c('Hab2') #20
bad.list$bad.IT = c('PD_Initial_Touch_Training_v3') #21
bad.list$bad.MT = c('PD_Must_Touch_Training_v3') #22
bad.list$bad.MI = c('PD_Must_Initiate_Training_v3') #23
bad.list$bad.punish = c('PD_Punish_Incorrect_Training_v3') #24

bad.list$bad.acq1 = c('PD_Acquisition_1_v3') #25
bad.list$bad.acq2 = c('PD_Acquisition_2_v3') #26
bad.list$bad.acq3 = c('PD_Acquisition_3_v3') #27

bad.list$bad.base1 = c('PD_Baseline_1_v3') #28
bad.list$bad.base2 = c('PD_Baseline_2_v3') #29
bad.list$bad.base3 = c('PD_Baseline_3_v3') #30

bad.list$bad.rev1 = c('PD_Reversal_1_v3') #31
bad.list$bad.rev2 = c('PD_Reversal_2_v3') #32
bad.list$bad.rev3 = c('PD_Reversal_3_v3') #33

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
good.list[25] = 'PD Acquisition 1'
good.list[26] = 'PD Acquisition 2'
good.list[27] = 'PD Acquisition 3'
good.list[28] = 'PD Baseline 1'
good.list[29] = 'PD Baseline 2'
good.list[30] = 'PD Baseline 3'
good.list[31] = 'PD Reversal 1'
good.list[32] = 'PD Reversal 2'
good.list[33] = 'PD Reversal 3'

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
  good.10m = as.vector(unlist(good.names[c(25:33)]))
  new.data = new.data[which(new.data[ ,7] == "4"), ]
  new.data = new.data[!new.data[ ,8] %in% good.10m, ]
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
  good.4m = as.vector(unlist(good.names[c(25)]))
  good.7m = as.vector(unlist(good.names[c(26)]))
  good.10m = as.vector(unlist(good.names[c(27)]))
  new.data.4m = new.data[which(new.data[ ,7] == "4"), ]
  new.data.7m = new.data[which(new.data[ ,7] == "7"), ]
  new.data.10m = new.data[which(new.data[ ,7] == "10"), ]
  new.data.4m = new.data.4m[!new.data.4m[ ,8] %in% c(good.7m,good.10m), ]
  new.data.7m = new.data.7m[!new.data.7m[ ,8] %in% c(good.4m,good.10m), ]
  new.data.10m = new.data.10m[!new.data.10m[ ,8] %in% c(good.4m,good.7m), ]
  new.data = rbind(new.data.4m, new.data.7m)
  new.data = rbind(new.data,new.data.10m)
  new.data = new.data[which(((new.data[ ,15] == 30 ) & (new.data[ ,14] < 3600)) | ((new.data[ ,15] < 30 ) & (new.data[ ,14] > 3599))), ]
  final.data = as.data.frame(matrix(nrow=0,ncol=ncol(new.data)))
  for(a in 1:length(initial.ID.list)){
    temp.id = initial.ID.list[a]
    temp.data = new.data[which(new.data[ ,2] == temp.id), ]
    temp.4m = temp.data[which(temp.data[ ,7] == "4"), ]
    temp.7m = temp.data[which(temp.data[ ,7] == "7"), ]
    temp.10m = temp.data[which(temp.data[ ,7] == "10"), ]
    criteria.check.func = function(dataset){
      date.org = dataset[order(dataset[ ,9]), ]
      if(nrow(dataset) > 0){
        checked.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
        criteria.count = 0
        for(b in 1:nrow(date.org)){
          if(criteria.count == 2){
            break
          }
          current.trials = date.org[b,15]
          current.acc = date.org[b,17]
          if(isTRUE((current.trials == 30) & (current.acc >= 80))){
            criteria.count = criteria.count + 1
          }else{
            criteria.count = 0
          }
          checked.data = rbind(checked.data,date.org[b, ])
        }
        return(checked.data)
      }else{
        return(date.org)
      }
    }
    temp.4m = criteria.check.func(temp.4m)
    temp.7m = criteria.check.func(temp.7m)
    temp.10m = criteria.check.func(temp.10m)
    temp.data = rbind(temp.4m,temp.7m)
    temp.data = rbind(temp.data,temp.10m)
    final.data = rbind(final.data,temp.data)
  }
  
  return(final.data)
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
  good.4m = as.vector(unlist(good.names[c(28,31)]))
  good.7m = as.vector(unlist(good.names[c(29,32)]))
  good.10m = as.vector(unlist(good.names[c(30,33)]))
  new.data.4m = new.data[which(new.data[ ,7] == "4"), ]
  new.data.7m = new.data[which(new.data[ ,7] == "7"), ]
  new.data.10m = new.data[which(new.data[ ,7] == "10"), ]
  new.data.4m = new.data.4m[!new.data.4m[ ,8] %in% c(good.7m,good.10m), ]
  new.data.7m = new.data.7m[!new.data.7m[ ,8] %in% c(good.4m,good.10m), ]
  new.data.10m = new.data.10m[!new.data.10m[ ,8] %in% c(good.4m,good.7m), ]
  new.data = rbind(new.data.4m, new.data.7m)
  new.data = rbind(new.data,new.data.10m)
  new.data[ ,14] = as.character(new.data[ ,14])
  new.data[ ,14] = gsub(",", "", new.data[ ,14])
  new.data[ ,14] = as.numeric(new.data[ ,14])
  new.data2 = new.data[which(((new.data[ ,15] == 30 ) & (new.data[ ,14] < 3600)) | ((new.data[ ,15] < 30 ) & (new.data[ ,14] > 3599))), ]
  return(new.data)
}

# Fix Naming Rules - Main & Acq #
Naming.MainAcq.Function = function(dataset, good.names=good.list, bad.names=bad.list){
  colnames(dataset)[1:10] = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day')
  fixed.colnames = colnames(dataset)[12:81]
  lat.start = c(18,50)
  for(a in 1:length(lat.start)){
    lat.start[a] = lat.start[a] - 11
  }
  lat.sub = FALSE
  lat.num = 0
  for(a in 1:length(fixed.colnames)){
    fixed.colnames[a] = gsub("End.Summary...", "", fixed.colnames[a], ignore.case=FALSE)
    fixed.colnames[a] = gsub("..s.", ".", fixed.colnames[a],ignore.case=FALSE)
    if(is.element(a,lat.start) == TRUE){
      lat.sub = TRUE
      lat.num = 1
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
  colnames(dataset)[12:81] = fixed.colnames
  col.list.spacefix = c(1,2,3,4,5,6,7,8,10,11)
  for(a in col.list.spacefix){
    dataset[ ,a] = gsub(" ", "", dataset[ ,a])
  }
  if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) > 2)){
    geno.long = TRUE
    position.vec = c(1:5,8:13,14:33)
  }else if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) == 2)){
    geno.long = FALSE
    position.vec = c(1:5,6:7,14:33)
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
  colnames(dataset) = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day')
  col.list.spacefix = c(1,2,3,4,5,6,7,8,9,10)
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

# Fix Date #
Datefix.Function <- function(data, colnum){
  a<-list()
  formats = c('%Y%m%d','%m/%d/%Y','%Y-%m-%d %H:%M')
  for(i in 1:length(formats)){
    a[[i]]<- as.Date(as.character(data[ ,colnum]),format=formats[i])
    a[[1]][!is.na(a[[i]])]<-a[[i]][!is.na(a[[i]])]
  }
  data[ ,colnum] = a
  return(data)
}

# Fix Latency Calc #
LatFix.MainAcq.Function = function(dataset,IQD.num){
  mean.colnums = c(48,80)
  for(a in 1:length(mean.colnums)){
    lat.raw.cols = c((mean.colnums[a] - 30):(mean.colnums[a] - 1))
    mean.col = mean.colnums[a]
    std.col = mean.colnums[a] + 1
    for(b in 1:nrow(dataset)){
      lat.data = as.vector(as.numeric(dataset[b,lat.raw.cols]))
      lat.iqr = IQR(lat.data,na.rm=TRUE)
      lat.mean = mean(lat.data, na.rm=TRUE)
      lat.iqd.upper = lat.mean + (IQD.num * lat.iqr)
      lat.iqd.lower = lat.mean - (IQD.num * lat.iqr)
      lat.data = lat.data[lat.data > lat.iqd.lower]
      lat.data = lat.data[lat.data < lat.iqd.upper]
      lat.newmean = mean(lat.data, na.rm=TRUE)
      lat.newsd = sd(lat.data, na.rm=TRUE)
      dataset[b,mean.col] = lat.newmean
      dataset[b,std.col] = lat.newsd
    }
    mean.avg.data = as.vector(as.numeric(dataset[ ,mean.col]))
    mean.iqr = IQR(mean.avg.data, na.rm=TRUE)
    mean.mean = mean(mean.avg.data, na.rm=TRUE)
    mean.iqd.upper = mean.mean + (IQD.num * mean.iqr)
    mean.iqd.lower = mean.mean - (IQD.num & mean.iqr)
    mean.avg.data[mean.avg.data < mean.iqd.lower] = NA
    mean.avg.data[mean.avg.data > mean.iqd.upper] = NA
    dataset[ ,mean.col] = mean.avg.data
  }
  return(dataset)
}


############################################################################################

## Collect Initial Dataset ##
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PD\\Weston Data PD Pretraining.csv')
raw.data.acquisition = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PD\\Weston Data PD Acquisition.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\PD\\Weston Data PD Baseline Reversal.csv')

## Run Namecheck on Data ##
name.data.pretrain = Naming.Pretrain.Function(raw.data.pretrain)
name.data.acquisition = Naming.MainAcq.Function(raw.data.acquisition)
name.data.main = Naming.MainAcq.Function(raw.data.main)

## Fix Dates for Files ##
date.data.pretrain = Datefix.Function(name.data.pretrain,9)
date.data.acquisition = Datefix.Function(name.data.acquisition,9)
date.data.main = Datefix.Function(name.data.main,9)

## Run QC Analysis
qc.data.pretrain = QC.Pretrain.Function(date.data.pretrain)

qc.idlist.pretrain = as.vector(unique(as.character(qc.data.pretrain$AnimalID)))

qc.data.acquisition = QC.Acq.Function(date.data.acquisition)

qc.data.main = QC.Main.Function(date.data.main)

## Run LatFix ##
qc.data.acquisition.lat = LatFix.MainAcq.Function(qc.data.acquisition,3)
qc.data.main.lat = LatFix.MainAcq.Function(qc.data.main,3)


## Save Raw Data Files ##
write.csv(qc.data.pretrain, "Weston PD Pretrain QC Oct 3 2017 NEW LATENCY.csv")
write.csv(qc.data.acquisition, "Weston PD Acquisition QC Oct 3 2017 NEW LATENCY.csv")
write.csv(qc.data.main, "Weston PD Basline Reversal QC Oct 3 2017 NEW LATENCY.csv")
