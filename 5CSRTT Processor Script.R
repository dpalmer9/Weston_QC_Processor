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

bad.list$bad.app = c('APPPS1', ' APPPS1', 'APP') #3
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
bad.list$bad.IT = c('5CSRTT_Initial_Touch_Training_v3') #21
bad.list$bad.MT = c('5CSRTT_Must_Touch_Training_v2') #22
bad.list$bad.MI = c('5CSRTT_Must_Initiate_Training_v1', ' 5CSRTT_Must_Initiate_training', '5CSRTT_Must_Initiate_training') #23
bad.list$bad.punish = c('5CSRTT_Punish_Incorrect_Training_v3') #24

bad.list$bad.t4 = c('4000ms','5CSRTT_4000ms_Var1') #25
bad.list$bad.t2 = c('2000ms', '5CSRTT_2000ms_Var1') #26

bad.list$bad.pr15 = c('1500ms') #27
bad.list$bad.pr10 = c('1000ms') #28
bad.list$bad.pr08 = c('0800ms') #29
bad.list$bad.pr06 = c('0600ms') #30



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
good.list[25] = '4'
good.list[26] = '2'
good.list[27] = '1.5'
good.list[28] = '1'
good.list[29] = '0.8'
good.list[30] = '0.6'


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
  good.10m = as.vector(unlist(good.names[c(25:30)]))
  new.data = new.data[which(new.data[ ,7] == "4"), ]
  new.data = new.data[!new.data[ ,8] %in% good.10m, ]
  return(new.data)
}

# QC Function PAL Acq #
QC.Acq.Function = function(dataset, good.names=good.list){
  dataset = dataset[which(dataset[ ,7] != '13_15'), ]
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
  new.data = new.data[which(((new.data[ ,14] == 50 ) & (new.data[ ,13] < 3600)) | ((new.data[ ,14] < 50 ) & (new.data[ ,13] > 3599))), ]
  final.data = as.data.frame(matrix(nrow=0,ncol=ncol(new.data)))
  for(a in 1:length(initial.ID.list)){
    temp.id = initial.ID.list[a]
    temp.data = new.data[which(new.data[ ,2] == temp.id), ]
    temp.4sec = temp.data[which(temp.data[ ,8] == "4"), ]
    temp.2sec = temp.data[which(temp.data[ ,8] == "2"), ]
    criteria.check.func = function(dataset,trainstim){
      date.org = dataset[order(dataset[ ,9]), ]
      if(nrow(dataset) > 0){
        checked.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
        criteria.count = 0
        for(b in 1:nrow(date.org)){
          if(criteria.count == 3){
            break
          }
          current.trials = date.org[b,14]
          current.acc = date.org[b,15]
          current.omit = date.org[b,16]
          if(trainstim == 4){
            if(isTRUE((current.trials >= 30) & (current.acc >= 80) & (current.omit <= 20))){
              criteria.count = criteria.count + 1
            }else{
              criteria.count = 0
            }
          }else if(trainstim == 2){
            if(isTRUE((current.trials == 50) & (current.acc >= 80) & (current.omit <= 20))){
              criteria.count = criteria.count + 1
            }else{
              criteria.count = 0
            }
          }
          checked.data = rbind(checked.data,date.org[b, ])
        }
        return(checked.data)
      }else{
        return(date.org)
      }
    }
    temp.4sec = criteria.check.func(temp.4sec,4)
    temp.2sec = criteria.check.func(temp.2sec,2)
    temp.data = rbind(temp.4sec,temp.2sec)
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
  new.data2 = new.data[which(((new.data[ ,14] == 50 ) & (new.data[ ,13] < 3600)) | ((new.data[ ,14] < 50 ) & (new.data[ ,13] > 3599))), ]
  final.data = as.data.frame(matrix(ncol=ncol(new.data2),nrow=0))
  colnames(final.data) = colnames(new.data2)
  for(a in 1:length(initial.ID.list)){
    temp.data = new.data2[which(new.data2[ ,2] == initial.ID.list[a]), ]
    temp.4m = temp.data[which(temp.data[ ,7] == "4"), ]
    temp.7m = temp.data[which(temp.data[ ,7] == "7"), ]
    temp.10m = temp.data[which(temp.data[ ,7] == "10"), ]
    probe.check.func = function(dataset){
      probe.vec = c('1.5','1','0.8','0.6')
      count.vec = c()
      count.vec[1] = nrow(dataset[which(dataset[ ,8] == probe.vec[1]), ])
      count.vec[2] = nrow(dataset[which(dataset[ ,8] == probe.vec[2]), ])
      count.vec[3] = nrow(dataset[which(dataset[ ,8] == probe.vec[3]), ])
      count.vec[4] = nrow(dataset[which(dataset[ ,8] == probe.vec[4]), ])
      return.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
      colnames(return.data) = colnames(dataset)
      for(b in 1:4){
        if(count.vec[b] > 2){
          probe.data = dataset[which(dataset[ ,8] == probe.vec[b]), ]
          probe.data = probe.data[order(probe.data$Date), ]
          probe.data = probe.data[c(1:2), ]
        }else{
          probe.data = dataset[which(dataset[ ,8] == probe.vec[b]), ]
        }
        return.data = rbind(return.data,probe.data)
      }
      return(return.data)
    }
    check.4m = probe.check.func(temp.4m)
    check.7m = probe.check.func(temp.7m)
    check.10m = probe.check.func(temp.10m)
    final.id.data = rbind(check.4m,check.7m)
    final.id.data = rbind(final.id.data,check.10m)
    final.data = rbind(final.data,final.id.data)
  }
  return(final.data)
}

# Fix Naming Rules - Main & Acq #
Naming.MainAcq.Function = function(dataset, good.names=good.list, bad.names=bad.list,type){
  probelist = c('Probe','probe','p',2,'P')
  dataset[ ,1] = NULL
  if(is.element(type,probelist) == TRUE){
    dataset[ ,c(175,227)] = NULL
  }
  colnames(dataset)[1:10] = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Stimulus.Length','Date','Day')
  fixed.colnames = colnames(dataset)[12:226]
  lat.start = c(19,71,123,175)
  for(a in 1:length(lat.start)){
    lat.start[a] = lat.start[a] - 11
  }
  lat.sub = FALSE
  lat.num = 0
  for(a in 1:length(fixed.colnames)){
    fixed.colnames[a] = gsub("Threshold...", "", fixed.colnames[a], ignore.case=FALSE)
    fixed.colnames[a] = gsub("Trial.Analysis...", "", fixed.colnames[a], ignore.case=FALSE)
    if(is.element(a,lat.start) == TRUE){
      lat.sub = TRUE
      lat.num = 1
    }
    if(lat.sub == TRUE){
      if((lat.num <= 50) & (lat.num > 9)){
        fixed.colnames[a] = substr(fixed.colnames[a],1,(nchar(fixed.colnames[a]) - 3))
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(lat.num), sep = ".")
        lat.num = lat.num + 1
      }else if((lat.num <= 9) & (lat.num > 1)){
        fixed.colnames[a] = substr(fixed.colnames[a],1,(nchar(fixed.colnames[a]) - 2))
        fixed.colnames[a] = paste(fixed.colnames[a], as.character(lat.num), sep = ".")
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
  colnames(dataset)[12:226] = fixed.colnames
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
  dataset[ ,c(1,12,13)] = NULL
  colnames(dataset) = c('Database','AnimalID','TestSite','Mouse.Strain','Genotype','Sex','Age.Months','Task','Date','Day')
  dataset = dataset[which(dataset[ ,8] != '5CSRTT_Habituation_1_'), ]
  dataset = dataset[which(dataset[ ,8] != '5CSRTT_Habituation_2_'), ]
  col.list.spacefix = c(1,2,3,4,5,6,7,8,10)
  for(a in col.list.spacefix){
    dataset[ ,a] = gsub(" ", "", dataset[ ,a])
  }
  if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) > 2)){
    geno.long = TRUE
    position.vec = c(1:5,8:13,14:30)
  }else if(isTRUE(length(as.vector(unique(as.character(dataset[ ,5])))) == 2)){
    geno.long = FALSE
    position.vec = c(1:5,6:7,14:30)
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
  formats = c('%Y%m%d','%m/%d/%Y','%Y-%m-%d %H:%M','%Y-%m-%d')
  for(i in 1:length(formats)){
    a[[i]]<- as.Date(as.character(data[ ,colnum]),format=formats[i])
    a[[1]][!is.na(a[[i]])]<-a[[i]][!is.na(a[[i]])]
  }
  data[ ,colnum] = a
  return(data)
}

Datefix.Conversion.Function = function(dataset,colnum){
  bad.start.date = as.Date('1899-12-30',format = "%Y-%m-%d")
  dataset[ ,colnum] = as.character(dataset[ ,colnum])
  new.list = c()
  for(a in 1:nrow(dataset)){
    curr.value = dataset[a,colnum]
    if(isTRUE(grepl("/",curr.value))){
      dataset[a,colnum] = as.character(curr.value)
    }else{
      #fixed.value = bad.start.date + as.numeric(curr.value)
      dataset[a,colnum] = as.character(as.Date.numeric(as.numeric(curr.value),origin=bad.start.date))
    }
  }
  return(dataset)
}
# Fix Latency Decimal #
LatFix.Decimal.MainAcq.Function = function(dataset,fix.cols){
  dataset[ ,fix.cols] = apply(dataset[ ,fix.cols],c(1,2), function(x) as.numeric(x)/100000)
  return(dataset)
}

# Fix Latency Calc #
LatFix.Acq.Function = function(dataset,IQD.num, good.names=good.list){
  mean.colnums = c(173,225)
  for(a in 1:length(mean.colnums)){
    lat.raw.cols = c((mean.colnums[a] - 50):(mean.colnums[a] - 1))
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
    new.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
    colnames(new.data) = colnames(dataset)
    for(b in 1:2){
      curr.site = good.names[(b)]
      for(c in 1:6){
        curr.geno = good.names[(c + 7)]
        for(d in 1:2){
          curr.sex = good.names[(d + 13)]
          for(e in 1:3){
            curr.age = good.names[(e + 15)]
            for(f in 1:4){
              curr.probe = good.names[(f + 26)]
              section.data = dataset[which((dataset[ ,3] == curr.site) & (dataset[ ,5] == curr.geno) & (dataset[ ,6] == curr.sex) & (dataset[ ,7] == curr.age) & (dataset[ ,8] == curr.probe)), ]
              for(a in 1:length(mean.colnums)){
                mean.col = mean.colnums[a]
                mean.avg.data = as.vector(as.numeric(section.data[ ,mean.col]))
                mean.iqr = IQR(mean.avg.data, na.rm=TRUE)
                mean.iqd.multiplied = mean.iqr * IQD.num
                mean.iqr.quartile = quantile(mean.avg.data,na.rm=TRUE)
                mean.avg.data[mean.avg.data < (mean.iqr.quartile[2] - mean.iqd.multiplied)] = NA
                mean.avg.data[mean.avg.data > (mean.iqr.quartile[4] + mean.iqd.multiplied)] = NA
                section.data[ ,mean.col] = mean.avg.data
              }
              new.data = rbind(new.data,section.data)
            }
          }
        }
      }
    }
  }
  return(dataset)
}

LatFix.Probe.Function = function(dataset,IQD.num,good.names=good.list){
  mean.colnums = c(69,121)
  for(a in 1:length(mean.colnums)){
    lat.raw.cols = c((mean.colnums[a] - 50):(mean.colnums[a] - 1))
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
    new.data = as.data.frame(matrix(nrow=0,ncol=ncol(dataset)))
    colnames(new.data) = colnames(dataset)
    for(b in 1:2){
      curr.site = good.names[(b)]
      for(c in 1:6){
        curr.geno = good.names[(c + 7)]
        for(d in 1:2){
          curr.sex = good.names[(d + 13)]
          for(e in 1:3){
            curr.age = good.names[(e + 15)]
            for(f in 1:4){
              curr.probe = good.names[(f + 26)]
              section.data = dataset[which((dataset[ ,3] == curr.site) & (dataset[ ,5] == curr.geno) & (dataset[ ,6] == curr.sex) & (dataset[ ,7] == curr.age) & (dataset[ ,8] == curr.probe)), ]
              for(a in 1:length(mean.colnums)){
                mean.col = mean.colnums[a]
                mean.avg.data = as.vector(as.numeric(section.data[ ,mean.col]))
                mean.iqr = IQR(mean.avg.data, na.rm=TRUE)
                mean.iqd.multiplied = mean.iqr * IQD.num
                mean.iqr.quartile = quantile(mean.avg.data,na.rm=TRUE)
                mean.avg.data[mean.avg.data < (mean.iqr.quartile[2] - mean.iqd.multiplied)] = NA
                mean.avg.data[mean.avg.data > (mean.iqr.quartile[4] + mean.iqd.multiplied)] = NA
                section.data[ ,mean.col] = mean.avg.data
              }
              new.data = rbind(new.data,section.data)
            }
          }
        }
      }
    }
  }
  return(new.data)
}

Vigilance.FixRaw.Function = function(dataset){
  total.tri.col = 14
  acc.range = 15:64
  omission.range = 67:116
  for(a in 1:nrow(dataset)){
    acc.vec = as.vector(as.numeric(dataset[a,acc.range]))
    acc.vec = acc.vec[!is.na(acc.vec)]
    omission.vec = as.vector(as.numeric(dataset[a,omission.range]))
    omission.vec = omission.vec[!is.na(omission.vec)]
    total.trials = as.numeric(dataset[a,total.tri.col])
    acc.len = length(acc.vec)
    omission.len = length(omission.vec)
    if(isTRUE(acc.len < total.trials)){
      acc.discrepancy = total.trials - acc.len
      acc.zerofill = c()
      for(b in 1:acc.discrepancy){
        acc.zerofill[b] = 0
      }
      acc.vec = c(acc.zerofill,acc.vec)
      if(length(acc.vec) < 50){
        na.discrepancy = 50 - length(acc.vec)
        na.fill = c()
        for(b in 1:na.discrepancy){
          na.fill[b] = NA
        }
        acc.vec = c(acc.vec,na.fill)
      }
      dataset[a,acc.range] = acc.vec
    }
    if(isTRUE(omission.len < total.trials)){
      omission.discrepancy = total.trials - omission.len
      omission.zerofill = c()
      for(b in 1:omission.discrepancy){
        omission.zerofill[b] = 0
      }
      omission.vec = c(omission.zerofill,omission.vec)
      if(length(omission.vec) < 50){
        na.discrepancy = 50 - length(omission.vec)
        na.fill = c()
        for(b in 1:na.discrepancy){
          na.fill[b] = NA
        }
        omission.vec = c(omission.vec,na.fill)
      }
      dataset[a,omission.range] = omission.vec
    }
  }
  return(dataset)
}

Vigilance.Calc.Block.Function = function(dataset,binsize){
  acc.start = 15
  omission.start = 67
  total.bins = 50 / binsize
  new.data = as.data.frame(matrix(nrow=nrow(dataset),ncol=(total.bins * 2)))
  bin.num = 1
  for(a in 1:total.bins){
    acc.col = a
    colnames(new.data)[acc.col] = paste("Accuracy.Block",a,sep=".")
    om.col = a + total.bins
    colnames(new.data)[om.col] = paste("Omission.Block",a,sep=".")
    for(b in 1:nrow(dataset)){
      acc.list = as.vector(dataset[b,c(acc.start:(acc.start + binsize - 1))])
      om.list = as.vector(dataset[b,c(omission.start:(omission.start + binsize - 1))])
      acc.bin = 0
      om.bin = 0
      na.count = 0
      if(a == 1){
        start.acc = 0
        start.om = 0
      }else{
        start.acc = as.numeric(dataset[b,(acc.start - 1)])
        start.om = as.numeric(dataset[b,(omission.start - 1)])
      }
      acc.track = c()
      om.track = c()
      for(c in 1:(binsize)){
        curr.acc = as.numeric(acc.list[c])
        curr.om = as.numeric(om.list[c])
        if(isTRUE(is.na(curr.acc) & is.na(curr.om))){
          acc.bin = acc.bin
          om.bin = om.bin
          na.count = na.count + 1
        }
        if(isTRUE(curr.acc > start.acc)){
          acc.bin = acc.bin + 1
        }else if(isTRUE(curr.om > start.om)){
          om.bin = om.bin + 1
        }
        start.acc = curr.acc
        start.om = curr.om
      }
      acc.bin = (acc.bin / (binsize - om.bin)) * 100
      om.bin = (om.bin / binsize) * 100
      if(na.count == binsize){
        acc.bin = NA
        om.bin = NA
      }
      new.data[b,acc.col] = acc.bin
      new.data[b,om.col] = om.bin
    }
    acc.start = (acc.start + binsize)
    omission.start = omission.start + binsize
  }
  final.data = cbind(dataset[ ,c(1:10,13:14)],new.data)
  return(final.data)
}

Vigilance.Calc.Average.Function = function(dataset,binsize){
  acc.start = 15
  omission.start = 67
  total.bins = 50 / binsize
  new.data = as.data.frame(matrix(nrow=nrow(dataset),ncol=(total.bins * 2)))
  bin.num = 1
  for(a in 1:total.bins){
    acc.col = a
    colnames(new.data)[acc.col] = paste("Accuracy.Block",a,sep=".")
    om.col = a + total.bins
    colnames(new.data)[om.col] = paste("Omission.Block",a,sep=".")
    for(b in 1:nrow(dataset)){
      acc.list = as.vector(as.numeric(dataset[b,c(acc.start:(acc.start + binsize - 1))]))
      om.list = as.vector(as.numeric(dataset[b,c(omission.start:(omission.start + binsize - 1))]))
      new.data[b,acc.col] = mean(acc.list, na.rm=TRUE)
      new.data[b,om.col] = mean(om.list, na.rm=TRUE)
    }
    acc.start = acc.start + binsize
    omission.start = omission.start + binsize
  }
  final.data = cbind(dataset[ ,c(1:10,13:14)],new.data)
  return(final.data)
}

Vigilance.Calc.Last.Function = function(dataset,binsize){
  acc.start = 15
  omission.start = 67
  total.bins = 50 / binsize
  new.data = as.data.frame(matrix(nrow=nrow(dataset),ncol=(total.bins * 2)))
  bin.num = 1
  for(a in 1:total.bins){
    acc.col = a
    colnames(new.data)[acc.col] = paste("Accuracy.Block",a,sep=".")
    om.col = a + total.bins
    colnames(new.data)[om.col] = paste("Omission.Block",a,sep=".")
    for(b in 1:nrow(dataset)){
      last.acc = as.numeric(dataset[b,(acc.start + binsize - 1)])
      last.om = as.numeric(dataset[b,(omission.start + binsize - 1)])
      if(isTRUE(is.na(last.acc))){
        for(c in (binsize - 1):1){
          last.acc = as.numeric(dataset[b,(acc.start + c - 1)])
          if(is.na(last.acc) == FALSE){
            break
          }
        }
      }
      if(isTRUE(is.na(last.om))){
        for(c in (binsize - 1):1){
          last.om = as.numeric(dataset[b,(omission.start + c - 1)])
          if(is.na(last.om) == FALSE){
            break
          }
        }
      }
      new.data[b,acc.col] = last.acc
      new.data[b,om.col] = last.om
    }
    acc.start = acc.start + binsize
    omission.start = omission.start + binsize
  }
  final.data = cbind(dataset[ ,c(1:10,13:14)],new.data)
  return(final.data)
}
############################################################################################

## Collect Initial Dataset ##
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\5CSRTT\\Weston Data 5CSRTT Pretraining.csv')
raw.data.acquisition = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\5CSRTT\\Weston Data 5CSRTT Acquisition.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_QC_Processor\\Data\\5CSRTT\\Weston Data 5CSRTT Probe.csv')

## Run Namecheck on Data ##
name.data.pretrain = Naming.Pretrain.Function(raw.data.pretrain)
name.data.acquisition = Naming.MainAcq.Function(raw.data.acquisition,type=1)
name.data.main = Naming.MainAcq.Function(raw.data.main,type=2)

## Fix Dates for Files ##
date.data.pretrain = Datefix.Function(name.data.pretrain,9)
date.data.acquisition = Datefix.Conversion.Function(name.data.acquisition,9)
date.data.acquisition = Datefix.Function(date.data.acquisition,9)
date.data.main = Datefix.Conversion.Function(name.data.main,9)
date.data.main = Datefix.Function(date.data.main,9)

## Run QC Analysis
qc.data.pretrain = QC.Pretrain.Function(date.data.pretrain)

qc.data.acquisition = QC.Acq.Function(date.data.acquisition)

qc.data.main = QC.Main.Function(date.data.main)

## Fix Latency Decimals ##
qc.data.acq.latfix = LatFix.Decimal.MainAcq.Function(qc.data.acquisition,c(123:226))
qc.data.main.latfix = LatFix.Decimal.MainAcq.Function(qc.data.main,c(123:226))

## Split Main Data ##
qc.data.mainprobe = qc.data.main.latfix[ ,c(1:18,123:226)]
qc.data.mainvigilance = qc.data.main.latfix[ ,c(1:14,19:122)]
## Run LatFix ##

qc.data.acquisition.lat = LatFix.Acq.Function(qc.data.acq.latfix,3)
qc.data.mainprobe.lat = LatFix.Probe.Function(qc.data.mainprobe,3)

## Aggregate Probe Data ##
agg.list = list(qc.data.mainprobe.lat$Database,qc.data.mainprobe.lat$AnimalID,qc.data.mainprobe.lat$TestSite, qc.data.mainprobe.lat$Mouse.Strain, qc.data.mainprobe.lat$Genotype, qc.data.mainprobe.lat$Sex,qc.data.mainprobe.lat$Age.Months, qc.data.mainprobe.lat$Stimulus.Length)
qc.data.mainprobe.agg = aggregate(qc.data.mainprobe.lat,by=agg.list, FUN=mean, na.rm=TRUE)
qc.data.mainprobe.agg[ ,9:16] = qc.data.mainprobe.agg[ ,1:8]
qc.data.mainprobe.agg[ ,1:8] = NULL


## Fix Vigilance Import Issues) ##
qc.data.mainvigilance.fixed = Vigilance.FixRaw.Function(qc.data.mainvigilance)

## Vigilance Calculation - Three Methods ##
qc.data.vigilance.blocked = Vigilance.Calc.Block.Function(qc.data.mainvigilance.fixed,10)
qc.data.vigilance.averaged = Vigilance.Calc.Average.Function(qc.data.mainvigilance.fixed,10)
qc.data.vigilance.last = Vigilance.Calc.Last.Function(qc.data.mainvigilance.fixed,10)


## Aggregate Final Data ## 
agg.list.blocked = list(qc.data.vigilance.blocked$Database,qc.data.vigilance.blocked$AnimalID,qc.data.vigilance.blocked$TestSite, qc.data.vigilance.blocked$Mouse.Strain, qc.data.vigilance.blocked$Genotype, qc.data.vigilance.blocked$Sex,qc.data.vigilance.blocked$Age.Months, qc.data.vigilance.blocked$Stimulus.Length)
agg.list.averaged = list(qc.data.vigilance.averaged$Database,qc.data.vigilance.averaged$AnimalID,qc.data.vigilance.averaged$TestSite, qc.data.vigilance.averaged$Mouse.Strain, qc.data.vigilance.averaged$Genotype, qc.data.vigilance.averaged$Sex,qc.data.vigilance.averaged$Age.Months, qc.data.vigilance.averaged$Stimulus.Length)
agg.list.last = list(qc.data.vigilance.last$Database,qc.data.vigilance.last$AnimalID,qc.data.vigilance.last$TestSite, qc.data.vigilance.last$Mouse.Strain, qc.data.vigilance.last$Genotype, qc.data.vigilance.last$Sex,qc.data.vigilance.last$Age.Months, qc.data.vigilance.last$Stimulus.Length)

qc.data.vigilance.blocked.agg = aggregate(qc.data.vigilance.blocked,by=agg.list.blocked, FUN=mean, na.rm=TRUE)
qc.data.vigilance.blocked.agg[ ,9:16] = qc.data.vigilance.blocked.agg[ ,1:8]
qc.data.vigilance.blocked.agg[ ,1:8] = NULL

qc.data.vigilance.averaged.agg = aggregate(qc.data.vigilance.averaged,by=agg.list.averaged, FUN=mean, na.rm=TRUE)
qc.data.vigilance.averaged.agg[ ,9:16] = qc.data.vigilance.averaged.agg[ ,1:8]
qc.data.vigilance.averaged.agg[ ,1:8] = NULL

qc.data.vigilance.last.agg = aggregate(qc.data.vigilance.last,by=agg.list.last, FUN=mean, na.rm=TRUE)
qc.data.vigilance.last.agg[ ,9:16] = qc.data.vigilance.last.agg[ ,1:8]
qc.data.vigilance.last.agg[ ,1:8] = NULL
## Remove Non-Completers ##

qc.idlist.final = as.vector(unique(as.character(qc.data.mainprobe.lat$AnimalID)))
qc.data.acquisition.final = qc.data.acquisition.lat[qc.data.acquisition.lat$AnimalID %in% qc.idlist.final, ]
qc.data.pretrain.final = qc.data.pretrain[qc.data.pretrain$AnimalID %in% qc.idlist.final, ]

## Sum & Aggregate Acquisition Sessions ##
qc.data.acq.summed = qc.data.acquisition.final[ ,c(1:8,10)]
qc.data.acq.summed[ ,9] = 1
qc.data.acq.agg = aggregate(Day ~ AnimalID + TestSite + Mouse.Strain + Genotype + Sex + Age.Months + Stimulus.Length, FUN=sum, na.rm=TRUE, data=qc.data.acq.summed)

## Aggregated Pretraining Data ##
qc.data.pretrain.final = qc.data.pretrain.final[ ,c(2:10)]
qc.data.pretrain.final[ ,9] = 1
qc.data.pretrain.final[ ,8] = NULL
qc.data.pretrain.agg = aggregate(Day ~ AnimalID + TestSite + Mouse.Strain + Genotype + Sex + Age.Months + Task, FUN=sum, na.rm=TRUE, data=qc.data.pretrain.final)


## Save Raw Data Files ##
write.csv(qc.data.pretrain.agg, "Weston 5CSRTT Pretrain QC Oct 12 2017 Updated.csv")

write.csv(qc.data.acquisition.final, "Weston 5CSRTT Acquisition QC Oct 12 2017 Updated.csv")
write.csv(qc.data.acq.agg, "Weston 5CSRTT Acquisition Aggregated QC Oct 12 2017 Updated.csv")

write.csv(qc.data.mainprobe.lat, "Weston 5CSRTT Probe QC Oct 12 2017 Updated.csv")
write.csv(qc.data.mainprobe.agg, "Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv")

write.csv(qc.data.vigilance.averaged.agg, "Weston 5CSRTT Probe Aggregated Vigilance AVERAGED METHOD QC Oct 12 2017 Updated.csv")
write.csv(qc.data.vigilance.last.agg, "Weston 5CSRTT Probe Aggregated Vigilance LAST METHOD QC Oct 12 2017 Updated.csv")
write.csv(qc.data.vigilance.blocked.agg, "Weston 5CSRTT Probe Aggregated Vigilance BLOCK BY BLOCK METHOD QC Oct 12 2017 Updated.csv")
