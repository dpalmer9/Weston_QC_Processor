raw.data = file.choose()

Vigilance.Calc.Block.Function = function(dataset,binsize){
  acc.start = 15
  omission.start = 65
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
        if(isTRUE(((curr.om > start.om) | ((curr.om == 100) & (start.om == 100))))){
          om.bin = om.bin + 1
        }else if(isTRUE((curr.acc > start.acc) | ((curr.acc == 100) & (start.acc == 100)))){
          acc.bin = acc.bin + 1
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
  final.data = cbind(dataset[ ,c(1:14)],new.data)
  return(final.data)
}

vig.data = Vigilance.Calc.Block.Function(raw.data,10)

write.csv(vig.data,'VAChT-KD 5-CSRTT Data.csv',row.names=FALSE)
