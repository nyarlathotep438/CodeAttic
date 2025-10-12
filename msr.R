#需要处理的文件放在哪里，请在下方引号内填写
ADDRESS = "D:/Ready"
setwd(ADDRESS) 
library("msd")

FILES = list.files(ADDRESS)
FILE_NUM = length(FILES)
for (j in 1:FILE_NUM) {
  file1 = FILES[j]
  
  #msd 计算核心部分
  dat1=read.csv(file1,header=T) 
  data1=as.matrix(dat1)
  data1 <- apply(data1,2,as.numeric)
  msd1 = msd(data = data1,misfit = TRUE)           
  
  #GAMA值计算
  item1 = msd1[["item_measures"]]
  person1 = msd1[["person_measures"]]
  n_max = length(msd1[["person_measures"]])
  m_max = length(msd1[["item_measures"]])
  G = list()                                
  for (n in 1:n_max) {                        
    for (m in 1:m_max) {                        
      t = m + m_max*(n-1)                       
      gama = person1[n]-item1[m]
      gama ->G[t]
    }
  }
  
  #probabililties
  f1 <- function(r){                        
    1/(1+exp(-r))
  }
  PR = list()  
  lr = length(msd1[["thresholds"]])+1      
  GG = as.numeric(G)
  for (t in 1:(m*n)) {                     
    for (rr in 1:lr) {
      u = rr + lr*(t-1)
      if(rr == 1){
        pr = f1(-(GG[t]-msd1[["thresholds"]][rr]))
        pr -> PR[u]
      }else if(rr == lr){
        pr = 1-f1(-(GG[t]-msd1[["thresholds"]][rr-1]))
        pr -> PR[u]
      }else {
        pr = f1(-(GG[t]-msd1[["thresholds"]][rr]))-f1(-(GG[t]-msd1[["thresholds"]][rr-1]))
        pr -> PR[u]
      }
    }
  }
  
  
  PRR = list(NULL)
  length(PRR) <- t
  PRPR = NULL
  for (i in 1:lr) {
    PRPR[[i]] <- PRR
  }
  
  #sorting 按顺序将PR值分组PRR
  for (u in 1:(t*lr)) {                       
    for (t in 1:(m*n)) {                      
      for (i in 1:lr) {
        if(u==(i+lr*(t-1))){
          PR[u] -> PRPR[[i]][t]
        }
      }
    }
  }
  
  #移除影响作图的NA值
  NA4 = list()
  b = 0
  for (t in 1:(m*n)) {
    if(is.na(G[t])){
      b = b+1
      t -> NA4[b]
    }
  }
  
  G1 = G
  PRPR2 = PRPR
  if(length(NA4) != 0){
    for (c in b:1) {
      d=as.numeric(NA4[c])
      G1[d] = NULL
    }
    for (i in 1:lr) {
      for (c in b:1) {
        d=as.numeric(NA4[c])
        PRPR2[[i]][d] = NULL
      }
    }
  }
  
  #对文件进行操作
  file_sample = substring(file1,0,18)
  x = paste(file_sample,"...result",sep = "")
  x1 = paste(x,".tiff",sep = "")
  x2_1 =paste(x,"$item_measures",".csv",sep = "")
  x2_2 =paste(x,"$person_measures",".csv",sep = "")
  x2_3 =paste(x,"$thresholds",".csv",sep = "")
  x2_4 =paste(x,"$item_std_errors",".csv",sep = "")
  x2_5 =paste(x,"$person_std_errors",".csv",sep = "")
  x2_6 =paste(x,"$threshold_std_errors",".csv",sep = "")
  x2_7 =paste(x,"$item_reliability",".csv",sep = "")
  x2_8 =paste(x,"$person_reliability",".csv",sep = "")
  x2_9 =paste(x,"$infit_items",".csv",sep = "")
  x2_10 =paste(x,"$infit_persons",".csv",sep = "")
  x2_11 =paste(x,"$outfit_items",".csv",sep = "")
  x2_12 =paste(x,"$outfit_persons",".csv",sep = "")
  dir.create(x)
  result_adress_tiff = paste(ADDRESS,x,x1,sep = "/")
  result_adress_csv1 = paste(ADDRESS,x,x2_1,sep = "/")
  result_adress_csv2 = paste(ADDRESS,x,x2_2,sep = "/")
  result_adress_csv3 = paste(ADDRESS,x,x2_3,sep = "/")
  result_adress_csv4 = paste(ADDRESS,x,x2_4,sep = "/")
  result_adress_csv5 = paste(ADDRESS,x,x2_5,sep = "/")
  result_adress_csv6 = paste(ADDRESS,x,x2_6,sep = "/")
  result_adress_csv7 = paste(ADDRESS,x,x2_7,sep = "/")
  result_adress_csv8 = paste(ADDRESS,x,x2_8,sep = "/")
  result_adress_csv9 = paste(ADDRESS,x,x2_9,sep = "/")
  result_adress_csv10 = paste(ADDRESS,x,x2_10,sep = "/")
  result_adress_csv11 = paste(ADDRESS,x,x2_11,sep = "/")
  result_adress_csv12 = paste(ADDRESS,x,x2_12,sep = "/")
  write.csv(msd1[["item_measures"]],result_adress_csv1)
  write.csv(msd1[["person_measures"]],result_adress_csv2)
  write.csv(msd1[["thresholds"]],result_adress_csv3)
  write.csv(msd1[["item_std_errors"]],result_adress_csv4)
  write.csv(msd1[["person_std_errors"]],result_adress_csv5)
  write.csv(msd1[["threshold_std_errors"]],result_adress_csv6)
  write.csv(msd1[["item_reliability"]],result_adress_csv7)
  write.csv(msd1[["person_reliability"]],result_adress_csv8)
  write.csv(msd1[["infit_items"]],result_adress_csv9)
  write.csv(msd1[["infit_persons"]],result_adress_csv10)
  write.csv(msd1[["outfit_items"]],result_adress_csv11)
  write.csv(msd1[["outfit_persons"]],result_adress_csv12)
  
  #绘图
  tiff(filename = result_adress_tiff)
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
  plot(G,PRPR[[1]],pch=".",col="#ed1299",xlab="Gamma",ylab="probabililty")
  lines(smooth.spline(G1,PRPR2[[1]]),col="#ed1299",lwd="2")
  for (i in (2:lr)) {
    points(G,PRPR[[i]],pch=".",col=cb_palette[i])
    lines(smooth.spline(G1,PRPR2[[i]]),col=cb_palette[i],lwd="2")
  }
  legend("right",legend=c(1:lr),col=cb_palette,lty=1,lwd=2)
  dev.off()
}
