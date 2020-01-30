require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(gtools)
require(LSD)

# =============== Define the main path and read he plate map ============================== 
here_path = "/Users/anaritabrochado/Documents/TypasPC/TypasLab/LabWork/Collaborations/Pneumo/Triclo/Uploaded/"
source(paste0(here_path,"Aux_functions.R",collapse=NULL))
#save.image("~/Documents/Typas Lab/LabWork/Collaborations/Pneumo/Triclo/Prestwick/R-scripts/Prestwick_environment.RData")
#load(paste0(here_path,"R-scripts/Prestwick_environment.RData"))
here_path = "/Users/anaritabrochado/Documents/TypasPC/TypasLab/LabWork/Collaborations/Pneumo/Triclo/Uploaded/"

Load_dir = paste0(here_path)
Out_dir = paste0(here_path)

Map = read.table(file=paste0(Load_dir,"PrestwickMap.txt"),header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=T)
Map = as.data.frame(cbind(row.names(Map),Map))
names(Map)[1:2] = c("well_id","Chemical_Name")

PrestwickData = read.table(file=paste0(Load_dir,"InputData_Prestwick.txt"),header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=T)

Plates = unique(as.character(PrestwickData$plate_id))

#Get all the plate names
Plates_OD = Plates[grep("OD",Plates)]
Plates_Lux = Plates[grep("Lux",Plates)]

Plates_OD_data = list()
Plates_Lux_data = list()

for(p in 1:length(Plates_OD))
{
  plate_id_OD = Plates_OD[p]
  dataset_OD = PrestwickData[grep(plate_id_OD,PrestwickData$plate_id),]
  dataset_OD = dataset_OD[-1]
  if(length(grep(T,is.na(dataset_OD[1,])))>0)
  {dataset_OD = dataset_OD[grep(F,is.na(dataset_OD[1,]))]}
  Plates_OD_data[[length(Plates_OD_data)+1]] = dataset_OD
  
  plate_id_Lux = Plates_Lux[p]
  dataset_Lux = PrestwickData[grep(plate_id_Lux,PrestwickData$plate_id),]
  dataset_Lux = dataset_Lux[-1]
  if(length(grep(T,is.na(dataset_Lux[1,])))>0)
  {dataset_Lux = dataset_Lux[grep(F,is.na(dataset_Lux[1,]))]}
  Plates_Lux_data[[length(Plates_Lux_data)+1]] = dataset_Lux
}
names(Plates_OD_data) = Plates_OD
names(Plates_Lux_data) = Plates_Lux

#============== Plot growth curves and lux raw data ================
pdf(paste0(Out_dir,"All Plates OD",".pdf"))
#for(p in 1:length(Plates_OD_data))
for(p in 1:6)  
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  for(q in 1:4)
  {
    data_OD = quadr(q,dataset_OD)
    par(mfrow=c(8,12),mar=c(1,1,1,1))
    for(w in 2:length(data_OD))
    {
      plot(data_OD[[1]],data_OD[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",ylim=c(0,0.5),xlim=c(0,7))
      text(x=1,y=0.45,cex=1,labels=names(data_OD)[w],adj=0)
    }
    title(paste0(plate_id," Quadrant ",q),outer=T,line=-1)
  }
}
for(p in 7:8)
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  
  data_OD = dataset_OD
  par(mfrow=c(8,12),mar=c(1,1,1,1))
  for(w in 2:length(data_OD))
  {
    plot(data_OD[[1]],data_OD[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",ylim=c(0,0.5),xlim=c(0,7))
    text(x=1,y=0.45,cex=1,labels=names(data_OD)[w],adj=0)
  }
  title(paste0(plate_id),outer=T,line=-1)

}
dev.off()

pdf(paste0(Out_dir,"All Plates Lux",".pdf"))
#for(p in 1:length(Plates_Lux_data))
for(p in 1:6)
{
  plate_id = names(Plates_Lux_data)[p]
  dataset_Lux = Plates_Lux_data[[p]]
  for(q in 1:4)
  {
    data_Lux = quadr(q,dataset_Lux)
    #chemicals = get_feature(names(data_OD[2:length(data_OD)]),feature = "Chemical_Name",Plate_database = Map)
    
    par(mfrow=c(8,12),mar=c(1,1,1,1))
    for(w in 2:length(data_Lux))
    {
      plot(data_Lux[[1]],data_Lux[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",xlim=c(0,7),ylim=c(0,500),col="red3")
      text(x=1,y=450,cex=1,labels=names(data_Lux)[w],adj=0)
    }
    title(paste0(plate_id," Quadrant ",q),outer=T,line=-1)
  }
}
for(p in 7:8)
{
  plate_id = names(Plates_Lux_data)[p]
  dataset_Lux = Plates_Lux_data[[p]]
  data_Lux = dataset_Lux
  
  par(mfrow=c(8,12),mar=c(1,1,1,1))
  for(w in 2:length(data_Lux))
  {
    plot(data_Lux[[1]],data_Lux[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",xlim=c(0,7),ylim=c(0,500),col="red3")
    text(x=1,y=450,cex=1,labels=names(data_Lux)[w],adj=0)
  }
  title(paste0(plate_id),outer=T,line=-1)
}
dev.off()

# ========== Replicate correlation ==========
pdf(paste0(Out_dir,"Prestwick Analysis_updated",".pdf"),useDingbats=FALSE)
time_cutoff = 4.8
for(p in 1:length(Plates_OD_data))
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  
  time_point = grep(T,time_cutoff<dataset_OD[[1]])[1] - 1
  #plate_data = as.data.frame(t(dataset_OD[time_point,]))
  
  plate_data = as.data.frame((Trape_rule(dataset_OD[1:time_point,])))
  names(plate_data) = plate_id
  
  if(p==7 || p==8)
  {
    plate_data_format = rep(NA,length(All_Plates_OD[[1]]))
    plate_data_format[match(row.names(plate_data),row.names(All_Plates_OD))]=plate_data[[1]]
    plate_data_format = as.data.frame(plate_data_format)
    names(plate_data_format) = plate_id
    row.names(plate_data_format) = row.names(All_Plates_OD)
    plate_data=plate_data_format
  }
  
  if(p==1)
  {All_Plates_OD = plate_data} else
  {All_Plates_OD = as.data.frame(cbind(All_Plates_OD,plate_data))}
}
for(p in 1:length(Plates_Lux_data))
{
  plate_id = names(Plates_Lux_data)[p]
  dataset_Lux = Plates_Lux_data[[p]]
  
  time_point = grep(T,time_cutoff<dataset_Lux[[1]])[1] - 1
  #plate_data = as.data.frame(t(dataset_Lux[time_point,]))
  plate_data = as.data.frame((Trape_rule(dataset_Lux[1:time_point,])))
  names(plate_data) = plate_id
  
  if(p==7 || p==8)
  {
    plate_data_format = rep(NA,length(All_Plates_OD[[1]]))
    plate_data_format[match(row.names(plate_data),row.names(All_Plates_OD))]=plate_data[[1]]
    plate_data_format = as.data.frame(plate_data_format)
    names(plate_data_format) = plate_id
    row.names(plate_data_format) = row.names(All_Plates_OD)
    plate_data=plate_data_format
  }
  
  
  if(p==1)
  {All_Plates_Lux = plate_data} else
  {All_Plates_Lux = as.data.frame(cbind(All_Plates_Lux,plate_data))}
}

pairs(All_Plates_OD,lower.panel = panel.cor,cex=0.3,col="grey",pch=19,cex.cor=1.5,main="AUC OD")
pairs(All_Plates_Lux,lower.panel = panel.cor,cex=0.3,pch=19,cex.cor=1.5,col="red3",main="AUC Lux")

cols = c("darkgreen","darkgreen","red3","red3","blue","blue","orange","orange")
control_wells = as.character(Map[grep("Control",Map[[2]]),1])
par(mfrow=c(1,1),mar=c(10,4,4,4))
boxplot(All_Plates_OD,ylim=c(0,0.4),ylab="AUC OD",las=2, border=cols,pch=19,cex=0.5,
        main=paste0("AUC for OD, ",time_cutoff,"h"))
boxplot(All_Plates_Lux,ylim=c(0,1000),ylab="AUC Lux",las=2,border=cols,pch=19,cex=0.5,
        main=paste0("AUC for Lux, ",time_cutoff,"h"))

boxplot(All_Plates_Lux/All_Plates_OD,ylim=c(0,5000),ylab="AUC Lux",border=cols,las=2,pch=19,cex=0.5,
        main=paste0("AUC for Lux/AUC OD, ",time_cutoff,"h"))

All_Plates_OD[All_Plates_OD<0] = 0
All_Plates_Lux[All_Plates_Lux<0] = 0
#============== Plot lux against growth to get hits ================
lux_cutoff = 0.3
growth_cutoff = 0.7
par(mfrow=c(2,2),mar=c(4,4,4,4))
for(p in 1:length(All_Plates_OD))
{
  plate_id = names(All_Plates_OD)[p]
  p_OD = All_Plates_OD[[p]]
  p_Lux = All_Plates_Lux[[p]]
  plib = strsplit(plate_id,split = "_")[[1]][1]
  
  plot(p_OD,p_Lux,cex=0.5,pch=19,col="grey",ylab="Lux",xlab="OD",main=plate_id)
  fit=lm(p_Lux~p_OD)
  abline(fit)
  points(c(growth_cutoff*median(p_OD,na.rm=T),growth_cutoff*median(p_OD,na.rm=T)),c(0,1000),type="l")
  points(c(0,1),c(lux_cutoff*median(p_Lux,na.rm=T),lux_cutoff*median(p_Lux,na.rm=T)),type="l")
  
  resids = residuals(fit)
  resids_cutoff = quantile(resids,probs = 0.075)
  
  c1 = as.numeric(names(residuals(fit))[grep(T,residuals(fit)<resids_cutoff)])
  c2 = grep(T,p_OD > growth_cutoff*median(p_OD,na.rm=T))
  c3 = grep(T,p_Lux < lux_cutoff*median(p_Lux,na.rm=T))
  
  c1_c2 = c1[c1 %in% c2]
  c1_c2_c3 = c1_c2[c1_c2 %in% c3]
  points(p_OD[c1_c2_c3],p_Lux[c1_c2_c3],pch=19,cex=0.7,col="red3")
  #points(p_OD[c1],p_Lux[c1],pch=19,cex=0.5,col="green")
  
  low_hits = c1_c2_c3
  high_hits = c() #dummy
  
  low_hits_wells = row.names(All_Plates_OD)[(low_hits)]
  high_hits_wells = row.names(All_Plates_OD)[(high_hits)]
  
  
  Plate_hits = as.data.frame(cbind(paste0(plib,"_",c(low_hits_wells,high_hits_wells)),
                                   c(rep("Low lux",length(low_hits)),rep("High lux",length(high_hits)))))
  Chemicals = get_feature(as.character(Plate_hits[[1]]),feature = "Chemical_Name",Plate_database = Map)
  Therap_effect = get_feature(as.character(Plate_hits[[1]]),feature = "Therapeutic_Class",Plate_database = Map)
  Plate_hits = as.data.frame(cbind(Plate_hits,Chemicals,Therap_effect,rep(plate_id,length(Plate_hits[[1]]))))
  names(Plate_hits)[c(1,2,5)] = c("well_id","Lux signal","Plate_id")
  
  if(p==1)
  {
    All_hits = Plate_hits
  } else
  {
    All_hits = as.data.frame(rbind(All_hits,Plate_hits))
  }
}

nr_replicates = table(as.character(All_hits[[1]]))
bad_wells = names(nr_replicates)[nr_replicates<2]
All_hits = All_hits[-match(bad_wells,as.character(All_hits[[1]])),]
All_hits = All_hits[-5]
All_hits = All_hits[grep(F,duplicated(All_hits)),]

for(h in 1:length(All_hits[[1]]))
{
  hit_well = as.character(All_hits[[1]][h])
  well = strsplit(hit_well,split = "_")[[1]][2]
  plate = strsplit(hit_well,split = "_")[[1]][1]
  
  OD_data = as.numeric(All_Plates_OD[match(well,row.names(All_Plates_OD)),grep(plate,names(All_Plates_OD))])
  lux_data = as.numeric(All_Plates_Lux[match(well,row.names(All_Plates_Lux)),grep(plate,names(All_Plates_Lux))])
  
  well_data = as.data.frame(t(c(OD_data,lux_data)))
  names(well_data) = c("OD AUC rep 1","OD AUC rep 2","Lux AUC rep 1","Lux AUC rep 2")
  
  if(h==1)
  {All_hits_data = well_data} else
  {All_hits_data = as.data.frame(rbind(All_hits_data,well_data))}
}

All_hits = as.data.frame(cbind(All_hits,All_hits_data))

file_id = paste0(Out_dir,"Hits_updated.txt")
write.table(All_hits,file_id,quote = F,sep = "\t",row.names=F)
dev.off()

#get All_AUC

compiled_row_names = paste0(rep(c("P1_","P2_","P3_","P4_"),each=384),row.names(All_Plates_Lux))

compiled_Plates_Lux = as.data.frame(cbind(as.vector(as.matrix(All_Plates_Lux[seq(1,length(All_Plates_Lux),by=2)])),
                                          as.vector(as.matrix(All_Plates_Lux[seq(2,length(All_Plates_Lux),by=2)]))))
compiled_Plates_OD = as.data.frame(cbind(as.vector(as.matrix(All_Plates_OD[seq(1,length(All_Plates_OD),by=2)])),
                                         as.vector(as.matrix(All_Plates_OD[seq(2,length(All_Plates_OD),by=2)]))))
All_AUC_out = as.data.frame(cbind(compiled_Plates_Lux,compiled_Plates_OD))
row.names(All_AUC_out) = compiled_row_names

All_AUC_out = All_AUC_out[grep(F,is.na(All_AUC_out[1])),]
names(All_AUC_out) = c("Lux - rep 1","Lux - rep 2","OD - rep 1","OD - rep 2")

All_AUC_out = as.data.frame(cbind(get_feature(row.names(All_AUC_out),feature = "Chemical_Name",Plate_database = Map),
                                  get_feature(row.names(All_AUC_out),feature = "Therapeutic_Class",Plate_database = Map),
                                  get_feature(row.names(All_AUC_out),feature = "Therapeutic_Effect",Plate_database = Map),
                                  All_AUC_out
))
names(All_AUC_out) = c("Chemical_Name","Therapeutic_Class","Therapeutic_Effect","Lux - rep 1","Lux - rep 2","OD - rep 1","OD - rep 2")

file_id = paste0(Out_dir,"All_AUC.txt")
write.table(All_AUC_out,file_id,quote = F,sep = "\t",row.names=T)
# ========= get interesting compunds in individual plots =========
pdf(paste0(Out_dir,"Individual hits.pdf"))
unique_Hits = as.character(unique(All_hits[[1]]))
par(mfrow=c(1,1),mar=c(6,6,6,6),lwd=2)
for(h in 1:length(unique_Hits))
{
  hit = unique_Hits[h]
  plate = strsplit(hit,split="_")[[1]][1]
  well = strsplit(hit,split="_")[[1]][2]
  compound = get_feature(hit,feature = "Chemical_Name",Plate_database = Map)
  
  plate_ids_OD = names(Plates_OD_data)[grep(plate,names(Plates_OD_data))]
  plate_ids_Lux = names(Plates_Lux_data)[grep(plate,names(Plates_Lux_data))]
  
  for(p in 1:length(plate_ids_OD))
  {
    plate_id_OD = plate_ids_OD[p]
    p_OD = Plates_OD_data[[match(plate_id_OD,names(Plates_OD_data))]]
    
    if(p==1)
    {
      plot(p_OD[[1]],p_OD[[match(well,names(p_OD))]],ylab="",xlab="Time",ylim=c(0,0.4),type="l",main=paste0(hit, ", ",compound),
           axes=F,frame=T,col="#8CA336",xlim=c(0,7))
      axis(2,seq(0,0.4,by=0.1),las=2)
      mtext("OD",side=2,line=2.5,las=2)
      axis(1,seq(0,7,by=1),las=2)
    } else
    {points(p_OD[[1]],p_OD[[match(well,names(p_OD))]],type="l",lty=2,col="#8CA336")}
  }
  
  for(p in 1:length(plate_ids_Lux))
  {
    plate_id_Lux = plate_ids_Lux[p]
    p_Lux = Plates_Lux_data[[match(plate_id_Lux,names(Plates_Lux_data))]]
    
    par(new=T)
    if(p==1)
    {
      plot(p_Lux[[1]],p_Lux[[match(well,names(p_Lux))]],type="l",ylim=c(0,500),axes = F,ylab="",xlab="",col="orange",xlim=c(0,7))
      axis(4,seq(0,500,by=100),las=2)
      mtext("Lux",side=4,line=2.5,las=2)
    } else
    {plot(p_Lux[[1]],p_Lux[[match(well,names(p_Lux))]],type="l",lty=2,ylim=c(0,500),axes = F,ylab="",xlab="",col="orange",xlim=c(0,7))}
  }
  legend("topleft",legend=c("rep 1","rep 2"),lty=c(1,2),bty="n")  
  
}
dev.off()

# ========= All prestwick drugs =========
pdf(paste0(Out_dir,"All prestwick drugs.pdf"))
All_wells = as.character(Map[[1]])
par(mfrow=c(2,2),mar=c(4,4,4,4),lwd=2)
for(h in 1:length(All_wells))
#for(h in 1:16)
{
  hit = All_wells[h]
  plate = strsplit(hit,split="_")[[1]][1]
  well = strsplit(hit,split="_")[[1]][2]
  compound = get_feature(hit,feature = "Chemical_Name",Plate_database = Map)
  
  plate_ids_OD = names(Plates_OD_data)[grep(plate,names(Plates_OD_data))]
  plate_ids_Lux = names(Plates_Lux_data)[grep(plate,names(Plates_Lux_data))]
  
  for(p in 1:length(plate_ids_OD))
  {
    plate_id_OD = plate_ids_OD[p]
    p_OD = Plates_OD_data[[match(plate_id_OD,names(Plates_OD_data))]]
    
    if(p==1)
    {
      plot(p_OD[[1]],p_OD[[match(well,names(p_OD))]],ylab="",xlab="Time",ylim=c(0,0.4),type="l",main=paste0(hit, ", ",compound),
           axes=F,frame=T,col="#8CA336",xlim=c(0,7))
      axis(2,seq(0,0.4,by=0.1),las=2)
      mtext("OD",side=2,line=2.5,las=2)
      axis(1,seq(0,7,by=1),las=2)
    } else
    {points(p_OD[[1]],p_OD[[match(well,names(p_OD))]],type="l",lty=2,col="#8CA336")}
  }
  
  for(p in 1:length(plate_ids_Lux))
  {
    plate_id_Lux = plate_ids_Lux[p]
    p_Lux = Plates_Lux_data[[match(plate_id_Lux,names(Plates_Lux_data))]]
    
    par(new=T)
    if(p==1)
    {
      plot(p_Lux[[1]],p_Lux[[match(well,names(p_Lux))]],type="l",ylim=c(0,500),axes = F,ylab="",xlab="",col="orange",xlim=c(0,7))
      axis(4,seq(0,500,by=100),las=2)
      mtext("Lux",side=4,line=2.5,las=2)
    } else
    {plot(p_Lux[[1]],p_Lux[[match(well,names(p_Lux))]],type="l",lty=2,ylim=c(0,500),axes = F,ylab="",xlab="",col="orange",xlim=c(0,7))}
  }
  legend("topleft",legend=c("rep 1","rep 2"),lty=c(1,2),bty="n")  
  
}
dev.off()

