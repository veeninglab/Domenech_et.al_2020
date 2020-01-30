require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(gtools)
require(LSD)

# =============== Define the main path and read he plate map ============================== 
here_path = "Define!!!!"
source(paste0(here_path,"Aux_functions.R",collapse=NULL))

Load_dir = paste0(here_path)
Out_dir = paste0(here_path)

Map = read.table(file=paste0(Load_dir,"TailoredMap.txt"),header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=T)


TailoredData = read.table(file=paste0(Load_dir,"InputData_Tailored.txt"),header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=T)

Plates = unique(as.character(TailoredData$plate_id))

#Get all the plate names
Plates_OD = Plates[grep("OD",Plates)]
Plates_Lux = Plates[grep("Lux",Plates)]

Plates_OD_data = list()
Plates_Lux_data = list()

for(p in 1:length(Plates_OD))
{
  plate_id_OD = Plates_OD[p]
  dataset_OD = TailoredData[grep(plate_id_OD,TailoredData$plate_id),]
  dataset_OD = dataset_OD[-1]
  if(length(grep(T,is.na(dataset_OD[1,])))>0)
  {dataset_OD = dataset_OD[grep(F,is.na(dataset_OD[1,]))]}
  Plates_OD_data[[length(Plates_OD_data)+1]] = dataset_OD
  
  plate_id_Lux = Plates_Lux[p]
  dataset_Lux = TailoredData[grep(plate_id_Lux,TailoredData$plate_id),]
  dataset_Lux = dataset_Lux[-1]
  if(length(grep(T,is.na(dataset_Lux[1,])))>0)
  {dataset_Lux = dataset_Lux[grep(F,is.na(dataset_Lux[1,]))]}
  Plates_Lux_data[[length(Plates_Lux_data)+1]] = dataset_Lux
}
names(Plates_OD_data) = Plates_OD
names(Plates_Lux_data) = Plates_Lux


#============== Plot growth curves and lux raw data ================
pdf(paste0(Out_dir,"All Plates OD",".pdf"))
for(p in 1:length(Plates_OD_data))
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  for(q in 1:4)
  {
    data_OD = quadr(q,dataset_OD)
    par(mfrow=c(8,12),mar=c(1,1,1,1))
    for(w in 2:length(data_OD))
    {
      plot(data_OD[[1]],data_OD[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",ylim=c(0,0.5))
      text(x=1,y=0.45,cex=1,labels=names(data_OD)[w],adj=0)
    }
    title(paste0(plate_id," Quadrant ",q),outer=T,line=-1)
  }
}
dev.off()

pdf(paste0(Out_dir,"All Plates Lux",".pdf"))
for(p in 1:length(Plates_Lux_data))
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
      plot(data_Lux[[1]],data_Lux[[w]],type="b",pch=19,cex=0.5,axes=F,frame=T,xlab="",ylab="",ylim=c(0,500),col="red3")
      text(x=1,y=450,cex=1,labels=names(data_Lux)[w],adj=0)
    }
    title(paste0(plate_id," Quadrant ",q),outer=T,line=-1)
  }
}
dev.off()

rm(data_Lux,data_OD)

# ========== Plot wells to decide on time point ==========

control_well = "N6"
time_cutoff = 5
cols = c("black","red3","orange","blue")
for(p in 1:length(Plates_OD_data))
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  
  if(p==1)
  {
    plot(dataset_OD[[1]],dataset_OD[[match(control_well,names(dataset_OD))]],type="b",ylab="OD",xlab="Time",main=control_well,
         pch=19,col=cols[p],ylim=c(0,0.5))
  } else
  {points(dataset_OD[[1]],dataset_OD[[match(control_well,names(dataset_OD))]],type="b",pch=19,col=cols[p])}
  
}
points(c(time_cutoff,time_cutoff),c(-1,1),type="l")
legend("topleft",legend=names(Plates_OD_data),pch=19,col=cols)

cols = c("black","red3","orange","blue")
for(p in 1:length(Plates_Lux_data))
{
  plate_id = names(Plates_Lux_data)[p]
  dataset_Lux = Plates_Lux_data[[p]]
      
  if(p==1)
  {
    plot(dataset_Lux[[1]],dataset_Lux[[match(control_well,names(dataset_Lux))]],type="b",ylab="OD",xlab="Time",main=control_well,
         pch=19,col=cols[p],ylim=c(0,500))
  } else
  {points(dataset_Lux[[1]],dataset_Lux[[match(control_well,names(dataset_Lux))]],type="b",pch=19,col=cols[p])}
}
points(c(time_cutoff,time_cutoff),c(-200,1000),type="l")
legend("topleft",legend=names(Plates_Lux_data),pch=19,col=cols)

# ========== Replicate correlation ==========
pdf(paste0(Out_dir,"Tailored analysis_updated.pdf"),useDingbats=FALSE)
time_cutoff = 5
for(p in 1:length(Plates_OD_data))
{
  plate_id = names(Plates_OD_data)[p]
  dataset_OD = Plates_OD_data[[p]]
  
  time_point = grep(T,time_cutoff<dataset_OD[[1]])[1] - 1
  #plate_data = as.data.frame(t(dataset_OD[time_point,]))
  plate_data = as.data.frame((Trape_rule(dataset_OD[1:time_point,])))
  names(plate_data) = plate_id
  
  if(p==1)
  {All_Plates_OD = plate_data} else
  {All_Plates_OD = as.data.frame(cbind(All_Plates_OD,plate_data))}
}
for(p in 1:length(Plates_Lux_data))
{
  plate_id = names(Plates_Lux_data)[p]
  dataset_Lux = Plates_Lux_data[[p]]
  
  plate_data = as.data.frame((Trape_rule(dataset_Lux[1:time_point,])))
  names(plate_data) = plate_id
  
  if(p==1)
  {All_Plates_Lux = plate_data} else
  {All_Plates_Lux = as.data.frame(cbind(All_Plates_Lux,plate_data))}
}

pairs(All_Plates_OD,lower.panel = panel.cor,cex=0.3,col="grey",pch=19,cex.cor=1.5,main="AUC OD")
pairs(All_Plates_Lux,lower.panel = panel.cor,cex=0.3,pch=19,cex.cor=1.5,col="red3",main="AUC Lux")

control_wells = as.character(Map[grep("Control",Map[[2]]),1])
par(mfrow=c(1,1),mar=c(10,4,4,4))
cols = c("black","red3","orange","blue")
boxplot(All_Plates_OD[match(control_wells,row.names(All_Plates_OD)),],border=cols,ylim=c(0,0.5),ylab="AUC OD",las=2,
        main=paste0("AUC of OD of control wells ",time_cutoff,"h"))
boxplot(All_Plates_Lux[match(control_wells,row.names(All_Plates_Lux)),],border=cols,ylim=c(0,1000),ylab="AUC Lux",las=2,
        main=paste0("AUC of Lux of control wells ",time_cutoff,"h"))
boxplot(All_Plates_Lux[match(control_wells,row.names(All_Plates_Lux)),]/All_Plates_OD[match(control_wells,row.names(All_Plates_OD)),],
        border=cols,ylim=c(0,2000),ylab="AUC Lux/AUC OD",las=2,
        main=paste0("AUC of Lux/AUC of OD of control wells ",time_cutoff,"h"))

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
  
  fit_high = fit
  var_tresh = sd(p_Lux[match(control_wells,row.names(All_Plates_Lux))])
  fit_high$coefficients[[1]] = fit$coefficients[[1]] + 7*var_tresh
  theor_lux_high = p_OD*fit_high$coefficients[[2]]+fit_high$coefficients[[1]]
  high_hits = grep(T,p_Lux>theor_lux_high)
  points(p_OD[high_hits],p_Lux[high_hits],pch=19,col="orange",cex=0.7)
  
  low_hits_wells = row.names(All_Plates_OD)[(low_hits)]
  high_hits_wells = row.names(All_Plates_OD)[(high_hits)]
  
  Plate_hits = as.data.frame(cbind(c(low_hits_wells,high_hits_wells),
                                   c(rep("Low lux",length(low_hits)),rep("High lux",length(high_hits)))))
  Chemicals = get_feature(as.character(Plate_hits[[1]]),feature = "Drug_conc",Plate_database = Map)
  Plate_hits = as.data.frame(cbind(Plate_hits,Chemicals,rep(plate_id,length(Plate_hits[[1]]))))
  names(Plate_hits)[c(1,2,4)] = c("well_id","Lux signal","Plate_id")
  
  if(p==1)
  {
    All_hits = Plate_hits
  } else
  {
    All_hits = as.data.frame(rbind(All_hits,Plate_hits))
  }
}

All_hits = All_hits[-grep("High",All_hits$`Lux signal`),]

nr_replicates = table(as.character(All_hits[[1]]))
bad_wells = names(nr_replicates)[nr_replicates<2]
All_hits = All_hits[-match(bad_wells,as.character(All_hits[[1]])),]

for(h in 1:length(All_hits[[1]]))
{
  well = as.character(All_hits[[1]][h])
  plate = as.character(All_hits[[4]][h])
  
  OD_data = as.numeric(All_Plates_OD[match(well,row.names(All_Plates_OD)),grep_exact(plate,names(All_Plates_OD))])
  lux_data = as.numeric(All_Plates_Lux[match(well,row.names(All_Plates_Lux)),grep_exact(plate,names(All_Plates_OD))])
  
  well_data = as.data.frame(t(c(OD_data,lux_data)))
  names(well_data) = c("OD AUC","Lux AUC")
  
  if(h==1)
  {All_hits_data = well_data} else
  {All_hits_data = as.data.frame(rbind(All_hits_data,well_data))}
}

All_hits = as.data.frame(cbind(All_hits,All_hits_data))

file_id = paste0(Out_dir,"Hits_updated.txt")
write.table(All_hits,file_id,quote = F,sep = "\t",row.names=F)

par(mfrow=c(1,1),mar=c(10,4,4,4))
All_Low_Lux_drugs = get_feature(as.character(All_hits[grep("Low",All_hits[[2]]),1]),feature="Drug",Plate_database = Map)
barplot(sort(table(All_Low_Lux_drugs)),las=2,col="red3",cex.names=0.8,ylab="Nr wells",main="Competence inhibitors")

#===== Plot all AUC data =======

All_AUC_out = as.data.frame(cbind(All_Plates_Lux,All_Plates_OD,
                                  get_feature(row.names(All_Plates_Lux),feature = "Drug",Plate_database = Map),
                                  get_feature(row.names(All_Plates_Lux),feature = "Conc_full",Plate_database = Map),
                                  get_feature(row.names(All_Plates_Lux),feature = "Conc_quarter",Plate_database = Map)
                                  ))
file_id = paste0(Out_dir,"All_AUC.txt")
write.table(All_AUC_out,file_id,quote = F,sep = "\t",row.names=T)

dev.off()

#====== Plot Hits =======
pdf(paste0(Out_dir,"Hits_updated.pdf"))

base_cols = colorRampPalette(c("red3","white"))(14)
base_cols = base_cols[seq(1,length(base_cols),by=2)]

par(mfrow=c(2,2),mar=c(4,4,4,4))
Drugs = unique(All_Low_Lux_drugs)
for(d in 1:length(Drugs))
{
  drug=Drugs[d]
  #plot(dataset_OD,dataset_Lux,cex=0.5,pch=19,col="grey",ylab="Lux",xlab="OD",main="The positive controls - inhibiting")
  plot(as.vector(as.matrix(All_Plates_OD)),as.vector(as.matrix(All_Plates_Lux)),
       cex=0.5,pch=19,col="grey",ylab="Lux",xlab="OD",main=drug,ylim=c(0,800))
  points(c(growth_cutoff*median(as.matrix(All_Plates_OD)),growth_cutoff*median(as.matrix(All_Plates_OD))),c(-100,1000),type="l") #cutoff growth
  points(c(-10,1),c(lux_cutoff*median(as.matrix(All_Plates_Lux)),lux_cutoff*median(as.matrix(All_Plates_Lux))),type="l") # cutoff lux
  
  wells = as.character(Map[grep(drug,Map[[2]]),1])
  cols = c(base_cols[1:4],base_cols[1:4],base_cols[3:6],base_cols[3:6])
  pchs = c(rep(19,4),rep(17,4),rep(19,4),rep(17,4))
  pchs2 = c(rep(1,4),rep(2,4),rep(1,4),rep(2,4))
  points(as.vector(as.matrix(All_Plates_OD[match(wells,row.names(All_Plates_OD)),])),
         as.vector(as.matrix(All_Plates_Lux[match(wells,row.names(All_Plates_Lux)),])),
         pch=pchs,col=cols,cex=1.5)
  points(as.vector(as.matrix(All_Plates_OD[match(wells,row.names(All_Plates_OD)),])),
         as.vector(as.matrix(All_Plates_Lux[match(wells,row.names(All_Plates_Lux)),])),
         pch=pchs2,col="red3",cex=1.5,lwd=2)
  legend("topleft",legend = c("1-rep1","2-rep1","3-rep1","4-rep1",
                              "1-rep2","2-rep2","3-rep2","4-rep2",
                              "3-rep1","4-rep1","5-rep1","6-rep1",
                              "3-rep2","4-rep2","5-rep2","6-rep2"),pch=pchs,col=cols,bty="n",y.intersp = 0.5,x.intersp = 0.5)
  legend("topleft",legend = c("1-rep1","2-rep1","3-rep1","4-rep1",
                              "1-rep2","2-rep2","3-rep2","4-rep2",
                              "3-rep1","4-rep1","5-rep1","6-rep1",
                              "3-rep2","4-rep2","5-rep2","6-rep2"),pch=pchs2,col="red3",bty="n",y.intersp = 0.5,x.intersp = 0.5)
  
}

par(mfrow=c(2,2),mar=c(4,4,4,4),lwd=2)
for(d in 1:length(Drugs))
{
  drug = Drugs[d]
  wells = as.character(Map[grep_exact(drug,Map[[3]]),1])
  concs = as.character(Map[match(wells,Map[[1]]),4])
  cols1 = c("dummy", base_cols[1:4])
  cols2 = c("dummy", base_cols[3:6])
  concs1 = c("conc 1","conc 2","conc 3","conc 4")
  concs2 = c("conc 3","conc 4","conc 5","conc 6")
  
  for(p in 1:length(Plates_OD_data))
  {
    plate_id_OD = names(Plates_OD_data)[p]
    p_OD = Plates_OD_data[[p]]
    p_Lux = Plates_Lux_data[[p]]
    
    p_OD_controls = p_OD[c(1,match(control_wells,names(p_OD)))]
    p_Lux_controls = p_Lux[c(1,match(control_wells,names(p_Lux)))]
    
    p_OD = p_OD[c(1,match(wells,names(p_OD)))]
    p_Lux = p_Lux[c(1,match(wells,names(p_Lux)))]
    
    if(p==1||p==2)
    {
      cols = cols1
      concs = paste0(as.character(Map[match(wells,Map[[1]]),5])," µg/ml")
    } else
    {
      cols = cols2
      concs = paste0(as.character(Map[match(wells,Map[[1]]),6])," µg/ml")
    }
    
    #plot the control wells
    med_OD_controls = sapply(as.data.frame(t(p_OD_controls[2:length(p_OD_controls)])),FUN=median)
    sd_OD_controls = sapply(as.data.frame(t(p_OD_controls[2:length(p_OD_controls)])),FUN=sd)
    plot(p_OD_controls[[1]],med_OD_controls,ylab="",xlab="Time",ylim=c(0,0.4),type="l",main=paste0(drug,", ",plate_id_OD),
         axes=F,frame=T,xlim=c(0,7),col="grey",lty=3)
    #arrows(x0=p_OD_controls[[1]],y0=(med_OD_controls+sd_OD_controls),x1=p_OD_controls[[1]],y1=(med_OD_controls-sd_OD_controls),
     #      length=0.05, angle=90, code=3,col="grey",lwd=1)
    axis(2,seq(0,0.4,by=0.1),las=2)
    mtext("OD",side=2,line=2.5)
    axis(1,seq(0,7,by=1),las=2)
    
    #plot the rest of the wells
    for(c in 2:length(p_OD))
    {points(p_OD[[1]],p_OD[[c]],type="l",lty=3,col=cols[c])}
    
    for(c in 2:length(p_OD))
    {
      par(new=T)
      if(c==2)
      {
        plot(p_Lux[[1]],p_Lux[[c]],ylab="",xlab="Time",ylim=c(0,500),type="l",main=paste0(drug,", ",plate_id_OD),
             axes=F,frame=T,col=cols[c],xlim=c(0,7),lty=1)
        axis(4,seq(0,500,by=100),las=2)
        mtext("Lux",side=4,line=2.5)
      } else
      {points(p_Lux[[1]],p_Lux[[c]],type="l",lty=1,col=cols[c])}
    }
    
    #plot the control wells
    med_Lux_controls = sapply(as.data.frame(t(p_Lux_controls[2:length(p_Lux_controls)])),FUN=median)
    sd_Lux_controls = sapply(as.data.frame(t(p_Lux_controls[2:length(p_Lux_controls)])),FUN=sd)
    
    points(p_Lux_controls[[1]],med_Lux_controls,type="l",col="grey",lty=1)
    #arrows(x0=p_Lux_controls[[1]],y0=(med_Lux_controls+sd_Lux_controls),x1=p_Lux_controls[[1]],y1=(med_Lux_controls-sd_Lux_controls),
     #      length=0.05, angle=90, code=3,col="grey",lwd=1)
    
    if(p==1||p==2)
    {legend("topleft",legend=c("no drug",concs),lty=1,col=c("grey",cols[2:5]),bty="n",y.intersp=0.7,x.intersp=0.5)} else
    {legend("topleft",legend=c("no drug",concs),lty=1,col=c("grey",cols[2:5]),bty="n",y.intersp=0.7,x.intersp=0.5)}
  }  #for(p in 1:length(Plates_OD_data)) 
}

dev.off()

#====== Plot all drugs =======

pdf(paste0(Out_dir,"All drugs.pdf"))
Drugs = unique(as.character(Map[[3]]))
par(mfrow=c(2,2),mar=c(4,4,4,4),lwd=2)
for(d in 1:length(Drugs))
{
  drug = Drugs[d]
  wells = as.character(Map[grep_exact(drug,Map[[3]]),1])
  concs = as.character(Map[match(wells,Map[[1]]),4])
  cols1 = c("dummy", base_cols[1:4])
  cols2 = c("dummy", base_cols[3:6])
  concs1 = c("conc 1","conc 2","conc 3","conc 4")
  concs2 = c("conc 3","conc 4","conc 5","conc 6")
  
  for(p in 1:length(Plates_OD_data))
  {
    plate_id_OD = names(Plates_OD_data)[p]
    p_OD = Plates_OD_data[[p]]
    p_Lux = Plates_Lux_data[[p]]
    
    p_OD_controls = p_OD[c(1,match(control_wells,names(p_OD)))]
    p_Lux_controls = p_Lux[c(1,match(control_wells,names(p_Lux)))]
    
    p_OD = p_OD[c(1,match(wells,names(p_OD)))]
    p_Lux = p_Lux[c(1,match(wells,names(p_Lux)))]
    
    if(p==1||p==2)
    {
      cols = cols1
      concs = paste0(as.character(Map[match(wells,Map[[1]]),5])," µg/ml")
    } else
    {
      cols = cols2
      concs = paste0(as.character(Map[match(wells,Map[[1]]),6])," µg/ml")
    }
    
    #plot the control wells
    med_OD_controls = sapply(as.data.frame(t(p_OD_controls[2:length(p_OD_controls)])),FUN=median)
    sd_OD_controls = sapply(as.data.frame(t(p_OD_controls[2:length(p_OD_controls)])),FUN=sd)
    plot(p_OD_controls[[1]],med_OD_controls,ylab="",xlab="Time",ylim=c(0,0.4),type="l",main=paste0(drug,", ",plate_id_OD),
         axes=F,frame=T,xlim=c(0,7),col="grey",lty=3)
    #arrows(x0=p_OD_controls[[1]],y0=(med_OD_controls+sd_OD_controls),x1=p_OD_controls[[1]],y1=(med_OD_controls-sd_OD_controls),
    #      length=0.05, angle=90, code=3,col="grey",lwd=1)
    axis(2,seq(0,0.4,by=0.1),las=2)
    mtext("OD",side=2,line=2.5)
    axis(1,seq(0,7,by=1),las=2)
    
    #plot the rest of the wells
    for(c in 2:length(p_OD))
    {points(p_OD[[1]],p_OD[[c]],type="l",lty=3,col=cols[c])}
    
    for(c in 2:length(p_OD))
    {
      par(new=T)
      if(c==2)
      {
        plot(p_Lux[[1]],p_Lux[[c]],ylab="",xlab="Time",ylim=c(0,500),type="l",main=paste0(drug,", ",plate_id_OD),
             axes=F,frame=T,col=cols[c],xlim=c(0,7),lty=1)
        axis(4,seq(0,500,by=100),las=2)
        mtext("Lux",side=4,line=2.5)
      } else
      {points(p_Lux[[1]],p_Lux[[c]],type="l",lty=1,col=cols[c])}
    }
    
    #plot the control wells
    med_Lux_controls = sapply(as.data.frame(t(p_Lux_controls[2:length(p_Lux_controls)])),FUN=median)
    sd_Lux_controls = sapply(as.data.frame(t(p_Lux_controls[2:length(p_Lux_controls)])),FUN=sd)
    
    points(p_Lux_controls[[1]],med_Lux_controls,type="l",col="grey",lty=1)
    #arrows(x0=p_Lux_controls[[1]],y0=(med_Lux_controls+sd_Lux_controls),x1=p_Lux_controls[[1]],y1=(med_Lux_controls-sd_Lux_controls),
    #      length=0.05, angle=90, code=3,col="grey",lwd=1)
    
    if(p==1||p==2)
    {legend("topleft",legend=c("no drug",concs),lty=1,col=c("grey",cols[2:5]),bty="n",y.intersp=0.7,x.intersp=0.5)} else
    {legend("topleft",legend=c("no drug",concs),lty=1,col=c("grey",cols[2:5]),bty="n",y.intersp=0.7,x.intersp=0.5)}
  }  #for(p in 1:length(Plates_OD_data)) 
}
dev.off()



