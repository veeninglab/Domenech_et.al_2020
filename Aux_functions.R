#List of functions defined here


#======================== Basics ====================

#This function takes an object (x) and a vector (v) and returns the first position of a word contain x in v, when x is found in v. 
grep_uniq <- function(x,v)
{
  match_vec = grep(x, v, fixed = TRUE)
  uniq_match = c()
  
  #To make sure that well_index points to the exact one and only one well
  if (length(match_vec) > 0)
  {
    for(i in 1:length(match_vec))
    {
      match_vec_temp = match_vec[i]
      if (nchar(match_vec_temp) > 0 && x == v[match_vec_temp])
      {uniq_match = c(uniq_match, match_vec_temp)}
    }
  }
  return(uniq_match[1])
}

#Like grep, but returns only exact match
grep_exact <- function(x,v)
{
  match_vec = grep(x, v, fixed = TRUE)
  uniq_match = c()
  
  #To make sure that well_index points to the exact one and only one well
  if (length(match_vec) > 0)
  {
    for(i in 1:length(match_vec))
    {
      match_vec_temp = match_vec[i]
      if (nchar(match_vec_temp) > 0 && x == v[match_vec_temp])
      {uniq_match = c(uniq_match, match_vec_temp)}
    }
  }
  return(uniq_match)
}

#Gets quadrant "quad" from the dataset. Dataset = 384 well plate data. quad = 1,2,3 or 4
quadr <- function(quad,dataset)  #NEW
{
  all_wells = names(dataset)[2:length(dataset)]
  sel_wells = c()

  nr_rows = 16
  nr_cols = length(all_wells)/nr_rows
  c = 1
  for(r in 1:nr_rows)
  {
    while(c <= (nr_cols*r))
    {
      if(quad == 1 && odd(r) && odd(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 2 && odd(r) && even(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 3 && even(r) && odd(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 4 && even(r) && even(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      c = c+1
    }
      
  }
  
  red_dataset = Red_dataset(dataset,sel_wells)
  return(red_dataset)
}

Red_dataset <- function(dataset, sel_wells)
{
  red_dataset = as.vector(dataset[1])
  for (well in 1:length(sel_wells))
  {
    well_name <- sel_wells[well]
    well_index <- grep_uniq(well_name,names(dataset))
    red_dataset <- cbind(red_dataset,dataset[well_index])
  }
  return (red_dataset)
}  

#Function returns a vector of 2 integers, ready to be the dimentions of a par
dimPar <- function(x,y)
{
  #x is the number of plots to be made
  #y is the number of the columns of the par
  if(x <= y)
  {
    nr_cols = x
    nr_rows = 1
  } else
  {
    if(x > y)
    {
      nr_cols = y
      nr_rows = ceiling(x/y)
    }
  }
  
  if (x == 96)
  {
    nr_cols = 12
    nr_rows = 8
  }
  
  c(nr_rows,nr_cols)
}

panel.empty <- function(x,y)
{
  xx_coord = max(x)/2
  yy_coord = max(y)/2
  text(xx_coord, yy_coord, "white", col = "white")
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="complete"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#Check: http://www.math.ist.utl.pt/~calves/courses/integra/capiii32.html
Trape_rule <- function(dataset)
{
  areas = c()
  
  xData = as.numeric(as.vector(dataset[[1]]))
  nr_points = length(xData)
  h = (xData[length(xData)]-xData[1])/(nr_points-1)
  
  for(i in 2:length(dataset))
  {
    yData = as.numeric(as.vector(dataset[[i]]))
    
    t = (yData[1] + yData[length(yData)])/2
    for(j in 2:(length(yData)-1))
    {
      t = t + yData[j]
    }
    
    Tn = t*h
    areas = c(areas,Tn)
  }
  names(areas) = names(dataset)[2:length(dataset)]
  return(areas)
  
}

NonLin_Rate <- function(dataset,min_growth=0.08,max_lag=9,growth_window=c(0,0.45),time_window=c(0,10),plotting=T,time_cutoff=0)
{
  xData = dataset[[1]]
  
  if(plotting==T)
  {par(mfrow=c(8,12),mar=c(1,1,1,1))}
  
  for(w in 2:length(dataset))
  {
    yData = dataset[[w]]
    if(plotting==T)
    {
      plot(xData,yData,ylim=growth_window,xlim=time_window,col="black",pch=19,cex=0.7,
           frame=T,axes=F,ylab="",xlab="")
      text(x=2,y=(growth_window[2]-0.05),labels=names(dataset)[w])
      points(c(time_cutoff,time_cutoff),growth_window,type="l",lty=2)
    } #if(plotting==T)
    
    #Check for growth
    if(max(yData[grep(F,xData > max_lag)]) <= min_growth)
    {parameters=c(0,0,NA,0,NA,NA,0,0)} else
    {
      nonlin_fit <- SummarizeGrowth(xData, yData)
      
      #nonlin_fit <- gcFitModel(xData, yData, gcID = "undefined",control=grofit.control(model.type=c("gompertz")))
      
      if(is.na(nonlin_fit$vals[[1]])==FALSE)
      {
        fitted_curve=predict(nonlin_fit$model)
        fitted_point=0
        fitted_area=0
        if(time_cutoff>0)
        {
          fitted_point = fitted_curve[length(grep(T,xData<time_cutoff))+1]
          fitted_area = nonlin_fit$vals$auc_e
        }
        
        if(plotting==T)
        {points(xData,fitted_curve,col="red3",type="l",lwd=2)}#ylim=growth_window,xlim=time_window)} #if(plotting==T)
        
        R = cor(yData,fitted_curve)
        R2 = summary(lm(yData~fitted_curve))$r.squared
        parameters=c(nonlin_fit$vals$r,nonlin_fit$vals$k,NA,nonlin_fit$vals$auc_l,
                     R2,R,fitted_point,fitted_area)
      } else
      {parameters=c(0,0,0,0,0,2,0,0)}
    }
    
    parameters = as.data.frame(parameters)
    row.names(parameters) = c("Growth_rate_h-1","plateauOD","lag-phase","AUCintegral","Rsq","PearsonCor","fitted_point","AUCtrap")
    names(parameters) = names(dataset)[w]
    if(w==2)
    {Parameters = as.data.frame(t(parameters))} else
    {Parameters = as.data.frame(rbind(Parameters,t(parameters)))}
  }
  
  return(Parameters)   
}

get_feature <- function(plate_ids,feature,Plate_database) #Feb 2015
{
  All_features = names(Plate_database)
  features = c()
  if(length(grep_exact(feature,All_features))>0)
  {
    column = grep_exact(feature,All_features)
    for(p in 1:length(plate_ids))
    {
      plate_id = plate_ids[p]
      features = c(features,as.character(Plate_database[grep_exact(plate_id,Plate_database[[1]]),column]))
    }
  }
  
  return(features)
}