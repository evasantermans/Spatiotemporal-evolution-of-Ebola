########################################################
# Modelling the Ebola outbreak in West Africa 
# SEIR model & MCMC estimation
########################################################
# Author: Eva Santermans
# Last update: July 17, 2015
########################################################
 
#######################
## DATA MANIPULATION ##
#######################


###### Extracting data and calculating input data.
###### NOTE! Data is also manipulated manually to correct typos. 


########################################################################################
## COUNTY CHOICE (Forecariah, Western Area Urb, Grand Cape Mount, Conakry)

GEOSEL = "Conakry"


########################################################################################
## LIBRARIES

library(xlsx)


########################################################################################
## DATA

## Reading original data
newdata=read.csv2("data_08-07-2015.csv")

## Subsetting district data 
subdata=subset(newdata,GEO_NAME==GEOSEL)

## Extracting data
vec1=subdata$Date[subdata$Status=="case"|subdata$Status=="Case"|subdata$Status=="CASE"]
dates1=as.Date(vec1,"%d/%m/%Y")
index1=is.na(dates1); tmp=as.Date(as.numeric(as.character(vec1[index1])),origin = "1899-12-30")
dates1[index1]=tmp

vec2=subdata$Date[subdata$Status=="Death"|subdata$Status=="DEATH"]
dates2=as.Date(vec2,"%d/%m/%Y")
index2=is.na(dates2); tmp=as.Date(as.numeric(as.character(vec2[index2])),origin = "1899-12-30")
dates2[index2]=tmp

reportingtimes=seq(as.Date("2014-05-23"),max(c(dates1,dates2)+56,na.rm=T),1)   
index1=reportingtimes%in%dates1
index2=reportingtimes%in%dates2

## Calculation new cases and deaths
cumcases=rep(NA,length(index1))
tmp=subdata$Count[subdata$Status=="case"|subdata$Status=="Case"|subdata$Status=="CASE"]
cumcases[index1]=tmp[!is.na(tmp)]
newcases=rep(NA,length(cumcases))
newcases[!is.na(cumcases)]=c(NA,diff(cumcases[!is.na(cumcases)]))
deaths=rep(NA,length(index2))
tmp=subdata$Count[subdata$Status=="Death"|subdata$Status=="DEATH"]
deaths[index2]=tmp[!is.na(tmp)]
newdeaths=rep(NA,length(deaths))
newdeaths[!is.na(deaths)]=c(NA,diff(deaths[!is.na(deaths)]))

## Dates
datesnum=c(reportingtimes-min(reportingtimes))
dates=reportingtimes

## Time and bookkeeping
time=seq(-1,max(datesnum,na.rm=T),precision)[-1]+(1-precision)
index1=rep(cumsum(!is.na(cumcases)),each=1/precision) # create an index vector to show what period newly reported cases refer to
index2=rep(cumsum(!is.na(deaths)),each=1/precision) # create an index vector to show what period newly reported deaths refer to
tmp1=table(index1)[2:max(index1)]*precision
tmp2=table(index2)[2:max(index2)]*precision

## Cutoff days effective reproduction number
days= pwc.int
max.date = max(c(dates1,dates2),na.rm=T); min.date = min(c(dates1,dates2),na.rm=T)
s.date = seq(from=min.date,to=max.date,by=days); l = length(s.date); date.cutoffs = s.date[2:l]
if (c(max.date-date.cutoffs[l-1])<=7){
  date.cutoffs=date.cutoffs[-c(l-1)]
}
n.breaks = length(date.cutoffs)+1
n.cutoffs = rep(NA,n.breaks)
for (i in 1:n.breaks){
  if(i==1){n.cutoffs[i] = sum((min(dates)+time)<date.cutoffs[i])}
  if(i>1&i<n.breaks){n.cutoffs[i] = sum((min(dates)+time)<date.cutoffs[i]&(min(dates)+time)>=date.cutoffs[i-1])}
  if(i==n.breaks){n.cutoffs[i] = sum((min(dates)+time)>=date.cutoffs[i-1])}
}


########################################################################################
## EXCLUDING OUTLIERS/TYPOS PER DISTRICT

nw = newcases[!is.na(newcases)]/tmp1; sort(nw,decreasing=F)[1:10]
cc = cumcases[!is.na(cumcases)]; cc 

nd = newdeaths[!is.na(newdeaths)]/tmp2; sort(nd,decreasing=F)[1:10]
cd = deaths[!is.na(deaths)]; cd


## Western Area Urb 
## ----------------
if (GEOSEL == "Western Area Urb"){
  cumcases[!is.na(cumcases)][c(109)]=2964
  newcases=rep(NA,length(cumcases))
  newcases[!is.na(cumcases)]=c(NA,diff(cumcases[!is.na(cumcases)]))
}

## Grand Cape Mount
## ----------------
if (GEOSEL == "Grand Cape Mount"){
  deaths[!is.na(deaths)][c(54)]=43; deaths[!is.na(deaths)][c(112:116)]=183
  newdeaths=rep(NA,length(deaths))
  newdeaths[!is.na(deaths)]=c(NA,diff(deaths[!is.na(deaths)]))
}

## Conakry
## -------

if (GEOSEL == "Conakry"){
  deaths[!is.na(deaths)][216]=302; 
  newdeaths=rep(NA,length(deaths))
  newdeaths[!is.na(deaths)]=c(NA,diff(deaths[!is.na(deaths)]))
}



########################################################################################
## MONOTONIZING USING THE  POOL ADJACENT VIOLATOR ALGORITHM
## A decrease in cumulative numbers is corrected by using the PAVA algorithm


pavit= function(pos=pos,tot=rep(1,length(pos))) { 
  pai1 = pai2 = pos/tot 
  N = length(pai1) 
  ni=tot 
  for(i in 1:(N - 1)) { 
    if(pai2[i] > pai2[i + 1]) { 
      pool = (ni[i]*pai1[i] + ni[i+1]*pai1[i + 1])/(ni[i]+ni[i+1]) 
      pai2[i:(i + 1)] = pool 
      k=i+1 
      for(j in (k - 1):1) { 
        if(pai2[j] > pai2[k]) { 
          pool.2 = sum(ni[j:k]*pai1[j:k])/(sum(ni[j:k])) 
          pai2[j:k] = pool.2 
        } 
      } 
    } 
  } 
  return(list(pai1=pai1,pai2=pai2)) 
}


cumcases[!is.na(cumcases)]=round(pavit(cumcases[!is.na(cumcases)])$pai2)
newcases=rep(NA,length(cumcases))
newcases[!is.na(cumcases)]=c(NA,diff(cumcases[!is.na(cumcases)]))

deaths[!is.na(deaths)]=round(pavit(deaths[!is.na(deaths)])$pai2)
newdeaths=rep(NA,length(deaths))
newdeaths[!is.na(deaths)]=c(NA,diff(deaths[!is.na(deaths)]))


