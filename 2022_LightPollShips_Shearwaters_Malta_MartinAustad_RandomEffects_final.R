#########################################################################################################################
##### RFID DATA OF YELKOUAN SHEARWATER IN MALTA ####################################
#########################################################################################################################
### original code provided by Martin Austad (MA)
### updated by Steffen Oppel on 8 Aug 2017
### updating for analysis of all years together started by MA on 05 Mar 2020
### updating started by MA on 21 Sep 2020 for purposes of further analysis and publication

library(data.table)
library(ggplot2)
library(ggforce)
library(reshape2)
library(dplyr)
library(tidyverse)
library(lubridate)
library(oce)
library(scales)
library(maptools)
library(sp)
library(fasttime)
library(MASS)
library(DHARMa)
library(lme4)
library(glmmTMB)

################################################################################################################
###################### YESH COLONY RFID DATA FROM MALTA   ######################################################
################################################################################################################

setwd("C:\\Users\\martin.austad\\Documents\\PHD\\RFID_Analysis_PHD")

## read in PIT deployments
deployments <-fread("PITdep2020.csv") #updated with 2020 deployments! Kept seperate cause added two additional columns- one bird has had three diff. pit tags (2nd replacment)

###Read in location data 
locs<-fread("Location_coordinates_SQM_RFID.csv")
locs<-SpatialPoints(locs[,c(6,5)], proj4string=CRS("+proj=longlat +datum=WGS84"))

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

###Read in RFID data 
#raw file in separate workspace file on github due to large size

raw <- fread("RFID_raw_data_all.txt", fill=TRUE, header=TRUE) #header: X added as header for first column; space in TAG ID removed and V1 to V3 added  
raw <- raw[,-c(10:12)]

#different uploads pasted into one file manually due to all the different lines created by PUTTY that don't go under the header as well as periods of overlap between files. All original putty files from which data was copied are stored in: C:\Users\Martin\Documents\working_folder\RFID\Analysis_ALL\uncombinedPuttyFiles
head(raw)

#breeding season 2017 was not in UTC; change to CEST was not done with a lag; UTC from deployment on 10 OCT 2017

raw17 <- raw %>%
  filter(Date>="2017-03-22"&Date<="2017-09-11") #filtering on character so to avoid specifying tz in POSIXct just yet

rawrest <- raw %>%
  filter(Date>="2017-10-10"&Date<="2020-08-01") #this approach is already removing garbage dates

raw17 <- raw17 %>%
  mutate(Time=trim.trailing(Time)) %>%
  mutate(DateTimelocal=as.POSIXct(paste(raw17$Date, raw17$Time, " "), format="%Y-%m-%d %H:%M:%S", tz="Europe/Berlin")) %>%
  mutate(DateTime=with_tz(DateTimelocal, tz = "UTC"))
#time was changed on RFID on 29-03-2017 at 15:50 UTC, so the conversion to UTC above is not correct for the period 26-03-2017 to 29-03-2017; but this period will not be used in analysis anyway

### INSERT THE MISSING TIME STAMPS ### still doesn't work for 1hr on 2017-03-26
raw17$DateTime[is.na(raw17$DateTime)]<-as.POSIXct(paste(raw17$Date[is.na(raw17$DateTime)], raw17$Time[is.na(raw17$DateTime)], " "), format="%Y-%m-%d %H:%M:%S", tz="UTC")

lag <- interval(ymd_hms("2017-03-26 01:00:00"), ymd_hms("2017-03-29 14:47:00"))
excess <- interval(ymd_hms("2017-03-26 02:00:00"), ymd_hms("2017-03-26 02:59:51"))

#lag needs plus 1 hour 
#excess needs - 1 hour

raw17$change <- "good"
raw17$change[4196:4306] <- "excesstime"
raw17$change[4307:8702] <- "lagtime" #8703 is an error and can be ignored

raw17$DateTimeS <- raw17$DateTime

#raw17 <- raw17 %>%
# if(change == "excesstime"){DateTimeS<-DateTime-dhours(1)}else(change == "lagtime"){DateTimeS<-DateTime+dhours(1)}		
for(i in 1:nrow(raw17)){
  if(raw17$change[i] == "excesstime"){
    raw17$DateTimeS[i]<-raw17$DateTime[i]-dhours(1)
  }
}  

for(i in 1:nrow(raw17)){
  if(raw17$change[i] == "lagtime"){
    raw17$DateTimeS[i]<-raw17$DateTime[i]+dhours(1)
  }
}  

raw17 <- raw17[,-c(10:12)] 
colnames(raw17)[10] <- 'DateTime'

rawrest$DateTime<-fastPOSIXct(paste(rawrest$Date, rawrest$Time, " "),tz="UTC",required.components=6L)		### only works if string is in yyyy mm dd order

raw <- rbind(raw17, rawrest)

rm(raw17)
rm(rawrest)

################################################################################################################
###################### MANIPULATE DATA AND DELETE UNNECESSARY ROWS   ###########################################
################################################################################################################

#old data with dates from 2017 breeding season seem to pop up in 2018; 2018 data in 2019; despite format made to reader
#check1 <- raw %>%
#filter(Detectiontype == "D") #to remove B - which are bad type codes 
#they are duplicates i.e first recorded in correct sequence but then appear later on out of seq 
#unique(raw$TagID) #a lot of garbage! 

deployments <- deployments %>% 
  filter(Age==4) ##keeps only adult birds and filters away any juveniles tagged

PITS <- as.factor(c(deployments$TagID_orig, deployments$TagID_2, deployments$TagID_3))

recordingperiod <-interval(ymd("2017-03-22", tz="UTC"),ymd("2020-08-01", tz = "UTC")) 


#checkX<- raw %>% filter(!DateTime %within% recordingperiod)
#checkX <- raw %>% filter(is.na(DateTime))

colnames(raw)[1] <- 'Detectiontype' #change name of first column

RFIDdata <- raw %>%
  filter(Detectiontype == "D") %>%
  filter(Type!="HW") %>% #HW are tags on birds; HA is the marker tag
  filter(TagID %in% unique(PITS)) %>% 
  filter(DateTime %within% recordingperiod) %>%
  filter(Ant!= "A0") %>%  #appears when testing RFID
  dplyr::select(DateTime, Detectiontype, Duration, Type, TagID, Ant, Count, Gap) %>%
  mutate(Gap=as.numeric(Gap))%>%
  arrange(DateTime)

head(RFIDdata)
tail(RFIDdata)

unique(RFIDdata$Ant)
length(unique(RFIDdata$TagID))

#########################################################################
####ENSURE THAT PIT TAG REPLACEMENTS ARE INCORPORATED INTO RECORDS#######
#########################################################################
#Modified from Steffen Oppel's CMR script

#PIT tags were attached to darvic rings using marine resin. If rings or resin were damaged, they were replaced.

#at least one PIT tag was re-used on different birds ("8000E1349EA964FE")
#PIT tag was taken out of broken resin and re-used 
#therefore replacement PIT needs to be used and not original; 
#but also to factor in replacement date; replace orig with new before replacement date

#as of 2020 one bird had its 3rd pit tag deployed - do this in two stages; replace orig with TAGID_2; then any TAGID_2 with TAGID_3

## CREATE REPLACEMENT LIST

deployments$ReplacementDate_1 <- as.Date(deployments$ReplacementDate_1, format = "%d/%m/%Y")
deployments$ReplacementDate_2 <- as.Date(deployments$ReplacementDate_2, format = "%d/%m/%Y")
deployments$DeployDate_orig <- as.Date(deployments$DeployDate_orig, format = "%d/%m/%Y")
RFIDdata$Date <- as.Date(RFIDdata$DateTime, tz="UTC")

replist<-deployments %>% 
  
  dplyr::filter(TagID_2!="") %>%
  
  rename(orig=TagID_orig, repl=TagID_2) %>%
  
  mutate(Period_orig=interval(DeployDate_orig,ReplacementDate_1)) %>%
  
  dplyr::select(orig,repl, Period_orig)


## UPDATE RECORDS

#yesh <- RFIDdata #backup before replacements

for(i in 1:nrow(RFIDdata)){
  if(RFIDdata$TagID[i] %in% replist$orig){
    if(RFIDdata$Date[i] %within% replist$Period_orig[match(RFIDdata$TagID[i],replist$orig)]){
      RFIDdata$TagID[i] <- replist$repl[match(RFIDdata$TagID[i],replist$orig)]}
  }
}

#####2nd time replacments (3rd pit on a bird; just one in 2020)

replist2<-deployments %>% 
  
  dplyr::filter(TagID_3!="") %>%
  
  rename(Secondtag=TagID_2, repl=TagID_3) %>% 
  
  mutate(Period_2tag=interval(DeployDate_orig, ReplacementDate_2)) %>%
  
  dplyr::select(Secondtag,repl, Period_2tag)

for(i in 1:nrow(RFIDdata)){
  if(RFIDdata$TagID[i] %in% replist2$Secondtag){
    if(RFIDdata$Date[i] %within% replist2$Period_2tag[match(RFIDdata$TagID[i],replist2$Secondtag)]){
      RFIDdata$TagID[i] <- replist2$repl[match(RFIDdata$TagID[i],replist2$Secondtag)]}
  }
}

###Filter out registrations done for testing and/or obtaining TAGID or just errors in RFID system######

deployments <- deployments %>%
  rename(ReplacementOne = `1stReplacement`, ReplacementTwo=`2ndReplacement`)

deployments$ReplacementOne[is.na(deployments$ReplacementOne)] <- 0
deployments$ReplacementTwo[is.na(deployments$ReplacementTwo)] <- 0

deployments <- deployments %>%
  mutate(CurrentTagID=dplyr::if_else(ReplacementOne==1, TagID_2, TagID_orig)) %>%
  mutate(DateUntil="01-08-2020") #change when including 2021 data

deployments <- deployments %>%
  mutate(CurrentTagID=dplyr::if_else(ReplacementTwo==1, TagID_3, CurrentTagID)) 

deployments$DateUntil <- as.Date(deployments$DateUntil, format = "%d-%m-%Y")

deployments <- deployments%>%
  mutate(depperiod=interval(DeployDate_orig,DateUntil)) #due to replaced tags having replaced the original, original deployment date needs to be used

RFIDdata$correct <- NA

RFIDdata[is.na(RFIDdata$TagID),]

RFIDdata[is.na(RFIDdata$Date),]

deployments[is.na(deployments$depperiod),]
deployments[is.na(deployments$CurrentTagID),]

currentPITS <- as.factor(deployments$CurrentTagID)

length(unique(RFIDdata$TagID))

RFIDdata <- RFIDdata %>%
  filter(TagID %in% currentPITS) 


###remove testing of tags during day 

RFIDdata[is.na(RFIDdata$DateTime),]

#introduce nightstarting first

midday<-as.POSIXct("12:00:00", format="%H:%M:%S", tz="UTC")   # create a reference time to split the day
midday<-format(midday, format="%H:%M:%S", tz="UTC")

RFIDdata<- RFIDdata %>%
  mutate(DateTime=with_tz(DateTime,tz="UTC")) %>%
  mutate(Time=format(DateTime, format="%H:%M:%S",tz = "UTC")) %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600))))

Sys.setenv(TZ = "UTC")

RFIDdata<- RFIDdata %>%
  #mutate(Time=format(DateTime, format="%H:%M:%S",tz = "UTC")) %>%
  #mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600)))) %>%
  mutate(sunset=sunriset(locs[1,], DateTime, direction=c("sunset"), POSIXct.out=T)[,2])%>%
  mutate(sunrise=sunriset(locs[1,], DateTime, direction=c("sunrise"), POSIXct.out=T)[,2])%>%
  mutate(duringday=if_else(Time>midday,if_else(DateTime<sunset,1,0),if_else(DateTime>sunrise,1,0)))

RFIDdata <- RFIDdata %>% 
  filter(duringday== 0) %>%
  dplyr::select(NightStarting,DateTime,Time,Detectiontype,Duration,Type,TagID,Ant,Count,Gap)

length(unique(RFIDdata$TagID))

################################################################################################################
###### RFID operating nights from marker tag  ##############################
##################################################################################################

RFIDoper <- raw %>%
  filter(Detectiontype == "D") %>%
  filter(Type=="HW") %>%
  filter(DateTime %within% recordingperiod) %>%
  dplyr::select(DateTime,Detectiontype,Duration,Type,TagID,Ant,Count,Gap) %>%
  filter(Ant!= "A0")   %>%                                       
  mutate(Gap=as.numeric(Gap))%>%
  arrange(DateTime)
head(RFIDoper)
tail(RFIDoper)

opernights<- RFIDoper %>%
  mutate(det=1) %>%
  mutate(Time=format(DateTime, format="%H:%M:%S", tz="UTC")) %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600))))%>%  ### 'ifelse' in base converts date to number!!
  group_by(NightStarting, Ant) %>%                              
  summarise(n_detections=sum(det)) %>%
  spread(key=Ant, value=n_detections) %>%
  arrange(NightStarting)


head(opernights)
tail(opernights)

length(opernights$NightStarting)
plot(opernights$NightStarting, opernights$A2) 
#OBS! Decrease in detection over season??? more birds entering so more tag collisions?
#RATHER - RFID is on for shorter along season - shorter nights - so less marker tag detections
#decrease in 2020 is due to change in setting - marker tag set at lower freq


################################################################################################################
###### SUMMARISE THE NUMBER OF INDIVIDUALS DETECTED AT EACH ANTENNA IN EACH NIGHT  ##############################
################################################################################################################

NIGHTDETECTIONS<- RFIDdata %>% #SUMMARISE THE MOVEMENTS FOR EACH NIGHT ACROSS ALL INDIVIDUALS
  mutate(det=1) %>%
  #mutate(Time=format(DateTime, format="%H:%M:%S")) %>%
  #mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600)))) %>%  ### 'ifelse' in base converts date to number!!
  group_by(NightStarting,Ant,TagID) %>%                              
  summarise(n_detections=sum(det)) %>%
  spread(key=Ant, value=n_detections) %>%
  arrange(NightStarting)


length(unique(NIGHTDETECTIONS$TagID))

NIGHTINDIVIDUALS<- NIGHTDETECTIONS %>%
  gather(key=Ant, value=n_detections, A1,A2,na.rm = T) %>%
  mutate(n_detections=1) %>%
  group_by(NightStarting,Ant) %>%
  summarise(N_individuals=sum(n_detections)) %>%
  spread(key=Ant, value=N_individuals) %>%
  arrange(NightStarting)

###################################################################################
##########merge marker tag data with number of PIT tag detections per Ant ##########
###################################################################################

#are nights were marker tag was not registered in one of the two Ant to be removed? 
#alternatively marker is only registered in one Ant on these nights cause of being moved out of correct position

opernights_detectionnights <- merge(opernights, NIGHTINDIVIDUALS, by = "NightStarting", all.x=TRUE, all.y = TRUE)
colnames(opernights_detectionnights)[2] <- 'A1_marker'
colnames(opernights_detectionnights)[3] <- 'A2_marker'
colnames(opernights_detectionnights)[4] <-  'A1_birds'
colnames(opernights_detectionnights)[5] <- 'A2_birds'

RFIDMETA <- raw %>%
  filter(Detectiontype=="E")%>%
  filter(Gap!="") %>%
  mutate(Time=format(DateTime, format="%H:%M:%S", tz="UTC")) %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600))))%>%
  group_by(NightStarting)%>%
  dplyr::select(NightStarting,DateTime,Detectiontype,Duration,Type,TagID,Ant)

#header names are for detections not for event rows so renaming needed
colnames(RFIDMETA)[4] <- 'Voltage'
colnames(RFIDMETA)[5] <- 'ClockAMP'
colnames(RFIDMETA)[6] <- 'A1AMP'
colnames(RFIDMETA)[7] <- 'A2AMP'

head(RFIDMETA)
tail(RFIDMETA)

unwanted <- c("period", "records", "change:", "DIP", "26")

RFIDAMPS <-  RFIDMETA %>%
  group_by(NightStarting) %>%
  filter(!ClockAMP %in% unwanted) %>%
  summarise(maxA1= max(A1AMP),maxA2=max(A2AMP), minvoltage=min(Voltage))


opernights_detectionnights <- merge(opernights_detectionnights, RFIDAMPS, by = "NightStarting", all.x=TRUE, all.y = TRUE)

#############################################################################
##############filter out nights when RFID was not working properly...########
#############################################################################

#...but night still has a row in opernights_detections

A1notworkingwell <- interval(ymd("2017-10-10", tz="UTC"), ymd("2017-12-09", tz="UTC")) #A1 not working well - just detection data in A2 not "in- out"

A2notworkingwell <- interval(ymd("2017-05-08", tz="UTC"), ymd("2017-05-10", tz="UTC")) #A2 with a few less bird detections

down <- interval(ymd("2018-12-01", tz="UTC"), ymd("2019-01-18", tz="UTC"))
                             
workingnights <- opernights_detectionnights %>%
  filter(!(NightStarting=="2017-03-22")) %>% #only a few hours on as testing
  filter(!(NightStarting %within% A2notworkingwell)) %>%
  filter(!(NightStarting=="2017-06-19")) %>% #failing voltage
  filter(!(NightStarting=="2017-08-24")) %>% #office
  filter(!(NightStarting=="2017-08-25")) %>%  #office
  #filter(!(NightStarting %within% A1notworkingwell) # %>%
  filter(!(NightStarting=="2018-01-23")) %>%
  filter(!(NightStarting=="2018-02-09")) %>%
  filter(!(NightStarting=="2018-02-13")) %>%
  filter(!(NightStarting %within% down)) %>%
  filter(!(NightStarting=="2020-07-14"))

################################################################################################################
########################### Monitoring nights ##############################
################################################################################################################


head(RFIDdata)

Sys.setenv(TZ = "UTC")

monitornights<- workingnights %>% #monitoring nights should not be based on RFIDdata or on markertag nights but workingnights 
  #group_by(NightStarting) %>%
  #filter(NightStarting %within% recordingperiod2) %>%
  #summarise(firstdet=min(DateTime), lastdet=max(DateTime))%>%
  mutate(NightEnding=NightStarting+days(1))%>%
  mutate(sunset=sunriset(locs[1,], as.POSIXct(NightStarting, tz="UTC"),direction=c("sunset"), POSIXct.out=T)[,2])%>%
  mutate(sunrise=sunriset(locs[1,], as.POSIXct(NightEnding, tz="UTC"),direction=c("sunrise"), POSIXct.out=T)[,2])%>%
  mutate(nightdur=interval(sunset,sunrise))%>%
  dplyr::select(NightStarting, NightEnding, sunset, sunrise, nightdur) %>%
  mutate(ships=0) 

################################################################################################################
###### DETERMINE ARRIVAL AND DEPARTURE FROM SEQUENTAL REGISTRATION AT A1 and A2   ##############################
################################################################################################################


### LOOP OVER EACH ANIMAL ###
YESH<-unique(RFIDdata$TagID)

### IN for A1 followed by A2 for same tag
RFIDdata[is.na(RFIDdata$DateTime),]

#errors created if an antenna is skipped (no registration).
#This error would be detected if on a respective night the antenna is not the same as the one on the previous night. 
#Problem: we don't necessarily know the date of the missed antenna and should therefore also mark the subsequent registration as an error

MOVES<-data.frame()

for (A in YESH){
  IND <- RFIDdata %>%
    filter(TagID == A) %>%
    arrange(DateTime) %>%
    
    #nextAntV<-RFIDdata$Ant[RFIDdata$TagID==A]
    
    mutate(nextAnt=c(Ant[-1],NA)) %>%
    mutate(nextTime=c(DateTime[-1],NA)) %>% ### insert column for date at next antenna to compare time  
    mutate(nextNight=c(NightStarting[-1],NA)) %>%
    mutate(Diffnight=ifelse(NightStarting==nextNight, 0, 1))%>%
    mutate(error=ifelse(Diffnight==1, ifelse(nextAnt==Ant, 0, 1),0)) %>%
    mutate(move=ifelse(Ant==nextAnt,0,1)) %>%
    mutate(direction=ifelse(move==1,ifelse(nextAnt=="A2","IN","OUT"),"STAY")) %>%
    mutate(timediff=difftime(nextTime,DateTime, units="secs")) %>% ### calculate time between two registrations
    mutate(MoveDateTime=if_else(timediff>3600,nextTime,DateTime))%>% ### insert condition that if time between A1 and A2 is > 60min then take the date for the next registration
    mutate(Time=format(DateTime, format="%H:%M:%S")) %>%
    mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600))))
  
  MOVES<-bind_rows(MOVES,IND)
  
}

head(MOVES)

################################################################################################################
###################### account for tagging during study ###########################################
################################################################################################################ 

n_deployed<- deployments %>%
  group_by(DeployDate_orig) %>%
  summarise(n=length(unique(CurrentTagID)))

# this function allows to calculate the proportion of birds that had been tagged by a given date
taggedby<-function(x){sum(n_deployed$n[n_deployed$DeployDate_orig<x])}

dim(RFIDdata)

################################################################################################
################filter everything to study period ##########################
################################################################################################
  
  per17 <- interval(ymd("2017-03-23", tz="UTC"), ymd("2017-06-30", tz="UTC"))
  per18 <- interval(ymd("2018-03-01", tz="UTC"), ymd("2018-06-30", tz="UTC"))
  per19 <- interval(ymd("2019-02-01", tz="UTC"), ymd("2019-06-30", tz="UTC")) 
  per20 <- interval(ymd("2020-02-01", tz="UTC"), ymd("2020-06-30", tz="UTC")) 

monitornightsper <- monitornights %>% 
    filter((NightStarting %within% per17) | (NightStarting %within% per18) | (NightStarting %within% per19) | (NightStarting %within% per20))

mopernights <- as.list(monitornightsper$nightdur)
monitornightsperlist <- as.list(monitornightsper$NightStarting)

length(unique(monitornightsper$NightStarting))

MOVES <- MOVES %>%
mutate(Year= year(NightStarting))

################################################################################################
#########load and manipulate ship data###########################################################
################################################################################################

#AIS data purchased from MarineTraffic for all the study period and for all vessels with AIS 
#AIS data was purchased for the area: 

#LAT 35.93797533610 and 35.97319135900, LON 14.30613217470 and 14.34027744020

AIS17 <- fread("AIS2017.csv")
AIS18 <- fread("AIS2018.csv")
AIS19 <- fread("AIS2019.csv")
AIS20 <- fread("AIS2020.csv")

AIS20 <- AIS20 %>%
  mutate(`TIMESTAMP UTC`=as.POSIXct(`TIMESTAMP UTC`, format="%d/%m/%Y %H:%M:%S", tz="UTC")) 

AIS <- rbind(AIS17, AIS18, AIS19)

AIS <- AIS %>%
  mutate(`TIMESTAMP UTC`=as.POSIXct(`TIMESTAMP UTC`, format="%Y-%m-%d %H:%M:%S", tz="UTC")) 

AIS <- rbind(AIS, AIS20)

colnames(AIS)[3] <- 'VESSELTYPE'
colnames(AIS)[10] <- 'DateTime'
colnames(AIS)[5] <- 'SPEEDKNOTS10'

### vessels for which the type was not clear where individually checked in the Marine Traffic online database

#shipcheck <- AIS %>%
  #filter(VESSELTYPE=="Fishery Patrol Vessel") #just replace with type needed to check
#unique(shipcheck$MMSI)

### only industrial ships kept (removing pleasure craft, sailing vessels and fishing vessels)
### passenger ships are cruiseliners, the 3 in the dataset which are over 90 meters long

Vesselsall <- c("Motor Hopper", "Asphalt/Bitumen Tanker", "Crude Oil Tanker", "Oil Products Tanker", "Bulk Carrier", "Container Ship", "General Cargo", "Tanker", "Bunkering Tanker", "Deck Cargo Ship", "LPG Tanker", "Oil/Chemical Tanker", "Vehicles Carrier", "Reefer", "Livestock Carrier", "Ro-Ro/Passenger Ship", "Chemical Tanker", "Hopper Dredger", "Inland, General Cargo maritime", "Cargo/Containership", "LPG/Chemical Tanker", "Ro-Ro Cargo", "Passenger Ship","Work Vessel", "Supply Vessel", "Offshore Supply Ship", "Tug/Supply Vessel")

AIS <- AIS %>%
  #mutate(DateTime=as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) %>%
  mutate(Time=format(DateTime, format="%H:%M:%S",tz = "UTC")) %>%
  mutate(Date=as.Date(DateTime)) %>%
  mutate(year=year(DateTime)) %>%
  filter(VESSELTYPE %in% Vesselsall)#%>%


ships<-data.frame()

ship <- unique(AIS$SHIPNAME)

for (A in ship){
  S <- AIS %>%
    filter(SHIPNAME == A) %>%
    arrange(DateTime) %>%
    mutate(nextTime=c(DateTime[-1],NA))%>% ###
    mutate(timediff=difftime(nextTime,DateTime, units="secs")) %>% ### calculate time between two registrations insert column for date at next record to compare time  
    mutate(nextLON=c(LON[-1],NA))%>%
    mutate(nextLAT=c(LAT[-1],NA)) %>%
    mutate(Time=format(DateTime, format="%H:%M:%S")) %>%
    mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600))))
  
  ships<-bind_rows(ships,S)
  
}   
dim(ships)

ships <- ships %>%
  mutate(nextLON=ifelse(is.na(nextLON),LON,nextLON))%>%
  mutate(nextLAT=ifelse(is.na(nextLAT),LAT,nextLAT))

ships[is.na(ships$nextLAT),]

p1 <- SpatialPoints(ships[,c(6,7)], proj4string=CRS("+proj=longlat +datum=WGS84"))
p2 <- SpatialPoints(ships[,c(16,17)], proj4string=CRS("+proj=longlat +datum=WGS84"))      

for(i in 1:nrow(ships)){
  ships$distance[i] <- spDists(p1[i],p2[i])
}

ships[is.na(ships$distance),]

td <- c(300:2000)
#ships[is.na(ships$timediff),]

#need to create loop based not only on ship but on visit ID (for vessels that visit more than once)
#within each visit ID take first and last time and pos when stationary 

nights <- as.list(monitornights$nightdur)

ships <- ships %>%
  group_by(SHIPNAME) %>%
  mutate(stationary=ifelse(timediff %in% td & distance<0.08, "STAY", "MOVED"))%>% #timediff %in% td  removed because even with stationary ships can have short time resolution proabably, due to change in heading? but then stricter movements #previously distance<0.1 & SPEEDKNOTS10<20
  mutate(stationary2=ifelse(distance<0.1 & SPEEDKNOTS10<15, "STAY", "MOVED"))%>%  #timediff %in% td  removed because even with stationary ships can have short time resolution proabably, due to change in heading? #previously distance<0.2 & SPEEDKNOTS10<35
  group_by(NightStarting, SHIPNAME) 

ships1 <- ships %>%
  filter(stationary=="STAY")%>%
  filter(DateTime %within% nights) %>% 
  group_by(SHIPNAME, NightStarting)  %>%  
  summarise(firstdet=min(DateTime), lastdet=max(DateTime))

ships2 <- ships %>%
  filter(stationary2=="STAY")%>%
  filter(DateTime %within% nights) %>% 
  group_by(SHIPNAME, NightStarting)  %>%  
  summarise(firstdet=min(DateTime), lastdet=max(DateTime))

ships1 <- as.data.frame(ships1)
ships2 <- as.data.frame(ships2)

ships1 <- ships1 %>% mutate(presence=interval(firstdet,lastdet))
ships2 <- ships2 %>% mutate(presence=interval(firstdet,lastdet))

#### COUNT THE NUMBER OF SHIPS IN EACH MONITORED NIGHT ###
monitornightsper$ships = 0
monitornightsper$ships2 <- 0

for (n in 1:length(monitornightsper$sunrise)){
  x<-monitornightsper$nightdur[n]
  shipsthisnight<-int_overlaps(ships1$presence,x)
  shipsthisnight2<-int_overlaps(ships2$presence,x)
  monitornightsper$ships[n]<- sum(shipsthisnight, na.rm=T)
  monitornightsper$ships2[n]<-sum(shipsthisnight2, na.rm=T)
}

length(unique(ships1$SHIPNAME))
length(unique(ships2$SHIPNAME))

#we will use ships2 

#ships <- ships2 
rm(ships1)
head(monitornightsper)
monitornightsper <- monitornightsper[,-c(6)]
names(monitornightsper)[6]<-"ships"

#### Present ship nights in graph/table for report
monitornightsper_summary <- monitornightsper %>%
  filter(ships>0) %>%
  dplyr::select(NightStarting,ships) %>%
  mutate(Year=year(NightStarting))


length(unique(monitornightsper_summary$NightStarting))

monitornightsper_summary <- monitornightsper_summary %>% 
  mutate(shipsN = ifelse(ships>1, 'Bunkering events', 'One ship for entire night duration'))

monitornightsper_summary%>%
  ggplot(aes(x=NightStarting, y=ships, width=1, fill= shipsN))+
  geom_col()+
  facet_wrap('Year', scales="free")+
  ylab("Number of ships stationary in front of shearwater colony") +
  scale_x_date(name="Night Starting Date", labels=date_format("%m-%d"), breaks=date_breaks(width = "1 week"))+ 
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y = element_text(size=12, color="black"), 
        axis.text.x = element_text(size=12, color="black", angle = 90, vjust = 0.5),
        axis.title=element_text(size=14),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)), 
        strip.text.x=element_text(size=12, color="black"),
        strip.text.y=element_text(size=12, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        legend.title = element_blank())

################################################################################################################
###### LOAD AND MANIPULATE LIGHT METER DATA   ##############################
################################################################################################################
#a sky quality meter (Unihedron) was used to measure light levels at the cliff face

sqm17<-fread("MT24_SQM_17.csv") #rm(sqm17)
sqm18<-fread("MT24_SQM_18.csv") #rm(sqm18)
sqm19<-fread("MT24_SQM_19.csv") #rm(sqm19)
sqm20<-fread("MT24_SQM_20.csv") #rm(sqm20)

sqm <- rbind(sqm17,sqm18,sqm19,sqm20)
names(sqm)[7]<-"MAG"

#####Adding candela per m2 to SQM data ####################

#10.8*10^4*10^(-0.4*22)

sqm$candelaperm2 <- 10.8*10^4*10^(-0.4*sqm$MAG)
10.8*10^4*10^(-0.4*16)

## need to remove all daytime measurements and Account for moon phase

sqmnotworkingwell <- interval(ymd("2020-01-15", tz="UTC"), ymd("2020-03-14", tz="UTC"))


SQM<-sqm %>%
  mutate(DateTime=dmy_hms(paste(date_UTC, time_UTC,sep=" "),tz = "UTC")) %>% 
  mutate(Time=format(DateTime, format="%H:%M:%S",tz = "UTC")) %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600)))) %>%
  filter(volts>3.5) %>%     ## remove locs when battery was dead
  mutate(sunset=crepuscule(locs[2,], DateTime,direction=c("dusk"), solarDep=18, POSIXct.out=T)[,2])%>%
  mutate(sunrise=crepuscule(locs[2,], DateTime,direction=c("dawn"), solarDep=18, POSIXct.out=T)[,2])%>%
  mutate(duringday=if_else(Time>midday,if_else(DateTime<sunset,1,0),if_else(DateTime>sunrise,1,0)))%>%
  filter(duringday==0) %>%
  filter(!(NightStarting=="2018-05-15")) %>%
  filter(!(NightStarting=="2018-05-16"))%>% #stones were covering the SQM on two nights 15 & 16 May, much darker than all other nights
  filter(!(NightStarting %within% sqmnotworkingwell))%>%
  mutate(prop.illuminated=moonAngle(t=DateTime, longitude=coordinates(locs)[2,1], latitude=coordinates(locs)[2,2])$illuminatedFraction)%>%### This is how 'full' the moon is - 1 is full moon, 0 is new moon
  mutate(moon.elevation=moonAngle(t=DateTime, longitude=coordinates(locs)[2,1], latitude=coordinates(locs)[2,2])$altitude) %>%
  mutate(moon.light=prop.illuminated*(moon.elevation/coordinates(locs)[2,2]))


max(SQM$moon.light)
min(SQM$moon.light)

SQM <- SQM %>% 
  mutate(Year=year(DateTime)) %>%
  filter((NightStarting %within% per17) | (NightStarting %within% per18) | (NightStarting %within% per19) | (NightStarting %within% per20))

length(unique(SQM$NightStarting))

################################################################################################################
###### QUESTION 1: DOES SHIP BUNKERING AFFECT LIGHT METER DATA?   ##############################
###############################################################################################################

## summarise the number of ships within the 10-minute intervals of light readings

SQM$ships<-0
x<-SQM$DateTime[1]-minutes(10)
x1<-SQM$DateTime[1]
sqmint<-interval(x,x1)
shipsthisinterval<-int_overlaps(ships$presence,sqmint)
SQM$ships[1]<- sum(shipsthisinterval, na.rm=T)

for (n in 2:length(SQM$DateTime)){
  x<-SQM$DateTime[n]-minutes(10) #originally   x<-SQM$DateTime[n-1] but this creates a large number for the first record of the day
  x1<-SQM$DateTime[n]
  sqmint<-interval(x,x1)
  shipsthisinterval<-int_overlaps(ships$presence,sqmint)
  SQM$ships[n]<- sum(shipsthisinterval, na.rm=T)
}

max(SQM$ships)

windowsFonts()
windowsFonts(sans=windowsFont("Helvetica"))
fonts <- list(sans = "Helvetica",
              mono = "TT Courier New",
              serif = "TT Times New Roman")

m <- SQM %>% 
  ggplot(aes(x=moon.light, y=candelaperm2))+
  geom_point(aes(size=ships, color=ships)) +
  scale_size(range=c(0,4), limits=c(0,14), breaks=c(0,5,10,14)) +
  scale_color_gradient(low="black", high="orange", limits=c(0,14), breaks=c(0,5,10,14))+
  guides(color=guide_legend(title="N. ships"), size = guide_legend(title="N. ships"))+
  geom_vline(aes(xintercept=0), color='black', size=1) +
  annotate("text", x=1, y=0.02, label= "moon in the sky", size=7, family="sans") +  #when using MAG use y=23
  annotate("text", x = -0.5, y=0.02, label = "before moonrise", size=7, family="sans")+ #when using MAG use y=23
  xlab("Moonlight (prop. of moon illuminated*degrees elevation)") + #moon light (prop. of moon illuminated*degrees elevation)
  ylab(bquote('Cliff face brightness '~(cd/m^2)))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        text=element_text(family="sans"),
        axis.text = element_text(size=18, color="black", family="sans"), 
        axis.title=element_text(size=22, color="black", family="sans"), 
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)),
        strip.text.x=element_text(size=18, color="black", family="sans"), 
        strip.background=element_rect(fill="white", colour="black"), 
        legend.background = element_rect(fill="white", colour="white"),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=16, family="sans"),
        legend.title= element_text(size=18, family="sans"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


m


#### STATISTICAL TEST OF BUNKERING EFFECT ON CLIFF BRIGHTNESS###

head(SQM)
SQM$bunker<-ifelse(SQM$ships>0,1,0)
hist(SQM$MAG)

hist(SQM$candelaperm2)
hist(log(SQM$candelaperm2))

SQM$HR <- as.factor(hour(SQM$DateTime))
SQM$Year <- as.factor(SQM$Year)
SQM$NightStarting <- as.factor(SQM$NightStarting)
SQM$MONTH <- as.factor(month(SQM$DateTime))

Q1m_candela<-glmmTMB(candelaperm2~moon.light*bunker+(1|Year)+(1|MONTH)+(1|NightStarting/HR), data=SQM, family=Gamma(link="log")) 
outQ1m_candela<-summary(Q1m_candela)
outQ1m_candela
outQ1m_candela$coefficients

full <- glmmTMB(candelaperm2~ships+moon.light*bunker+(1|Year)+(1|MONTH)+(1|NightStarting/HR), data=SQM, family=Gamma(link = "log")) 
outfull<-summary(full)
outfull
outfull$coefficients

anova(Q1m_candela, full, test="F")

simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)
testOutliers(simulationOutput, type="bootstrap")
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

################################################################################################################
###### QUESTION 2: ARE NUMBER OF MOVEMENTS RELATED TO SQM LIGHT LEVEL MEASUREMENTS? ###########
################################################################################################################

###previously this question used all "in"s by birds but this might contain some noise i.e. birds that move around in entrance

#we want to get the first 'in' by a bird during a night which is the arrival to the colony (which is not following an 'out' - bird could be inside during the day and go to entrance)

MIOr <-MOVES %>%
  #filter(error!=1) %>% ##DON't remove errors now - just moves first in or last out to another non-error move but might not be the last move
  filter(move==1)%>%
  group_by(TagID, NightStarting) %>%
  mutate(firstmovetime=min(DateTime))%>% #changing to Datetime
  mutate(firstmove=ifelse(DateTime>firstmovetime, 0, 1)) %>% #changing to Datetime to avoid duplicates
  mutate(lastmovetime=max(DateTime))%>% #changing to Datetime
  mutate(lastmove=ifelse(DateTime==lastmovetime, 1, 0))#changing to Datetime
#filter(firstmove==1 | lastmove ==1)

for(i in nrow(MIOr):2){   #####the move on the next night should not be used
  if(MIOr$error[i] == 1){
    MIOr$error[i+1] <-1
  }
}

MIOperall <- MIOr %>%
  filter(DateTime %within% mopernights)

###to obtain number of entrances with/without error detection
#MovesIN <- MIOperall %>% 
  #filter(direction=="IN") %>%
  #filter(firstmove==1) %>%
  #mutate(count=1)

#MovesINerr <- MIOperall %>% 
  #filter(direction=="IN") %>%
  #filter(firstmove==1) %>%
  #filter(error==1) %>%
  #mutate(count=1) 

#MovesINg <- MIOperall %>% 
  #filter(direction=="IN") %>%
  #filter(firstmove==1) %>%
  #filter(error!=1) %>%
 # mutate(count=1) 

#length(unique(MovesIN$TagID))
#length(unique(MovesINg$TagID))

#Q2 with a higher solar angle than Q1 so all birds entering cave are included
#therefore solardep changed to 6 - end of civil twilight

SQMSS<-sqm %>%
  mutate(DateTime=dmy_hms(paste(date_UTC, time_UTC,sep=" "),tz = "UTC")) %>% 
  mutate(Time=format(DateTime, format="%H:%M:%S",tz = "UTC")) %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600)))) %>%
  mutate(Year=year(DateTime)) %>%
  filter(volts>3.5) %>%     ## remove locs when battery was dead
  mutate(sunset=crepuscule(locs[2,], DateTime,direction=c("dusk"), solarDep=6, POSIXct.out=T)[,2])%>%
  mutate(sunrise=crepuscule(locs[2,], DateTime,direction=c("dawn"), solarDep=6, POSIXct.out=T)[,2])%>%
  mutate(duringday=if_else(Time>midday,if_else(DateTime<sunset,1,0),if_else(DateTime>sunrise,1,0)))%>%
  filter(duringday==0) %>%
  #filter(!(NightStarting=="2018-05-15")) %>% #DO FILTERING LATER
  #filter(!(NightStarting=="2018-05-16"))%>% #stones were covering the SQM on two nights 15 & 16 May, much darker than all other nights
  #filter(!(NightStarting %within% sqmnotworkingwell))%>%
  mutate(prop.illuminated=moonAngle(t=DateTime, longitude=coordinates(locs)[2,1], latitude=coordinates(locs)[2,2])$illuminatedFraction)%>%### This is how 'full the moon is - 1 is full moon, 0 is new moon
  mutate(moon.elevation=moonAngle(t=DateTime, longitude=coordinates(locs)[2,1], latitude=coordinates(locs)[2,2])$altitude) %>%
  mutate(moon.light=prop.illuminated*(moon.elevation/coordinates(locs)[2,2]))

#filter SMSS to when RFID is working

SQMSS <-SQMSS %>%
  filter(NightStarting %in% monitornightsperlist)

length(unique(SQMSS$NightStarting)) 

SQMSS$ships <- 0 
for (n in 2:length(SQMSS$DateTime)){
  x<-SQMSS$DateTime[n]-minutes(10)
  x1<-SQMSS$DateTime[n]
  sqmint<-interval(x,x1)
  shipsthisinterval<-int_overlaps(ships$presence,sqmint)
  SQMSS$ships[n]<- sum(shipsthisinterval, na.rm=T)
}

MOVESIN <- MIOperall %>% 
  filter(firstmove==1) %>%
  filter(direction=="IN") %>% 
  filter(error!=1)%>%
  mutate(count=1) 

SQMSS$bunker<-ifelse(SQMSS$ships>0,1,0)

SQMSSQ3 <- SQMSS #to use later for Q3

SQMSS <- SQMSS %>% #to filter away nights SQM was not working well 
  filter(!(NightStarting=="2018-05-15")) %>% 
  filter(!(NightStarting=="2018-05-16"))%>% #stones were covering the SQM on two nights 15 & 16 May, much darker than all other nights
  filter(!(NightStarting %within% sqmnotworkingwell))

length(unique(SQMSS$NightStarting)) 

SQMSShr<- SQMSS %>%
  mutate(HR=as.factor(hour(DateTime)))%>%
  mutate(MONTH=as.factor(month(NightStarting)))%>% #, label=T, abbr=T, used as.factor instead
  group_by(MONTH,NightStarting,HR) %>%
  summarise(MAG=mean(MAG),candelaperm2=mean(candelaperm2),prop.illuminated=mean(prop.illuminated), moon.elevation=mean(moon.elevation), ships=max(ships), bunker=max(bunker))%>%
  mutate(moon.light=prop.illuminated*(moon.elevation/coordinates(locs)[2,2]))

MovesIN_hr <- MOVESIN%>% 
  #filter(firstmove==1) %>%
  #filter(direction=="IN") %>%
  #mutate(count=1) %>%
  mutate(HR=hour(DateTime))%>%
  mutate(MONTH=month(NightStarting))%>% 
  group_by(MONTH,NightStarting,HR) %>%
  summarise(activity=length(unique(TagID)))

MOVE_SQMSS<-merge(SQMSShr,MovesIN_hr,by=c('MONTH','NightStarting','HR'), all.x=T)
MOVE_SQMSS[is.na(MOVE_SQMSS$activity),]
MOVE_SQMSS$activity[is.na(MOVE_SQMSS$activity)]<-0

levels(MOVE_SQMSS$HR)

SQMSS <- SQMSS %>%
  mutate(HR=hour(DateTime))%>%
  mutate(MONTH=month(NightStarting))

MOVE_SQMSS$HRordered <- factor(MOVE_SQMSS$HR, levels = c("17", "18", "19","20","21","22","23","0","1","2","3","4", "5")) 
SQMSS$HRordered <- factor(SQMSS$HR, levels = c("17", "18", "19","20","21","22","23","0","1","2","3","4", "5"))
SQMSS$MONTH <- factor(SQMSS$MONTH, levels=c("2", "3", "4", "5", "6"))

hr <- c(17,18,19,20,21,22,23,0,1,2,3,4)
MOVE_SQMSS$HRordered <- factor(MOVE_SQMSS$HR, levels = c("17", "18", "19","20","21","22","23","0","1","2","3","4"))

MOVE_SQMSS <- MOVE_SQMSS %>%
  filter(HR %in% hr)

sum(MOVE_SQMSS$activity) #no activity at 5 hrs

colnames(MOVE_SQMSS)[11] <- 'IN_activity'

hist(MOVE_SQMSS$IN_activity) 


MOVE_SQMSS$Year <- as.factor(year(MOVE_SQMSS$NightStarting))
MOVE_SQMSS$NightStarting <- as.factor(MOVE_SQMSS$NightStarting)
MOVE_SQMSS$MONTH <- as.factor(month(MOVE_SQMSS$NightStarting))

#negative binomial glmm
Q2m_candela_nb <- glmmTMB(IN_activity~-1+HRordered+MONTH+candelaperm2+(1|NightStarting/HR)+(1|MONTH)+(1|Year), data=MOVE_SQMSS, family=nbinom1)
outQ2m_candela_nb<-summary(Q2m_candela_nb)
outQ2m_candela_nb$coefficients
outQ2m_candela_nb

simulationOutput <- simulateResiduals(fittedModel = Q2m_candela_nb, plot = F)
plot(simulationOutput)
testOutliers(simulationOutput, type="bootstrap")
testDispersion(simulationOutput)
testZeroInflation(simulationOutput) 

#############################################
#############Q2 entrance activity############
#############################################

MovesE <- MIOperall %>% 
  filter(lastmove==0 & firstmove==0) %>%
  mutate(count=1) 

#there are no identified errors cause the errors identified were all in the first or last move...

length(unique(MovesE$TagID))

MovesENTRANCE_hr <- MIOperall %>%
  filter(lastmove==0 & firstmove==0) %>%
  mutate(count=1) %>%
  mutate(HR=hour(DateTime))%>%
  mutate(MONTH=month(NightStarting))%>% 
  group_by(MONTH,NightStarting,HR) %>%
  summarise(entrance_activity=length(unique(TagID))) # a large expected difference 

MOVE_SQMSSe<-merge(SQMSShr,MovesENTRANCE_hr,by=c('MONTH','NightStarting','HR'), all.x=T)
MOVE_SQMSSe[is.na(MOVE_SQMSSe$entrance_activity),]
MOVE_SQMSSe$entrance_activity[is.na(MOVE_SQMSSe$entrance_activity)]<-0
MOVE_SQMSSe$HRordered <- factor(MOVE_SQMSSe$HR, levels = c("17", "18", "19","20","21","22","23","0","1","2","3","4", "5")) 


#plot number of entrance movements per hour
MOVE_SQMSSe %>% 
  ggplot(aes(x=HRordered, y=entrance_activity))+
  #geom_point(colour="black") +
  ylab("N individuals moving in cave entrance per hour") +
  xlab("Hour of the night (UTC)") +
  geom_boxplot(colour="black") +
  #facet_wrap('HRordered',scales = "fixed", shrink = TRUE, ncol=2)+
  #geom_smooth(fill="lightblue", size=1.5, method='lm', se=T)+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=15, color="black"), 
        axis.title=element_text(size=15),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)), 
        strip.text.x=element_text(size=18, color="black"),
        strip.text.y=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

hist(MOVE_SQMSSe$entrance_activity)

MOVE_SQMSSe$Year <- as.factor(year(MOVE_SQMSSe$NightStarting))
MOVE_SQMSSe$NightStarting <- as.factor(MOVE_SQMSSe$NightStarting)

#nb mixed effect model
Q2m_candela_e <- glmmTMB(entrance_activity~-1+HRordered+MONTH+candelaperm2+(1|NightStarting/HR)+(1|Year), ziformula=~1, family=nbinom2, data=MOVE_SQMSSe) #(1|MONTH) very low variance leading to convergence issue
outQ2m_candela_e<-summary(Q2m_candela_e)
outQ2m_candela_e$coefficients
summary(outQ2m_candela_e)
outQ2m_candela_e

simulationOutput <- simulateResiduals(fittedModel = Q2m_candela_nb, plot = F)
plot(simulationOutput)
testOutliers(simulationOutput, type="bootstrap")
testDispersion(simulationOutput)
testZeroInflation(simulationOutput) 
################################################################################################################
###### MANUALLY SELECTING TIME INTERVALS AROUND BUNKERING EVENTS   ##############################
################################################################################################################

#####Find length of bunkering events#####
#ie consecutive nights when bunkering took place

events <- monitornightsper %>%
  filter(ships>1)%>%
  arrange(NightStarting) %>%
  mutate(diff = c(0,diff(NightStarting)),
         eventID = 1 + cumsum(diff > 1)) %>%
  group_by(eventID) %>%
  summarise(eventlength = 1+(last(NightStarting) - first(NightStarting)),
            firstNightStarting=first(NightStarting),
            lastNightStarting=last(NightStarting))

events_1ship <- monitornightsper %>% #nights with only 1 ship
  filter(ships==1)

Oneshipperiods <- events_1ship$nightdur

ONESHIP <-data.frame()
for(i in 1:length(Oneshipperiods)){
  x<- AIS %>%
    filter(DateTime %within% Oneshipperiods[i])
  ONESHIP<-rbind(ONESHIP,x)
}

ONESHIP <- ONESHIP %>%
  mutate(NightStarting=if_else(Time>midday,as.Date(DateTime),as.Date(DateTime-(24*3600)))) %>%
  group_by(NightStarting) %>%
  summarise(numberships=length(unique(SHIPNAME))) 
  

length(unique(ONESHIP$NightStarting)) 

#in practice this is a check on whether the filtering done to obtain ships df makes sense
#these nights indeed only have 1 ship for the whole duration of the night

withships<-as.numeric(rownames(monitornightsper)[monitornightsper$ships>0]) #still important to put zero here otherwise 1 ship nights migth be selected as the before or after event night
before<-withships-1
after<-withships+1
selection<-unique(c(withships,before,after))
selection<-selection[selection>0]

CaseControl<-monitornightsper%>%
  filter(row_number() %in% selection)#%>%
  #filter(ships!=1)#%>% #remove nights with only one ship
  #select(NightStarting,sunset,sunrise,ships,nightdur)#%>%

################################################################################################################
###### QUESTION 3: Are number of movements related to ship presence?  #####################
################################################################################################################

### remove the lines outside the Case Control time periods

SQMcc<-data.frame()
for(i in 1:length(selection)){
  x<-SQMSSQ3 %>% #includes twilights &nights when SQM was down but bunkering took place 2018-05-15
    filter(DateTime %within% CaseControl$nightdur[i])
  SQMcc<-rbind(SQMcc,x)
}
length(unique(SQMcc$NightStarting)) 

SQMcc <- SQMcc %>%
 mutate(HR=as.factor(hour(DateTime)))%>%
  mutate(MONTH=as.factor(month(NightStarting)))%>% #, label=T, abbr=T, used as.factor instead
  group_by(MONTH,NightStarting,HR) %>%
  summarise(MAG=mean(MAG),candelaperm2=mean(candelaperm2),prop.illuminated=mean(prop.illuminated), moon.elevation=mean(moon.elevation), ships=max(ships), bunker=max(bunker))%>%
  mutate(moon.light=prop.illuminated*(moon.elevation/coordinates(locs)[2,2]))


colnames(MovesIN_hr)[4] <- 'IN_activity'
Mcc<-merge(SQMcc,MovesIN_hr, by=c('MONTH','NightStarting','HR'), all.x=T)

Mcc$IN_activity[is.na(Mcc$IN_activity)]<-0

Mcc <- Mcc %>% 
  filter(HR %in% c(17,18,19, 20, 21, 22, 23, 0, 1, 2, 3,4)) 


#generate bunkering event ID for use as random effect in model
Mcc <- Mcc %>%
  arrange(NightStarting) %>%
  mutate(diff = c(0,diff(NightStarting)),
         eventID = 1 + cumsum(diff > 1))

oneshipnights <- as.list(ONESHIP$NightStarting)

Mcc <- Mcc %>%
  filter(!NightStarting %in% oneshipnights)

#remove manually nights before and after 1 ship nights that were not part of longer events
Mcc<- Mcc %>%
  filter(NightStarting!="2018-03-01") %>% #2018-03-01 first night functioning after RFID failed on bunkering night 18062020
  filter(NightStarting!="2018-04-04") %>%
  filter(NightStarting!="2018-04-06") %>%
  filter(NightStarting!="2018-05-02") %>%
  filter(NightStarting!="2020-02-05") %>%
  filter(NightStarting!="2020-02-07") %>%
  filter(NightStarting!="2020-02-20") %>%
  filter(NightStarting!="2020-02-22") 

#split the last event into 2 
Mcc$eventID[1215:1254] <- 27 

Mcc <- Mcc %>%
  mutate(bunker=ifelse(ships>0,'ships present','no ships'))

Mcc$HRordered <- factor(Mcc$HR, levels = c("17","18","19", "20","21","22","23","0","1","2","3", "4"))

Mcc$monthtxt <- NA

Mcc$monthtxt[Mcc$MONTH==2] <- "Feb"
Mcc$monthtxt[Mcc$MONTH==3] <- "Mar"
Mcc$monthtxt[Mcc$MONTH==4] <- "Apr"
Mcc$monthtxt[Mcc$MONTH==5] <- "May"
Mcc$monthtxt[Mcc$MONTH==6] <- "Jun"
  
Mcc$Monthord <- factor(Mcc$monthtxt, levels = c("Feb", "Mar", "Apr", "May","Jun"))

Mcc_medians<-Mcc %>%
  filter(HR %in% c(17,18, 19, 20, 21, 22, 23, 0, 1, 2, 3,4)) %>% 
  group_by(Monthord,bunker) %>%
  summarise(MedIN=median(IN_activity))


###plot number of birds entering at the colony, with and without ships
Mcc %>%
  filter(HRordered %in% c(17,18, 19, 20, 21, 22, 23, 0, 1, 2, 3,4)) %>%
  ggplot(aes(x=HRordered, y=IN_activity, width=1))+
  geom_boxplot(colour="black") +
  facet_grid(Monthord ~ bunker,scales = "fixed", shrink = TRUE)+
  ylab("N individuals entering the cave") +
  geom_hline(data=Mcc_medians,aes(yintercept=MedIN), color='gray66',  linetype="dashed")+
  scale_x_discrete(name="Hour", labels=c(17,18,19,20,21,22,23,0,1,2,3,4))+ 
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black", family="sans"), 
        axis.title=element_text(size=22, family="sans"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)), 
        strip.text.x=element_text(size=18, color="black", family="sans"),
        strip.text.y=element_text(size=18, color="black", family="sans"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

hist(Mcc$IN_activity)

Mcc$Year <- as.factor(year(Mcc$NightStarting))
Mcc$NightStarting <- as.factor(Mcc$NightStarting)
levels(Mcc$HR)
Mcc$eventID <- as.factor(Mcc$eventID)

Q3m_nb<-glmmTMB(IN_activity~-1+HRordered+MONTH+bunker*moon.light+(1|MONTH)+(1|NightStarting/HRordered)+(1|eventID)+(1|Year), data=Mcc, family=nbinom1)
outQ3m_nb<-summary(Q3m_nb)
outQ3m_nb$coefficients
outQ3m_nb

simulationOutput <- simulateResiduals(fittedModel = Q3m_nb, plot = F)
plot(simulationOutput)
testOutliers(simulationOutput, type="bootstrap")
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

#predict and estimate percentage change caused by ship presence on colony nest attendance
Mcc$bunker <- as.factor(Mcc$bunker)
newdata1 <- expand.grid(bunker=levels(Mcc$bunker), HRordered=c(17,18,19, 20, 21, 22, 23, 0, 1, 2, 3,4), NightStarting=levels(Mcc$NightStarting), MONTH=levels(Mcc$MONTH), eventID=levels(Mcc$eventID), moon.light = seq(from= -1, to=2.0, by=0.5), Year=levels(Mcc$Year))                            
newdata<- predict(Q3m_nb, newdata=newdata1, type= "response", re.form=NULL, allow.new.levels=T)  #removing se.fit=T which needs too much memory 
newdatc <- data.frame(newdata1, newdata)
nws <- spread(newdatc, bunker, newdata)

nws <- nws %>%
  mutate(perchange=((`ships present` - `no ships`)/`no ships`)*100)#calculate % change

mean(nws$perchange)
sd(nws$perchange)

###in case 'se.fit=T' is kept in predict the simpler model below runs:

Q3m_nbr<-glmmTMB(IN_activity~-1+HRordered+MONTH+bunker*moon.light+(1|eventID)+(1|Year), data=Mcc, family=nbinom1)
outQ3m_nbr<-summary(Q3m_nbr)
outQ3m_nbr$coefficients
outQ3m_nbr

simulationOutput <- simulateResiduals(fittedModel = Q3m_nbr, plot = F)
plot(simulationOutput)
testOutliers(simulationOutput, type="bootstrap")
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

#predict % change from neg. bin. model #with (1|Nightstarting/HRordered) newdata was ~9000000 and didn't run with se.fit=T
#In original model with all random effects; variance was
#Groups                  Name        Variance  Std.Dev. 
#MONTH                   (Intercept) 3.148e-14 1.774e-07
#HRordered:NightStarting (Intercept) 4.215e-09 6.492e-05
#NightStarting           (Intercept) 6.239e-02 2.498e-01
#eventID                 (Intercept) 5.215e-02 2.284e-01
#Year                    (Intercept) 7.601e-03 8.719e-02

Mcc$bunker <- as.factor(Mcc$bunker)
newdata1 <- expand.grid(bunker=levels(Mcc$bunker), HRordered=c(17,18,19, 20, 21, 22, 23, 0, 1, 2, 3,4), MONTH=levels(Mcc$MONTH), eventID=levels(Mcc$eventID), moon.light = seq(from= -1, to=2.0, by=0.5), Year=levels(Mcc$Year))                            
newdataA<- predict(Q3m_nbr, newdata=newdata1, type= "response", re.form=NULL, allow.new.levels=T, se.fit=T)  
newdat <- data.frame(newdata1, newdataA)
nw <- newdat %>% group_by(bunker) %>% summarise(mean=mean(fit), mean(se.fit))

###check if number of ships is sig in model ####
full3 <-glmmTMB(IN_activity~HRordered+MONTH+ships+bunker*moon.light+(1|MONTH)+(1|NightStarting/HRordered)+(1|eventID), data=Mcc, family=nbinom1)  #convergence issue with (1|Year) in this model
outfull3<-summary(full3)
outfull3$coefficients
outfull3
#negative binomial glmm #(1|Year) removed to be able to compare to full model with Anova
Q3m_nb<-glmmTMB(IN_activity~HRordered+MONTH+bunker*moon.light+(1|MONTH)+(1|NightStarting/HRordered)+(1|eventID), data=Mcc, family=nbinom1)
anova(Q3m_nb, full3) #number of ships is not significant


##############################################################################################
############Q3 for entrance movements########################################################
##############################################################################################
Mcce<-merge(SQMcc,MovesENTRANCE_hr, by=c('MONTH','NightStarting','HR'), all.x=T)

Mcce$entrance_activity[is.na(Mcce$entrance_activity)]<-0

#generate bunkering event ID for use as random effect in model
Mcce <- Mcce %>%
  arrange(NightStarting) %>%
  mutate(diff = c(0,diff(NightStarting)),
         eventID = 1 + cumsum(diff > 1))

oneshipnights <- as.list(ONESHIP$NightStarting)

Mcce <- Mcce %>%
  filter(!NightStarting %in% oneshipnights)

#remove manually nights before and after 1 ship nights that were not part of longer events
Mcce<- Mcce %>%
  filter(NightStarting!="2018-03-01") %>% #2018-03-01 first night functioning after RFID failed on bunkering night 18062020
  filter(NightStarting!="2018-04-04") %>%
  filter(NightStarting!="2018-04-06") %>%
  filter(NightStarting!="2018-05-02") %>%
  filter(NightStarting!="2020-02-05") %>%
  filter(NightStarting!="2020-02-07") %>%
  filter(NightStarting!="2020-02-20") %>%
  filter(NightStarting!="2020-02-22") 

#split the last event into 2 
Mcce$eventID[1228:1267] <- 27 
#################
Mcce <- Mcce %>%
  mutate(bunker=ifelse(ships>0,'ships present','no ships')) 

Mcce$HRordered <- factor(Mcce$HR, levels = c("17","18","19", "20","21","22","23","0","1","2","3", "4","5"))


Mcce$monthtxt <- NA

Mcce$monthtxt[Mcce$MONTH==2] <- "Feb"
Mcce$monthtxt[Mcce$MONTH==3] <- "Mar"
Mcce$monthtxt[Mcce$MONTH==4] <- "Apr"
Mcce$monthtxt[Mcce$MONTH==5] <- "May"
Mcce$monthtxt[Mcce$MONTH==6] <- "Jun"

Mcce$Monthord <- factor(Mcce$monthtxt, levels = c("Feb", "Mar", "Apr", "May","Jun"))

Mcce_medians<-Mcce %>%
  filter(HR %in% c(17,18, 19, 20, 21, 22, 23, 0, 1, 2, 3,4,5)) %>% 
  group_by(Monthord,bunker) %>%
  summarise(Medentrance=median(entrance_activity))

#plot number of ind moving in entrance 
Mcce %>%
  ggplot(aes(x=ships, y=entrance_activity, width=1))+
  geom_point(aes(color=moon.light)) +
  scale_color_gradient(low="black", high="yellow", limits=c(-1,1.9), breaks=c(-1,-0.5,0,0.5,1,1.5,1.9))+
  guides(color=guide_legend(title="Moon light"), size = guide_legend(title="N. ships"))+
  #facet_wrap("Monthord", scales = "fixed", shrink = FALSE, ncol=1)+
  ylab("N individuals moving in the cave entrance") +
  xlab("N ships")+ 
  #geom_hline(data=Mcce_medians,aes(yintercept=Medentrance), color='gray66',  linetype="dashed")+
  #scale_x_discrete(name="Hour", labels=c(17,18,19,20,21,22,23,0,1,2,3,4,5))+ 
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black", family="sans"), 
        axis.title=element_text(size=22, family="sans"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)), 
        strip.text.x=element_text(size=18, color="black", family="sans"),
        strip.text.y=element_text(size=18, color="black", family="sans"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

Mcce$Year <- as.factor(year(Mcce$NightStarting))
Mcce$eventID <- as.factor(Mcce$eventID)

#neg bin glm 
Q3m_e<-glmmTMB(entrance_activity~-1+HRordered+MONTH+bunker*moon.light+(1|Year)+(1|eventID), data=Mcce,ziformula=~1, family=nbinom2)
outQ3m_e<-summary(Q3m_e)
outQ3m_e$coefficients
outQ3m_e

#neg bin glm with number of ships
Q3m_e_sh<-glmmTMB(entrance_activity~-1+HRordered+MONTH+bunker*moon.light+ships+(1|Year)+(1|eventID), data=Mcce, ziformula=~1, family=nbinom2)
outQ3m_e_sh<-summary(Q3m_e_sh)
outQ3m_e_sh$coefficients
outQ3m_e_sh

anova(Q3m_e, Q3m_e_sh) #number of ships sig

simulationOutput <- simulateResiduals(fittedModel = Q3m_e_sh, plot = F)
plot(simulationOutput)

##############################################################################################
############combine all data into one graph#################################################
##############################################################################################

Moves_Night <- MIOperall %>%
  filter(firstmove==1) %>%
   filter(direction=="IN")  %>%
   filter(error!=1)%>%
    group_by(NightStarting) %>%
    summarise(activity=length(unique(TagID)))%>%
    mutate(n_deployed=sapply(NightStarting,taggedby)) %>%
    mutate(prop_moving=(activity/n_deployed)*100)  %>%
    arrange(NightStarting)


ACTSQM <- merge(SQM, Moves_Night, by = "NightStarting", all.x=TRUE, all.y=TRUE)

ACTSQMna <- ACTSQM %>% 
  filter(is.na(activity)) %>%
  filter(NightStarting %in% monitornightsperlist) %>%
  mutate(year= year(NightStarting))

unique(ACTSQMna$NightStarting)

ACTSQM <- ACTSQM %>%
  mutate(year= year(NightStarting)) 

ACTSQM <- anti_join(ACTSQM, ACTSQMna)

ACTSQMna$activity[is.na(ACTSQMna$activity)] <- 0
ACTSQMna$prop_moving[is.na(ACTSQMna$prop_moving)] <- 0

ACTSQM <- rbind(ACTSQM, ACTSQMna)

    range(ACTSQM$candelaperm2, na.rm=TRUE)
    range(ACTSQM$activity, na.rm=TRUE)
    range(ACTSQM$prop_moving, na.rm=TRUE)
    
    ylim.prim <- c(9.581285e-05, 1.859618e-02)
    ylim.sec <- c(0, 45.8333)
    
    b <- diff(ylim.prim)/diff(ylim.sec)
    a <- b*(ylim.prim[1] - ylim.sec[1])

    shipnights<-monitornightsper %>% 
      filter(ships>0)  
    
       ACTSQM %>%
  #filter((NightStarting %within% per17x) | (NightStarting %within% per18) | (NightStarting %within% per19) | (NightStarting %within% per20x))%>%
  ggplot(aes(x=NightStarting, y=candelaperm2))+geom_point(colour="gray", size=1.5, alpha=0.5, na.rm=TRUE) +
  facet_wrap('year', scales="free_x")+
  geom_line(aes(y=a+activity*b), size=0.8, colour="black", na.rm=TRUE)+
  geom_vline(aes(xintercept=shipnights$NightStarting[1]), color='red',  linetype="dashed", size=shipnights$ships[1]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[2]), color='red',  linetype="dashed", size=shipnights$ships[2]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[3]), color='gray66',  linetype="dotdash", size=shipnights$ships[3]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[4]), color='red',  linetype="dashed", size=shipnights$ships[4]/50) + 
  geom_vline(aes(xintercept=shipnights$NightStarting[5]), color='red',  linetype="dashed", size=shipnights$ships[5]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[6]), color='red',  linetype="dashed", size=shipnights$ships[6]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[7]), color='red',  linetype="dashed", size=shipnights$ships[7]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[8]), color='red',  linetype="dashed", size=shipnights$ships[8]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[9]), color='red',  linetype="dashed", size=shipnights$ships[9]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[10]), color='red',  linetype="dashed", size=shipnights$ships[10]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[11]), color='red',  linetype="dashed", size=shipnights$ships[11]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[12]), color='red',  linetype="dashed", size=shipnights$ships[12]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[13]), color='red',  linetype="dashed", size=shipnights$ships[13]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[14]), color='gray66',  linetype="dotdash", size=shipnights$ships[14]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[15]), color='gray66',  linetype="dotdash", size=shipnights$ships[15]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[16]), color='gray66',  linetype="dotdash", size=shipnights$ships[16]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[17]), color='gray66',  linetype="dotdash", size=shipnights$ships[17]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[18]), color='red',  linetype="dashed", size=shipnights$ships[18]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[19]), color='red',  linetype="dashed", size=shipnights$ships[19]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[20]), color='red',  linetype="dashed", size=shipnights$ships[20]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[21]), color='red',  linetype="dashed", size=shipnights$ships[21]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[22]), color='red',  linetype="dashed", size=shipnights$ships[22]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[23]), color='red',  linetype="dashed", size=shipnights$ships[23]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[24]), color='gray66',  linetype="dotdash", size=shipnights$ships[24]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[25]), color='gray66',  linetype="dotdash", size=shipnights$ships[25]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[26]), color='red',  linetype="dashed", size=shipnights$ships[26]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[27]), color='red',  linetype="dashed", size=shipnights$ships[27]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[28]), color='red',  linetype="dashed", size=shipnights$ships[28]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[29]), color='red',  linetype="dashed", size=shipnights$ships[29]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[30]), color='red',  linetype="dashed", size=shipnights$ships[30]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[31]), color='red',  linetype="dashed", size=shipnights$ships[31]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[32]), color='red',  linetype="dashed", size=shipnights$ships[32]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[33]), color='red',  linetype="dashed", size=shipnights$ships[33]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[34]), color='red',  linetype="dashed", size=shipnights$ships[34]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[35]), color='red',  linetype="dashed", size=shipnights$ships[35]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[36]), color='red',  linetype="dashed", size=shipnights$ships[36]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[37]), color='red',  linetype="dashed", size=shipnights$ships[37]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[38]), color='red',  linetype="dashed", size=shipnights$ships[38]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[39]), color='red',  linetype="dashed", size=shipnights$ships[39]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[40]), color='red',  linetype="dashed", size=shipnights$ships[40]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[41]), color='red',  linetype="dashed", size=shipnights$ships[41]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[42]), color='red',  linetype="dashed", size=shipnights$ships[42]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[43]), color='red',  linetype="dashed", size=shipnights$ships[43]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[44]), color='red',  linetype="dashed", size=shipnights$ships[44]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[45]), color='red',  linetype="dashed", size=shipnights$ships[45]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[46]), color='red',  linetype="dashed", size=shipnights$ships[46]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[47]), color='red',  linetype="dashed", size=shipnights$ships[47]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[48]), color='red',  linetype="dashed", size=shipnights$ships[48]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[49]), color='red',  linetype="dashed", size=shipnights$ships[49]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[50]), color='red',  linetype="dashed", size=shipnights$ships[50]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[51]), color='red',  linetype="dashed", size=shipnights$ships[51]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[52]), color='red',  linetype="dashed", size=shipnights$ships[52]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[53]), color='red',  linetype="dashed", size=shipnights$ships[53]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[54]), color='red',  linetype="dashed", size=shipnights$ships[54]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[55]), color='red',  linetype="dashed", size=shipnights$ships[55]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[56]), color='red',  linetype="dashed", size=shipnights$ships[56]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[57]), color='red',  linetype="dashed", size=shipnights$ships[57]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[58]), color='gray66',  linetype="dotdash", size=shipnights$ships[58]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[59]), color='gray66',  linetype="dotdash", size=shipnights$ships[59]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[60]), color='red',  linetype="dashed", size=shipnights$ships[60]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[61]), color='red',  linetype="dashed", size=shipnights$ships[61]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[62]), color='red',  linetype="dashed", size=shipnights$ships[62]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[63]), color='gray66',  linetype="dotdash", size=shipnights$ships[63]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[64]), color='red',  linetype="dashed", size=shipnights$ships[64]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[65]), color='red',  linetype="dashed", size=shipnights$ships[65]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[66]), color='red',  linetype="dashed", size=shipnights$ships[66]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[67]), color='red',  linetype="dashed", size=shipnights$ships[67]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[68]), color='red',  linetype="dashed", size=shipnights$ships[68]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[69]), color='red',  linetype="dashed", size=shipnights$ships[69]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[70]), color='red',  linetype="dashed", size=shipnights$ships[70]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[71]), color='red',  linetype="dashed", size=shipnights$ships[71]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[72]), color='red',  linetype="dashed", size=shipnights$ships[72]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[73]), color='red',  linetype="dashed", size=shipnights$ships[73]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[74]), color='red',  linetype="dashed", size=shipnights$ships[74]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[75]), color='red',  linetype="dashed", size=shipnights$ships[75]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[76]), color='red',  linetype="dashed", size=shipnights$ships[76]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[77]), color='red',  linetype="dashed", size=shipnights$ships[77]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[78]), color='gray66',  linetype="dotdash", size=shipnights$ships[78]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[79]), color='red',  linetype="dashed", size=shipnights$ships[79]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[80]), color='red',  linetype="dashed", size=shipnights$ships[80]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[81]), color='gray66',  linetype="dotdash", size=shipnights$ships[81]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[82]), color='gray66',  linetype="dotdash", size=shipnights$ships[82]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[83]), color='gray66',  linetype="dotdash", size=shipnights$ships[83]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[84]), color='red',  linetype="dashed", size=shipnights$ships[84]/50) +
  geom_vline(aes(xintercept=shipnights$NightStarting[85]), color='red',  linetype="dashed", size=shipnights$ships[85]/50) +
      xlab("Date")+
  scale_y_continuous(bquote('Cliff face brightness '~(cd/m^2)), 
                     sec.axis=sec_axis(~ (. - a)/b, name="Prop. of tagged individual shearwaters entering colony", breaks=seq(0,60,10), labels=as.character(seq(0,60,10))))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        text=element_text(family="sans"),
        axis.text.x = element_text(size=18, color="black", family="sans"),
        axis.text.y = element_text(size=14, color="black", family="sans"),
        axis.title=element_text(size=22, family="sans"), 
        strip.text.x=element_text(size=18, color="black", family="sans"), 
        strip.background=element_rect(fill="white", colour="black"),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.y.right=element_text(margin=margin(0,0,0,10)),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
    
#######################################################################################################################################
###Save workspace for replication and upload to github ###############
#######################################################################################################################################   
      
#IMP. to conform with MarineTraffic data policy dataframes containing ship names and call signs are removed from workspace
       
    rm(AIS)
    rm(AIS17)
    rm(AIS18)
    rm(AIS19)
    rm(AIS20) 
    rm(ships2)
    rm(ship)
    rm(ships)
    rm(S)
    rm(p1)
    rm(p2)
    rm(A)
    rm(Moves_Night)
    rm(nights)
    
    #However remaining dataframes contain the required summaries of number of ships per relevant timeframes
    #for any questions contact martinaustad93@gmail.com
    
    #reduce the size of workspace to upload to github
  rm(raw) #raw rfid data uploaded in separate workspace file: '2022_Shearwaters_RFID_Malta_Rawdata.RDATA'
  rm(RFIDMETA)
  rm(RFIDAMPS)
  rm(RFIDoper)
  
  rm(ACTSQM) 
  rm(ACTSQMna)
  rm(fonts) 
  rm(monitornightsperlist)
  
  
  
  #other intervals & dataframes can be derived from workspace:
    save.image("~/PHD/RFID_Analysis_PHD/2022_LightPollShips_Shearwaters_Malta_Rworkspace.RData")
    
