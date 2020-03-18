filename <- "MS1003_storyline.gpx"

#decimal to round coordinates to#
decimal <- 4

  #scan file#
  n <- scan(filename, character(0), sep = "\n")
  
  #create vectors for data#
  LonVec <- rep(NA,100000)
  LatVec <- rep(NA,100000)
  DateVec <- rep(NA,100000)
  TimeVec <- rep(NA,100000)
  DateTimeVec <- rep(NA,100000)
  TypeVec <- rep(NA,100000)
  c <- 1  # index for parsing of GPS code
  
  #find data in file and place it in vectors#
  for(i in 1:length(n)){
    if(grepl("<trkpt ",n[i])){
      
      lonStart <- regexpr('lon="',n[i])[1]
      lonEnd <- regexpr('" ',n[i])[1]
      x <- round(as.numeric(substr(n[i],lonStart+5,lonEnd-1)), digits = decimal)
      
      latStart <- regexpr('lat="',n[i])[1]
      latEnd <- regexpr('">',n[i])[1]
      y <- round(as.numeric(substr(n[i],latStart+5,latEnd-1)), digits = decimal)
      
      dateStart <- regexpr("<time>",n[i+1])[1]
      timeEnd <- regexpr("</time>",n[i+1])[1]
      dateEnd <- regexpr("T",n[i+1])[1]
      dataDate <- substr(n[i+1],dateStart+6,dateEnd-1)
      dataTime <- substr(n[i+1],dateStart+17,timeEnd-1)
      dataDateTime <- substr(n[i+1],dateStart+6,timeEnd-1)
      
      typeEnd <- regexpr("</type>",n[i+2])[1]
      dataType <- substr(n[i+2],15,(typeEnd-1))
      
      #add information to vectors#
        LonVec[c] <- x
        LatVec[c] <- y
        TypeVec[c] <- dataType
        DateVec[c] <- dataDate
        TimeVec[c] <- dataTime
        DateTimeVec[c] <- dataDateTime
        c <- c+1
    }
  }
  LonVec <- LonVec[!is.na(LonVec)]
  LatVec <- LatVec[!is.na(LatVec)]
  DateVec <- DateVec[!is.na(DateVec)]
  TimeVec <- TimeVec[!is.na(TimeVec)]
  DateTimeVec <- DateTimeVec[!is.na(DateTimeVec)]
  TypeVec <- TypeVec[!is.na(TypeVec)]
  
  
  #create matrix of all coordinates#
  positionMatrix = matrix(
    c(LonVec,LatVec),
    nrow=length(LonVec),
    ncol=2)
  
  #calculate timeframe#
  startDate <- DateVec[1]
  endDate <- DateVec[length(DateVec)]
  totalDays <- as.Date(endDate)-as.Date(startDate)+1
  dayCheck <- 0
  
  #calculate total number of minutes during the time of data collection#
  totalMinutes <- as.numeric(totalDays*1440)
  
  #find coordinates with duplicate minutes (exclude place trackpoints)#
  uniqueMinutes <- 0
  for(k in 1:length(TimeVec)){
    if(length(TimeVec) > 1&&grepl(substr(TimeVec[k],1,5),substr(TimeVec[k+1],1,5))&&!grepl("place",TypeVec[k])){
      positionMatrix[k,1]<-"remove"
      positionMatrix[k,2]<-"remove"
      DateTimeVec[k] <- "remove"
      DateVec[k] <- "remove"
      TypeVec[k] <- "remove"
    }
    else
      uniqueMinutes <- uniqueMinutes+1
  }
  
  
  #remove coordinates with duplicate minutes or thrown out#
  positionMatrix <- positionMatrix[positionMatrix[,1]!="remove",]
  DateTimeVec <-DateTimeVec[DateTimeVec!="remove"]
  DateVec <- DateVec[DateVec!="remove"]
  TypeVec <- TypeVec[TypeVec!="remove"]
  
  #create data frame ubdex variable#
  RE <- character(totalDays)
  dates <- character(totalDays)
  e <- 1 
  
  #calculate entropy for each day#
  date1 <- startDate
  for(b in 1:totalDays) {
    
    #initialize day matrix #
    dmsize = 1440
    dayMatrix = matrix(
      nrow=dmsize,
      ncol=2)
    
    #find data points with given date#
    dataLocation <- character(2000) # preallocating for max possible number of trackpoints in a day
    f <- 1
    for (t in 1:length(DateVec)){
      if(length(DateVec) > 0&&DateVec[t]==date1){
        dataLocation[f] <- t
        f <- f +1
      }
      if(length(DateVec) > 0&&DateVec[t]>date1)
        break
    }
    dataLocation <- dataLocation[!dataLocation %in% ""]
    
    beginDate <- as.POSIXct(as.character(date1))
    
    #place coordinates in proper timeline (each row of a matrix represents a minute of the timeframe of collected data)#
    if(length(dataLocation) >= 1)
    {
      
      #check for daylight savings#
      if(difftime(as.POSIXct(DateTimeVec[as.numeric(dataLocation[length(dataLocation)])],format = "%Y-%m-%dT%H:%M:%S"),as.POSIXct(as.character(date1)),units = "mins")>1440){
      dmsize <- 1500
      dayMatrix = matrix(
        nrow=dmsize,
        ncol=2
      )
      }
      
      #check for first day#
      if(as.Date(date1)==startDate){
        dmsize <- as.numeric(difftime(as.POSIXct(as.character((as.Date(date1)+1))),as.POSIXct(DateTimeVec[as.numeric(dataLocation[1])],format = "%Y-%m-%dT%H:%M"), units = "mins"))+1
        dayMatrix = matrix(
          nrow=dmsize,
          ncol=2
        )
        beginDate <- as.POSIXct(DateTimeVec[as.numeric(dataLocation[1])],format = "%Y-%m-%dT%H:%M")
      }
      
      startPoint <- 1
      for(u in 1:length(dataLocation)){
        tTest <- as.POSIXct(DateTimeVec[as.numeric(dataLocation[u])],format = "%Y-%m-%dT%H:%M:%S")
        timeInMinutes <- floor(as.numeric(difftime(tTest,beginDate,units = "mins")))+1
        
        if(startPoint==1){
          for(v in 1:(timeInMinutes-1)){
            dayMatrix[v,] <- positionMatrix[as.numeric(dataLocation[u]),]
          }
        }
        
        if(u==length(dataLocation)&&startPoint<=dmsize){
          for(w in startPoint:dmsize){
            dayMatrix[w,] <- positionMatrix[as.numeric(dataLocation[u]),]
          }
        }
        
        if(timeInMinutes-startPoint>0&&u!=1&&u!=length(dataLocation)){
          
          previous_lat <- as.numeric(positionMatrix[as.numeric(dataLocation[u])-1,1])
          current_lat <- as.numeric(positionMatrix[as.numeric(dataLocation[u]),1])
          previous_lon <- as.numeric(positionMatrix[as.numeric(dataLocation[u])-1,2])
          current_lon <- as.numeric(positionMatrix[as.numeric(dataLocation[u]),2])
          previous_time <- startPoint-1
          current_time <- timeInMinutes
          avg_lat_velocity <- (current_lat-previous_lat)/(current_time-previous_time)
          avg_lon_velocity <- (current_lon-previous_lon)/(current_time-previous_time)
          
          for(z in startPoint:(timeInMinutes-1)){
            dayMatrix[z,1] <- round((avg_lat_velocity*(z-previous_time)+previous_lat), digits = decimal)
            dayMatrix[z,2] <- round((avg_lon_velocity*(z-previous_time)+previous_lon), digits = decimal)
          }
        }
        
        dayMatrix[timeInMinutes,] <- positionMatrix[as.numeric(dataLocation[u]),]
        startPoint <- timeInMinutes+1  
      } 
    }
    
    #convert matrix to dataframe#
    df <- as.data.frame(dayMatrix,row.names=NULL)
    colnames(df) <- c("lon","lat")
    
    #calculate roaming entropy for that date and place it in RE vector#
    library(plyr)
    uniquePositions <- ddply(df,.(lon,lat),nrow)
    uniquePositions$p <- uniquePositions$V1/dmsize
    uniquePositions$plogp <- uniquePositions$p*log2(uniquePositions$p)
    sampleRE = (-sum(uniquePositions$plogp))/log(648000000)
    if(length(dataLocation)==0)
      sampleRE <- NA
    RE[e] <- sampleRE
    dates[e] <- date1
    e <- e+1
    date1 <- as.character(as.Date(date1)+1)
  }

re_matrix <- matrix(nrow = 0, ncol=0)
re_matrix <- cbind(dates,RE)
colnames(re_matrix) <- c("date","roaming entropy")

writeoutfilename <- paste0(substr(filename,1,nchar(filename)-4),"_re_matrix.csv")
write.csv(re_matrix,file = writeoutfilename, row.names = FALSE)