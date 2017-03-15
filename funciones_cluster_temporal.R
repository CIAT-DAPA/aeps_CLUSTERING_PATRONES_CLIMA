
#rm(list = ls())



###### Functions

#Procesar datos d

unifDatos <- function(listaInfo,namesC)
{
    FECHA          <- listaInfo[[1]][,1]
    listProces     <- lapply(listaInfo,function(x){x[,2]})  
    listCompleta   <- do.call(cbind,listProces)
    dfListCompleta <- data.frame(FECHA,as.data.frame(listCompleta))
    names(dfListCompleta) <- c("DATE",namesC)
    return(dfListCompleta)
}


joinEventsClim <- function(climStat,cropData,datCS = "DATE", sowDat= "sowDate", harvDat ="harvestDate",formDate = "%d/%m/%Y",delEvents=F)
{
    climStat[,datCS]   <- as.Date(climStat[,datCS],format = formDate)
    cropData[,sowDat]  <- as.Date(cropData[,sowDat],format = formDate)
    cropData[,harvDat] <- as.Date(cropData[,harvDat],format = formDate)
    
    if( min(cropData[,sowDat]) < min(climStat[,datCS]) |  max(cropData[,harvDat]) > max(climStat[,datCS])  )warning("There are events with dates out of station dates")
    
    
    if(delEvents){warning("these will be deleted");cropData=cropData[cropData[,sowDat] %in% climStat[,datCS] & cropData[,harvDat] %in% climStat[,datCS],]} ##### Eliminar los eventos que no estan en el rango de tiempo
    lotNam <- row.names(cropData)
    cropDates <- lapply(1:nrow(cropData),function(x){
        EVENT <- lotNam[x];
        DATE <- as.Date(cropData[,sowDat][x]:cropData[,harvDat][x],origin="1970-01-01"); 
        eventDate <- data.frame(EVENT,DATE); 
        merge(eventDate,climStat,by.x = "DATE" , by.y = datCS, all.x= T, all.y = F,sort = T) })
    do.call(rbind,cropDates)[,c(2,1,3:ncol(cropDates[[1]]))]
}



GroupingPeriod= function(climEvent,idVar = "EVENT",datVar="Date",vars = "ESOL",LimitDays=F,MaxDaysProcess=127, Period=1,VarAccum="RAIN")
{
    
    if(Period==1){
      return(climEvent)
    }
    else{
      return(data.frame(do.call(rbind,lapply(split(climEvent,climEvent[,idVar]),function(x){
        if(LimitDays){x=x[1:ifelse(MaxDaysProcess>=nrow(x),nrow(x),MaxDaysProcess),]}
        if(nrow(x)%%Period/Period>0.8){dif=0}else{dif=nrow(x)%%Period}
        x=x[1:(nrow(x)-dif),]
        x$Group=ceiling((1:nrow(x))/Period)
        varMeans=vars[!vars%in%VarAccum]
        varAccum=vars[vars%in%VarAccum]
        auxDF=data.frame(do.call(rbind,lapply(split(x,x$Group),function(y){
          Means=data.frame(sapply(y[varMeans], mean))
          Sums=data.frame(sapply(y[varAccum], sum))
          return(c(as.character(unique(x[,idVar])),as.numeric(Means[,1]),as.numeric(Sums[,1])))
        })))
        names(auxDF)=c(idVar,varMeans,varAccum)
        auxDF=auxDF[,c(idVar,vars)]
        auxDF[,vars]=sapply(auxDF[,vars], function(y){as.numeric(as.character(y))})
        return(auxDF)
      }))))
    }
}

procesData <- function(climEvent,idVar = "ID",NormMethod=1,vars="ESOL"){
    
    require(gtools)
    
    climaEventos0 <-  climEvent
    
    climaEventos0=climaEventos0[,c(idVar,vars)]
    
    #NORMALIZACION
    
    climDataNorm0<-  as.data.frame(do.call(cbind,lapply(vars,function(x){
        y <- climaEventos0[x][,1];
        if(NormMethod==1){return((y-mean(y,na.rm = T))/sd(y,na.rm = T))}else{return((y - min(y,na.rm = T))/(max(y,na.rm = T)-min(y, na.rm = T)))}})))
    
    colnames(climDataNorm0) <-  vars
    
    climDataNorm      <- climaEventos0
    climDataNorm[vars] <- climDataNorm0
    
    tempM <- colMeans(climaEventos0[vars])
    sdM   <- apply(climaEventos0[vars],2,sd)
    
    parClim     <- data.frame(rbind(tempM,sdM ))
    
    nlevents <- lapply(split(climaEventos0,climaEventos0[,idVar]),function(x){x$day=1:nrow(x);x=x[c(idVar,vars,"day")];return(x)})
    nleventsN <- lapply(split(climDataNorm,climDataNorm[,idVar]),function(x){x$day=1:nrow(x);x=x[c(idVar,vars,"day")];return(x)})
    
    allEvents <- do.call(rbind,nlevents)
    allEventsN <- do.call(rbind,nleventsN)
    
    sorData <- reshape(allEvents,idvar =idVar,timevar = "day", direction = "wide")
    sorDataN <- reshape(allEventsN,idvar =idVar,timevar = "day", direction = "wide")
    
    namSortData <- names(sorData)[-1]
    namSortDataN <- names(sorDataN)[-1]
    
    sorDataF=sorData[,c(idVar,mixedsort(namSortData))]
    sorDataFN <- sorDataN[,c(idVar,mixedsort(namSortDataN))]
    events=unique(climaEventos0[,idVar])
    row.names(sorDataF) <- events
    row.names(sorDataFN) <- events
    
    #write.csv(sorDataF,"climateToCluster.csv")
    #names(nlevents)=unique(climaEventos0[,idVar])
    evenF <- lapply(nlevents,function(x){x[,vars]})
    evenN <- lapply(nleventsN,function(x){x[,vars]})
    
    save(evenF,events,evenN,parClim,file = "listClimatEvent.RData")
    
}

#distMatrix <- distAllMatrix

hirarCluster <- function(distMatrix)
{
    require(dtwclust)


    ClusterEvents      <- hclust(distMatrix,method="average")

    upDate="Y"
    save(ClusterEvents,file="ClusterEvents.RData")
 
    
    #Condicional
    #disClust
    #totLines <- readLines(n = 1)
    
    
    while(upDate!="N")
    {
        print("Choose the maximun number of cluster to expand the graphic:")
        maxNumb <- readLines(n = 1)
        
        barplot(sort(ClusterEvents$height,decreasing=T)[1:maxNumb],names.arg=1:maxNumb,col="lightseagreen",main="Select the number of Cluster",xlab="Number of cluster",ylab="height") 
        abline(h=0)
        
        print("Do you want update the barplot? Y/N")
        upDate = readLines(n = 1)
    }
    
    tiff("barplotGraph.tiff",width = 740,height = 400)   
    barplot(sort(ClusterEvents$height,decreasing=T)[1:maxNumb],names.arg=1:maxNumb,col="lightseagreen",main="Select the number of Cluster",xlab="Number of cluster",ylab="height") 
    abline(h=0)  
    dev.off()
    
    print("Number of cluster:")
    
    nClust <- readLines(n = 1)
    
    membAll <- cutree(ClusterEvents,k=nClust)

    return(membAll)
}


distDtwMV <- function(listObje)
{
    require(dtw)
    len      <-  length(listObje)
    listDist <- matrix(0,nrow=len,ncol=len)
    
    disDtw <- array(0,len)
    # pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)
    pb <- txtProgressBar(min = 0, max = len, style = 3)
    for(i in 1:len){ 
        
        for(j in i:len){
            listDist[j,i]<- dtw(listObje[[i]],listObje[[j]])$distance
        }
        setTxtProgressBar(pb, i)
        #  info <- sprintf("%d%% done", round((i/(len)*100)))
        # setWinProgressBar(pb, i/(len)*100, label=info)
    }
    close(pb)
    #close(pb)
    rownames <- names(listObje)
    colnames <- names(listObje)
    return(as.dist(listDist))
}



multiplotClust= function(listG=subBas,NomVarx,NomVary,limts,textNum=FALSE,mainTitle=paste0("Cluster_",i),SaveMean=F){
    numVarx=which(grepl(NomVarx,names(listG[[1]])))
    plot(listG[[1]][,numVarx],listG[[1]][[NomVary]],xlim = c(listG[[1]][,numVarx][1],listG[[1]][,numVarx][nrow(listG[[1]])]),ylim=c(limts[1],limts[2]),col=0,ylab=NomVary,xlab=NomVarx,main=mainTitle,cex.axis=1,cex.lab=1)
    lapply(listG,function(x){lines(x[,numVarx],x[[NomVary]],col="gray")})
    if(textNum)text(x = listG[[1]][,numVarx][1],y = limts[1],labels =  paste0("Num. Ind.: ",length(listG)),col = 4,cex.lab=1,cex = 1)
    limitX=lapply(listG, nrow)
    limitXmin=min(do.call(cbind,limitX))
    Xaxis=apply(do.call(cbind,lapply(listG,function(x){x[,numVarx][1:limitXmin]})),1,mean,na.rm=T)
    Yaxix=apply(do.call(cbind,lapply(listG,function(x){x[[NomVary]][1:limitXmin]})),1,mean,na.rm=T)
    lines(Xaxis,Yaxix,col="red",lwd=2,cex.axis=1,cex.lab=1,ylab="")
    if(SaveMean){return(as.data.frame(cbind(Xaxis,Yaxix)))}
}
##############################################################
rsmnClustClima <- function(evenF,eventClasf,vars){
    
    if(!file.exists("AllCluster")){dir.create("AllCluster")}else{}
    
    listEventsDate=lapply(evenF,function(x){return(as.data.frame(cbind(
        Dates=rep(1:nrow(x)),x)))})
    
    membAll=eventClasf
    
    table(membAll)
    
    nlevents=listEventsDate
    
    
    minVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x,2,
                                                                     min,na.rm = T)})),2,min,na.rm = T)
    
    maxVal  <- apply(do.call(rbind,lapply(nlevents,function(x){apply(x,2,max,
                                                                     na.rm = T)})),2,max,na.rm = T)
    
    limts <- as.data.frame(rbind(minVal,maxVal))
    
    NomVarxAxix="Dates"
    nomVarClust=vars
    
    
    listMeans=list()
    
    for(i in 1:length(unique(membAll)))
    {
        tiff(paste0("AllCluster/DTWcluster",i,"_Line.tiff"),width = 8, 
             height = 12,res=200,units = 'in')
        
        layout(cbind(1:length(nomVarClust)))
        subBas <- nlevents[which(membAll==i)]
        
        
        for(j in nomVarClust)
        {
            limtsV=limts[[j]]
            listMeans=multiplotClust(listG = subBas,NomVarx=NomVarxAxix,
                                     NomVary = j,limts =limtsV,textNum = T,
                                     mainTitle = paste0("Cluster_",i),SaveMean = F)##Activar para guardar 
        }
        dev.off()
    }
}


whisk <- function(df,cond_col=1,val_col=2) {
    require(reshape2)
    condname <- names(df)[cond_col]
    names(df)[cond_col] <- "labGr" 
    names(df)[val_col] <- "Yield"
    b <- boxplot(Yield~labGr,data=df,plot=FALSE)
    df2 <- data.frame(b$stats[,1:ncol(b$stats)],c("min","lq","m","uq","max"))
    names(df2) <- c(levels(df$labGr),"pos")
    df2 <- melt(df2,id="pos",variable.name="labGr")
    df2 <- dcast(df2,labGr~pos)  
    names(df2)[1] <- condname
    df2
}

#baseConComparacion <- basePrueba
# vary <- "RENDIMIENTO_kg.ha"
#varx <- "Variedad" 

krusk.boxplot <- function(baseConComparacion,vary,varx,ylabs=expression(paste("Rice Yield (kg.",ha^-1,")")),xlabs="Cluster - (Number of observation)",maxVar=1,title="")
{  
    conteo <- table(baseConComparacion[,varx])
    
    #Filtrar por valores mayores a 10
    
    baseConComparacion <- baseConComparacion[baseConComparacion[,varx] %in% names(conteo)[conteo>maxVar],]
    
    baseConComparacion <- droplevels(baseConComparacion)
    
    summary(baseConComparacion)
    
    #Extraer y ordernar labels
    
    orderLabels <- with(baseConComparacion,tapply(get(vary),get(varx),median))
    
    nPerClust <- with(baseConComparacion,tapply(get(vary),get(varx),length))
    
    maxValue <- with(baseConComparacion,tapply(get(vary),get(varx),max))
    
    medianCount <- data.frame(varCon=names(orderLabels),
                              orderLabels,nPerClust,maxValue,labGr=paste(names(orderLabels),
                                                                         "-","(",nPerClust,")",sep=""))
    
    resumenDatos <- medianCount[order(medianCount$orderLabels,decreasing = T), ]
    
    labsGr <- resumenDatos$labGr
    
    baseConComparacionN <-merge(baseConComparacion,medianCount[c("varCon","labGr")],by.x=varx,by.y="varCon",all.x=T,all.y=F,sort = F)
    
    baseConComparacionN$labGr <- factor(baseConComparacionN$labGr,levels = labsGr)
    
    #Prueba de Kruskal
    
    kruskalTest <- kruskal(y=baseConComparacionN[,vary],trt=baseConComparacionN$labGr,group=T)$groups
    
    #Grafico
    
    myplot <- ggplot(baseConComparacionN,aes(x=labGr))+ theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        xlab(xlabs)+ylab(ylabs)+ 
        geom_errorbar(aes(ymin=min,ymax=max),data=whisk(baseConComparacionN[c("labGr",vary)]),width = 0.2)+geom_boxplot(aes(y=get(vary)))+
        annotate("text",x=1:nrow(kruskalTest),y=resumenDatos$maxValue*1.04,label=kruskalTest$M)+scale_y_continuous(breaks = seq(3000,12000,1000))+
        ggtitle(title)
    
    
    #tiff("//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/PUBLICACIONES/GRAFICOS_PUBLICACION_RESCATE/GRAFICOS_RESULTADOS/boxplot_yield_cluster.tiff",width  = 6,height = 4.2,res=200,units='in')
    print(myplot)
    # dev.off()
}



