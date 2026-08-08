#-----------------------------------------------------
#

arma_base<-function(video,cut_above = 0,cut_below = 0,nombres = rois){
  path <- paste("txts/",as.character(video),sep = "")
  x <- read.delim(path,header = FALSE, sep = " ")
  if (dim(x)[2]==1){
    x <- read.delim(path,header = FALSE, sep ="\t")
  }
  x <-  select(x,1,2,3,5,6,7)
  t <- 1:dim(x)[1]
  x <- cbind(t,x)
  ca <- cut_above + 1
  cb <- dim(x)[1]-cut_below
  x <- x[ca:cb, ]
  names(x)<-c("frame",nombres)
  return(x)
}


#---------------

crea_lista_movs_ <- function(vp, q = 0.71, nombre_video = NA, sr = 30, lagmax = 150) {
  if (!exists("zones")) stop("Objeto global 'zones' no definido.")
  if (!exists("logbook")) stop("Objeto global 'logbook' no definido.")

  resultados <- list()
  idx <- 1

  for (zona in zones) {
    # smoothing the MEA signals for patient and therapist
    temp1 <- lapply(vp, function(df) suavizar(
      as.numeric(df[, grep(paste0(zona, ".*patient"), names(df))]), sr))
    temp2 <- lapply(vp, function(df) suavizar(
      as.numeric(df[, grep(paste0(zona, ".*therapist"), names(df))]), sr))

    # detect segments with high movement
    altos1 <- detecta_altos(temp1, q)
    altos2 <- detecta_altos(temp2, q)
    comunes <- intersect(altos1, altos2)

    for (j in comunes) {
      df <- vp[[j]]
      mov_p <- temp1[[j]]
      mov_t <- temp2[[j]]
      start_time <- round(df$frame[1] / (60 * sr), 1)
      end_time <- round(df$frame[nrow(df)] / (60 * sr), 1)

      # Spearman cross-correlation
      ccf_s <- ccf(rank(mov_p), rank(mov_t), lag.max = lagmax, plot = FALSE)
      y <- ccf_s$acf[, 1, 1]
      lag_vals <- ccf_s$lag[, 1, 1]

      max_i <- which.max(y)
      min_i <- which.min(y)

      spearman_max <- y[max_i]
      spearman_min <- y[min_i]
      lag_max <- round(lag_vals[max_i] / sr, 2)
      lag_min <- round(lag_vals[min_i] / sr, 2)

      # get sex and type info from logbook
      x <- which(logbook$file == nombre_video)
      sexo <- ifelse(length(x) > 0, logbook$sex[x], NA)
      tipo <- ifelse(length(x) > 0, logbook$classification[x], NA)

      resultados[[idx]] <- data.frame(
        video = nombre_video,
        zone = zona,
        num_period = j,
        start = start_time,
        end = end_time,
        mov_patient = mean(mov_p),
        mov_therapist = mean(mov_t),
        mov_medio = mean(c(mean(mov_p), mean(mov_t))),
        spearman_max = spearman_max,
        lag_spearman_max = lag_max,
        spearman_min = spearman_min,
        lag_spearman_min = lag_min,
        #leads = lider,
        sex = sexo,
        type = tipo,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  final_df <- do.call(rbind, resultados)
  return(final_df)
}

crea_lista_movs <- function(vp,q=0.71,nombre_video=NA,sr=30,lagmax=150){
  if (!exists("zones")) stop("Objeto global 'zones' no definido.")
  if (!exists("logbook")) stop("Objeto global 'logbook' no definido.")

  selectos <-data.frame(matrix(ncol = 15,nrow = 0))
  colnames(selectos) <- c("video","zone","num_period","start","end",
                          "mov_patient","mov_therapist","mov_medio",
                          "spearman_max","lag_spearman_max",
                          "spearman_min","lag_spearman_min",
                          "leads",
                          "sex","type")

  for (k in zones){
    temp<-list()
    temp1<-list()
    temp2<-list()
    for (i in 1:length(vp)){
      temp[[i]] <- vp[[i]][ ,grep(k,names(vp[[i]]))]
      temp1[[i]] <- suavizar(
        as.numeric(temp[[i]][ ,grep("patient",names(temp[[i]]))]),30)
      temp2[[i]] <- suavizar(
        as.numeric(temp[[i]][ ,grep("therapist",names(temp[[i]]))]),30)
    }
    aux1 <- detecta_altos(temp1,q)
    aux2 <- detecta_altos(temp2,q)
    Aux<-intersect(aux1,aux2)
    for (j in Aux){
      renglon <- data.frame(matrix(ncol=ncol(selectos),nrow=1))
      names(renglon)<-names(selectos)
      renglon$video <- nombre_video
      renglon$zone <- k
      renglon$num_period <- j
      renglon$start <- round(vp[[j]]$frame[1]/(60*sr),1)
      renglon$end <- round(vp[[j]]$frame[nrow(vp[[j]])]/(60*sr),1)
      renglon$mov_patient <- mean(temp1[[j]])
      renglon$mov_therapist <- mean(temp2[[j]])
      renglon<-mutate(renglon,mov_medio=mean(c(mov_patient,mov_therapist)))
      cecefe_s <- ccf(rank(temp1[[j]]),rank(temp2[[j]]),lag.max = lagmax,plot=FALSE)

      y<-cecefe_s$acf[ ,1,1]
      x1<-which(y==max(y))
      x2<-which(y==min(y))

      renglon$spearman_max <- y[x1]
      renglon$spearman_min <- y[x2]
      renglon$lag_spearman_max <- round(cecefe_s$lag[x1,1,1]/sr,2)
      renglon$lag_spearman_min <- round(cecefe_s$lag[x2,1,1]/sr,2)
      renglon$leads<-"-"
      signo <- renglon$spearman_max-abs(renglon$spearman_min)
      if (signo>0 & renglon$lag_spearman_max>0.25){
        renglon$leads="therapist"
      }
      if (signo>0 & renglon$lag_spearman_max< -0.25 ){
        renglon$leads="patient"
      }
      if (signo<0 & renglon$lag_spearman_min>0.25){
        renglon$leads="therapist"
      }
      if (signo<0 & renglon$lag_spearman_min< -0.25 ){
        renglon$leads="patient"
      }
      x <- which(logbook$file == renglon$video)
      renglon$sex <- logbook$sex[x]
      renglon$type <- logbook$classification[x]
      selectos<-rbind(selectos,renglon)
    }
  }
  return(selectos)
}


#------------------------------

crear_lista <- function(video_partido, nombre=NA, sr=30, segundos = 300,
                        umbral_max=0.7,umbral_min=-0.7,umbral_mean=0.35,
                        lag_mayor = 150){
  selectos <- data.frame(matrix(ncol = 11,nrow = 0))
  colnames(selectos) <- c("video","zone","num_period","start","end",
                          #"acf_max","acf_promedio","lag_acf_max",
                          "spearman_max","lag_spearman_max",
                          "spearman_min","lag_spearman_min","sex","type")

  for (i in 1:length(video_partido)){
    for (j in zones){
      temp <- video_partido[[i]][ ,grep(j,names(video_partido[[i]]))]
      suave<-data.frame(suavizar(temp[ ,1],30),suavizar(temp[ ,2],30))
      names(suave)<-names(temp)
      temp<-suave
      if(sum(abs(temp[1]))>1000 && sum(abs(temp[2]))>1000){
        cecefe_s <- ccf(rank(temp[1]),rank(temp[2]),plot = FALSE,lag.max = 150,method="spearman")
        if (max(abs(cecefe_s$acf)>umbral_max) | abs(mean(cecefe_s$acf))>umbral_mean){
          cecefe_s<-ccf(rank(temp[1]),rank(temp[2]),lag.max = lag_mayor,
                        plot=FALSE,method="spearman")
          y<-cecefe_s$acf[ ,1,1]
          x1<-which(y==max(y))
          x2<-which(y==min(y))
          aux <- data.frame(matrix(ncol=ncol(selectos),nrow=1))
          names(aux)<-names(selectos)
          aux$video = nombre
          aux$zone = j
          aux$num_period <- i
          aux$start <- round(video_partido[[i]]$tiempo[1]/(sr*60),1)
          aux$end <- aux$start + round(segundos/60,2)
          aux$spearman_max <- max(cecefe_s$acf)
          aux$lag_spearman_max <- round(cecefe_s$lag[x1,1,1]/30,2)
          aux$spearman_min <- min(cecefe_s$acf)
          aux$lag_spearman_min <- round(cecefe_s$lag[x2,1,1]/30,2)
          x <- which(logbook$file == aux$video)
          aux$sex <- logbook$sex[x]
          aux$type <- logbook$classification[x]
          selectos<-rbind(selectos,aux)
        }
      }
      else{
      }
    }
  }
  return(selectos)
}


#-----------------

detecta_altos <- function(lista_datos,q=0.9){
  promedios <- numeric()
  for (i in 1:length(lista_datos)){
    promedios <- c(promedios,mean(lista_datos[[i]],na.rm=TRUE))
  }
  altos <- which(promedios>quantile(promedios,probs=q))
  return(altos)
}

#-----------------------------------------------------


partition_data <- function(datos,minutes=5,sr=30,pasos=1){
  partes <- list()
  largo_bloques <- minutes*sr*60
  num_partes <- pasos*floor(nrow(datos)/largo_bloques)
  for (i in 1:num_partes){
    inicia <- (i-1)*floor(largo_bloques/pasos)+1
    termina <- inicia + largo_bloques -1
    partes[[i]] <- datos[inicia:termina, ]
  }
  return(partes)
}

#-----------------------------------------------------
#

suavizar <- function(listado,n=7){
  aux<-listado
  for (j in n:length(listado)){
    aux[j] <- mean(listado[(j-n+1):j],na.rm = TRUE)
  }
  for (j in 1:(n-1)){
    aux[j] <- mean(listado[j:(j+n-1)],na.rm = TRUE)
  }
  return(aux)
}



#--------------------GENERATION OF VISUALS -----------------

corta_lista <- function(id,ls=selected_lapses){
  x <- which(ls$id==id)
  delta <- ls$end[x]-ls$start[x]
  video<-ls$video[x]
  period<-ls$num_period[x]
  aux <- partition_data(completa[[video]],minutes=delta)
  aux<-aux[[period]]
  aux <- aux[ ,grep(ls$zone[x],names(aux))]
  aux$suave_patient <- suavizar(aux[ ,1],30)
  aux$suave_therapist <- suavizar(aux[ ,2],30)
  aux$tiempo <- (1:dim(aux)[1]/(60*30))+ls$start[x]
  temp <- list()
  temp[[1]]<-aux
  temp[["video"]]<-paste("Video: ",ls$video[x])
  temp[["minutes"]]<-paste("From ",as.character(ls$start[x])," to ",
                           as.character(ls$end[x]),"min.")
  temp[["roi"]]<-paste("ROI: ",ls$zone[x])
  return(temp)
}

comparativa_visual <- function(id,min_ini=0,min_fin=100,que_lista=selected_lapses){
  temp <- corta_lista(id,ls=que_lista)[[1]]
  temp2 <- corta_lista(id,ls=que_lista)
  g<-ggplot()+geom_line(data=filter(temp,tiempo>min_ini,tiempo<min_fin),aes(tiempo,suave_patient),
                        color='blue',size=0.5)+geom_line(data=filter(temp,tiempo>min_ini,tiempo<min_fin),
                                                         aes(tiempo,suave_therapist),color='brown',size=0.5)+
    ggtitle(temp2[["video"]])+xlab("minute in video")+
    ylab("Blue: patient. Red: therapist.")
  show(g)
  return(g)
}

corr_cruzadas <- function(id,que_lista=selected_lapses){
  temp <- corta_lista(id,ls=que_lista)[[1]]
  cecefe <- ccf(rank(temp$suave_patient),rank(temp$suave_therapist),
                plot = FALSE,lag.max = 450)
  aux <- data.frame(segundos=cecefe$lag[ ,1,1]/30,correlacion=cecefe$acf[,1,1])
  temp2 <- corta_lista(id,ls=que_lista)
  g <- ggplot(aux,aes(segundos,correlacion))+geom_col(color='purple',fill='purple')+
    ggtitle(paste(temp2[["video"]],temp2[["minutes"]],temp2[["roi"]],sep = "                 "))+
    xlab("patient leads      <--   lag in seconds   -->  therapist leads")+
    ylab("cross correlation")

  show(g)
  return(g)
}

