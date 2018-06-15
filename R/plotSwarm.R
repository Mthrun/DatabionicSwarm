plotSwarm=function(Points,Cls=rep(1,nrow(Points)),xlab='X',ylab='Y',main="DataBots"){
  X=Points[,1]
  Y=Points[,2]
    #ColorSymbSequence <- DefaultColorSymbSequence()
    PlotSymbol <- 20#ColorSymbSequence[1]
  
    #DefaultColorSeq <- DefaultColorSequence()
    # if (missing(ColorSequence))
    ColorSequence <- DatabionicSwarm::DefaultColorSequence
      
      NormalizeCls_hlp <- function(Cls) {
        #E<-NormalizeCls(Cls);
        #NormalizedCls    <- E$normalizedCls      #    Cls consistently recoded to positive consecutive integers
        #NormalizedClasses<- E$normalizedClasses  #    the different class numbers in NormalizedCls
        #UniqueCls        <- E$uniqueClasses      #    the different class numbers in Cls such that 
        #AnzClasses       <- E$numberOfClasses    #    the number of different classes
        # 
        # Values in Cls are consistently recoded to positive consecutive integers
        # INPUT
        # Cls                  vector of class identifiers can be integers or
        #                      NaN's, need not be consecutive nor positive
        # OUTPUT list of 
        # normalizedCls           Cls consistently recoded to positive consecutive integers
        # normalizedClasses        the different class numbers in NormalizedCls
        # uniqueClasses            the different class numbers in Cls such that 
        #                           NormalizedCls(i) <-> UniqueCls(i)
        # numberOfClasses           the number of different classes
        
        # ALU 2014
        # 1.Editor:MT 2016
        
        uniqueClasses <- sort(na.last=T,unique(Cls))
        numberOfClasses <- length(uniqueClasses)
        unique2Cls <- NULL #  initializing the vector
        
        for (i in 1:length(Cls) ) { # calculating the indexes of elements of Cls in uniqueClasses
          unique2Cls <- c( unique2Cls, which(uniqueClasses == Cls[i]))
        } 
        
        if (numberOfClasses > 0) {
          normalizedClasses <- c(1: numberOfClasses)
          normalizedCls <- normalizedClasses[unique2Cls]
        }
        else {
          normalizedClasses <- Cls
        }
        
        return(list(normalizedCls = normalizedCls, normalizedClasses = normalizedClasses, uniqueClasses = uniqueClasses, numberOfClasses = numberOfClasses))
      }
    E <- NormalizeCls_hlp(Cls)
    NormalizedCls <- E$normalizedCls
    UniqueCls <- E$uniqueClasses
    AnzClasses <- E$numberOfClasses
  
  
    AnzColors = length(ColorSequence)
    ColorNR = c(1:AnzColors)
    if (AnzClasses > AnzColors) {
      ColorNR = (c(0:AnzClasses)%%AnzColors) + 1
    }
    MinX = min(X,na.rm=T)
    MaxX = max(X,na.rm=T)
    MinY = min(Y,na.rm=T)
    MaxY = max(Y,na.rm=T)
    xlim=c((MinX-abs(0.1*MinX)),(MaxX+abs(0.1*MaxX)))
    ylim=c((MinY-abs(0.1*MinY)),(MaxY+abs(0.1*MaxY)))
  
  #Initialisierungsplot
  plot.new()
  #par(usr=c(MinX,MaxX,MinY,MaxY))
  par(usr=c(xlim,ylim))
  par(xaxs='i')
  par(yaxs='i')
  # Plot der Punkte
    for (i in 1:AnzClasses) {
      C = UniqueCls[i]
      Ind = which(C == NormalizedCls)
      points(X[Ind], Y[Ind], pch = PlotSymbol,
             col = ColorSequence[ColorNR[i]], new = TRUE)
    }
  
  axis(1,xlim=xlim,col="black",las=1) #x-Achse
  axis(2,ylim=ylim,col="black",las=1) #y-Achse
  title(xlab=xlab,ylab=ylab,main=main)
  box()
}