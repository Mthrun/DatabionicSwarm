setPolarGrid <-
function(Lines,Columns,QuadOrHexa=T,PlotIt=F,global=T){
  #  setPolarGrid (Lines,Columns,QuadOrHexa=F,PlotIt=F)
  #INPUT
  # Lines         Integer, hast to be able to be divided by 2
  # Columns       Integer, with Columns>Lines
  # OPTINAL
  # quadOrHex     BOOL, If False Hexagonal grid, default quad grid
  # PlotIt
  # global        Bool, Wie sollen die moeglichen Radien bestimmt werden
  #OUTPUT list V with
 # V$GridRadii[Lines,Columns]                     Radii Matrix of all possible Positions of DataBots in Grid
  # V$GridAngle[Lines,Columns]                     Angle Matrix of all possible Positions of DataBots in Grid
  # V$AllallowedDBPosR0[Lines+1,Lines+1]          Matrix of radii in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  # V$AllallowedDBPosPhi0[Lines+1,Lines+1]        Matrix of angle in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  # author: MT 12/2014
  #achtung unter umstaenden darf auch Lines==Columns nicht gelten, muss ich noch pruefen
  if(Lines%%2!=0){stop('Lines has to be even')}
  if(Lines>Columns){stop('Number of Columns has to be higher than number of Lines')}
  
  if(QuadOrHexa){
    #Generierung alles Moeglichen Position bezueglch Urpsrung in polaren Koordinaten
    #################
    if(global){
      
    xgrade2=seq(by=1,from=-Columns/2,to=Columns/2)
    
    unitsx2=rep(xgrade2,Columns+1)
    
    ygrade2=seq(by=1,from=-Columns/2,to=Columns/2)
    
    unitsy2=sort(na.last=T,rep(ygrade2,Columns+1))
    
    #Bestimmung aller moeglichen Positionen innerhalb der Karte
    posR2=sqrt((unitsx2)^2+(unitsy2)^2)
    posPhi2=atan2(y=unitsy2,x=unitsx2)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
    posPhi2[which(posPhi2<0)]<-posPhi2[which(posPhi2<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
    
    PositionR0=matrix(posR2,nrow=Columns+1,ncol=Columns+1,byrow=TRUE)
    PositionPhi0=matrix(posPhi2,nrow=Columns+1,ncol=Columns+1,byrow=TRUE)
    }else{
        xgrade2=seq(by=1,from=-Lines/2,to=Lines/2)
        
        unitsx2=rep(xgrade2,Lines+1)
        
        ygrade2=seq(by=1,from=-Lines/2,to=Lines/2)
        
        unitsy2=sort(na.last=T,rep(ygrade2,Lines+1))
        
        #Bestimmung aller moeglichen Positionen innerhalb der Karte
        posR2=sqrt((unitsx2)^2+(unitsy2)^2)
        
        posPhi2=atan2(y=unitsy2,x=unitsx2)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
        posPhi2[which(posPhi2<0)]<-posPhi2[which(posPhi2<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
   
        PositionR0=matrix(posR2,nrow=Lines+1,ncol=Lines+1,byrow=TRUE)
        PositionPhi0=matrix(posPhi2,nrow=Lines+1,ncol=Lines+1,byrow=TRUE)
    }

    ##################
    #Berechnen das eigentliche Gitter im polaren Koordinaten
    xgrade=seq(by=1,from=0,to=Columns-1)
    
    unitsx=rep(xgrade,Lines)
    
    ygrade=seq(by=1,from=0,to=Lines-1)
    
    unitsy=sort(na.last=T,rep(ygrade,Columns))
    ##################
    #Bestimmung aller moeglichen Positionen innerhalb der Karte
    posR=sqrt((unitsx)^2+(unitsy)^2)
    
    posPhi=atan2(y=unitsy,x=unitsx)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
    posPhi[which(posPhi<0)]<-posPhi[which(posPhi<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
    ##Aus den Positionen werden drei Matrizen generiert 
    MatrixR=matrix(posR,nrow=Lines,ncol=Columns,byrow=TRUE)
    MatrixW=matrix(posPhi,nrow=Lines,ncol=Columns,byrow=TRUE)
    binary=matrix(0,nrow=Lines,ncol=Columns)
    
    
  }else{ #Or in Hexadiagonal
    #Generierung alles Moeglichen Position bezueglch Urpsrung in polaren Koordinaten
    #################
    if(global){
      if(Columns %%4 ==0){
      xgrade2=seq(by=1,from=-Columns/2+0.5,to=Columns/2+1)
      xungrade2=seq(by=1,from=-Columns/2,to=Columns/2+0.5)
      unitsx2=c(rep(c(xungrade2,xgrade2),Columns/2),xungrade2)
      
      ygrade2=seq(by=1,from=-Columns/2,to=Columns/2+0.5)

      unitsy2=sort(na.last=T,rep(ygrade2,Columns+1))
      }else{
        stop('Hexadiagonal is not Implemented for number of Column length, which cannot be divided by four')
#         xgrade2=seq(by=1,from=-Columns/2,to=Columns/2)
#         xungrade2=seq(by=1,from=-Columns/2-0.5,to=Columns/2+0.5)
#         unitsx2=c(rep(c(xungrade2,xgrade2),Columns/2-0.5),xungrade2)
#       
#         ygrade2=seq(by=1,from=-Columns/2+0.5,to=Columns/2+1)
#       
#         unitsy2=sort(na.last=T,rep(ygrade2,Columns))
      }
      #Bestimmung aller moeglichen Positionen innerhalb der Karte
      posR2=sqrt((unitsx2)^2+(unitsy2)^2)
      posPhi2=atan2(y=unitsy2,x=unitsx2)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
      posPhi2[which(posPhi2<0)]<-posPhi2[which(posPhi2<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen

      PositionR0=matrix(posR2,nrow=Columns+1,ncol=Columns+1,byrow=TRUE)
      PositionPhi0=matrix(posPhi2,nrow=Columns+1,ncol=Columns+1,byrow=TRUE)
    }else{
      stop('not implemented yet')
#       xgrade2=seq(by=1,from=-Lines/2,to=Lines/2)
#       
#       unitsx2=rep(xgrade2,Lines+1)
#       
#       ygrade2=seq(by=1,from=-Lines/2,to=Lines/2)
#       
#       unitsy2=sort(na.last=T,rep(ygrade2,Lines+1))
#       
#       #Bestimmung aller moeglichen Positionen innerhalb der Karte
#       posR2=sqrt((unitsx2)^2+(unitsy2)^2)
#       
#       posPhi2=atan2(y=unitsy2,x=unitsx2)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
#       posPhi2[which(posPhi2<0)]<-posPhi2[which(posPhi2<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
#       
#       PositionR0=matrix(posR2,nrow=Lines+1,ncol=Lines+1,byrow=TRUE)
#       PositionPhi0=matrix(posPhi2,nrow=Lines+1,ncol=Lines+1,byrow=TRUE)
    }
    
    ##################
    #Berechnen das eigentliche Gitter im polaren Koordinaten
    xgrade=seq(by=1,from=0.5,to=Columns)
    xungrade=seq(by=1,from=0,to=Columns-0.5)
    
    unitsx=rep(c(xungrade,xgrade),Lines/2)
    #unitsx=rep(xgrade,Lines)
    ygrade=seq(by=1,from=0,to=Lines-0.5)
    
    unitsy=sort(na.last=T,rep(ygrade,Columns))
    
    
    posR=sqrt((unitsx)^2+(unitsy)^2)
    
    posPhi=atan2(y=unitsy,x=unitsx)*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
    posPhi[which(posPhi<0)]<-posPhi[which(posPhi<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
        
    MatrixR=matrix(posR,nrow=Lines,ncol=Columns,byrow=TRUE)
    MatrixW=matrix(posPhi,nrow=Lines,ncol=Columns,byrow=TRUE)
    binary=matrix(0,nrow=Lines,ncol=Columns)
  }
  
  if(PlotIt){   
   # tryCatch({requireNamespace('plotrix')},error=function(ex) {})#zum plotten     
    requireNamespace('plotrix')
    Rdiag=sqrt((Lines)^2+(Columns)^2)
    gLines=Lines-0.25#nur fuer Rahmen der Karte=Rechteck
    gColumns=Columns-0.25#nur fuer Rahmen der Karte
    
    points=cbind(x=c(gColumns/2,gColumns,gColumns,gColumns,gColumns/2,0,0,0,gColumns/2),
                 y=c(0,0,gLines/2,gLines,gLines,gLines,gLines/2,0,0))      
    Rvec=sqrt((points[,'x'])^2+(points[,'y'])^2) 
    phivec=atan2(y=points[,'y'],x=points[,'x'])*180/(pi) #hier koennen auch negative Gradzahlen rauskommen, falls zwischen 180 und 360
    phivec[which(phivec<0)]<-phivec[which(phivec<0)]+360 #hier die negativen Gradzahlen in entsprechende zwischen 180 und 360 umrechnen
    
    oldpar<-plotrix::polar.plot(Rvec,phivec,main="Gitter",
                       radial.lim=c(0,Rdiag),point.symbols=4,start=0, rp.type='p',show.grid=TRUE) 
    
    oldpar<-plotrix::polar.plot(posR2,posPhi2,main="Gitter",
                       radial.lim=c(0,Rdiag),point.symbols=1,start=0, rp.type='s',show.grid=TRUE,add=T) 
    
    oldpar<-plotrix::polar.plot(posR,posPhi,main="Gitter",
                       radial.lim=c(0,Rdiag),point.symbols=1,start=0, rp.type='s',show.grid=TRUE,add=TRUE)   
  }
  
  return(list(GridRadii=MatrixR,GridAngle=MatrixW,AllallowedDBPosR0=PositionR0,AllallowedDBPosPhi0=PositionPhi0,Binary=binary))
  
}
