setRmin <-
function(AllallowedDBPosR0,Lines,Columns,DataBots,p=0.05){
# Rmin=setRmin(Lines,Columns,AllallowedDBPosR0,p=0.05)
# Bestimmt den minimalen Radius
#
# INPUT
# Lines                                   Default=50, x-value determining the size of the map, i.e. how many open places for DataBots will be available  on the 2-dimensional grid
#  				                                BEWARE: has to be able to be divided by 2
# Columns                                 Default=80, y-value  determining the size of the map, i.e. how many open places for DataBots will be available on the 2-dimensional grid
#					                                Columns>Lines 
# AllallowedDBPosR0[Lines+1,Lines+1]      Matrix of radii in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump# p                                   Prozentangabe als Zahl kleiner 1
# p                                       percent of gitterpositions, which should be considered
#
# OUTPUT
# Rmin                                     Minimum Radius         
# Autor: MT 01/2015
  
Positions=Lines*Columns
PosPerDataBot=Positions/DataBots
AnzahlinNaehe=p*DataBots
Rmax=Lines/2
Rmin=1
rvec=seq(from=Rmin,by=1,to=Rmax)
#Zaehle Anzahl der Plaetze innerhalb von Radius
verlauf=c()
k=1;
for(r in rvec){
  verlauf[k]=length(which(AllallowedDBPosR0<=r))
  k=k+1
}
#Es sollen mindestens p Prozent der Bots (AnzahlinNaehe) in der Naehe eines anderen Bots beim
# letzten Radius vorhanden sein, wenn die Bots auf dem gitter gleichverteilt waeren
ind=which.min(abs(verlauf-PosPerDataBot*AnzahlinNaehe))
rvec[ind]
}
