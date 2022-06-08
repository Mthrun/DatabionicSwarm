Pswarm = pswarmCpp = function(DataOrDistance,PlotIt=F,Cls=NULL,Silent=T,Debug=FALSE,LC=c(NULL,NULL),method='euclidean',...){
# bmus=pswarmCpp(DataOrDists,PlotIt=T,Cls)
# Laesst den Pswarm Schwarmalgorithmus ueber einen Datensatz laufen
# polar oriatated swarm
#
# INPUT
# DataOrDistance[1:n,1:d]   distance matrix, the  matrix has to be symmetric, 
#                           if its not symetric data matrix is assumed and Dists=DistanceMatrix(data); 
#                           data matrix is array of data: n cases in rows, d variables in columns, matrix is not symmetric
#
# OPTIONAL
# PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made. At the end of pswarm Plots of sum(stress) will be mad
# Cls                       vector, Klassifikation of Data if available, ClassPlots will be colorized
#
# Silent                    =FALSE: No print Output, =TRUE some print outs
# Debug                     =TRUE: Debigging Modus, Slow with alot of prints
# OUTPUT Liste V
# V$BestMatchingUnits[1:n,1:3]        n by 2 matrix containing  X and Y coordinates of the n BestMatches for each databot including unique key
#                           BestMatches need to be unique. Transformation from polar (R,phi) to kartesisch (x,y) is done automatically
# V$Grid                    vector aus c(Lines,Columns)
# V$Control                 Alle m?glichen Parameter, um nachzuvollziehen, wenn was schiefgeht
#  
# Autor: MT 01/2015
# Nota: im debugging modus sollten relativ differenzen des payoffs angegeben werden statt festen werten
  
#############################
## Not required anymore  
  # LC                        Vector of Lines and Colums, specific:
# Lines   Default=50, x-value determining the size of the map, i.e. how many open places for DataBots will be available  on the 2-dimensional grid
#  				BEWARE: must be divisible by two
#
# Columns Default=80, y-value  determining the size of the map, i.e. how many open places for DataBots will be available on the 2-dimensional grid
#					Columns>Lines.
#						                # BEWARE: If you choose hexagonal, Column number must be divisible by four 
#

  DataOrDistance=checkInputDistancesOrData(DataOrDistance,funname='Pswarm')
  
  # if (missing(DataOrDistance)) {
  #   stop('Pswarm: Distances are Missing.')
  # } else{
  #   if (is.list(DataOrDistance)) {
  #     stop('Pswarm: DataOrDistance is a list! It has to be a matrix.')
  #   }
  #   if (!is.matrix(DataOrDistance)) {
  #     {
  #       warning('Pswarm: DataOrDistance is not a matrix. Trying to circumvent...')
  #       if (is.data.frame(DataOrDistance))
  #         DataOrDistance = data.matrix(DataOrDistance)
  #       else
  #         DataOrDistance = as.matrix(DataOrDistance)
  #     }
  #   }
  # }
  if (!Silent)
    message('Operator: Setting options')
  
  QuadOrHexa = F #hesagonales Gitter => dichteste Kugelpackung
  # if(PlotIt)
  #   tryCatch({requireNamespace("plotrix")},error=function(ex) {})#zum plotten
  #
  if (isSymmetric(unname(DataOrDistance))) {
    DataDists = DataOrDistance
    AnzVar = ncol(DataOrDistance)
    AnzData = nrow(DataOrDistance)
  } else{
    #!isSymmetric
	if(!Silent)
		message('Distances are not in a symmetric matrix, Datamatrix is assumed and parallelDist::parDist() ist called')
    
	if (!requireNamespace('parallelDist',quietly = TRUE)) {
		message(
		  'Subordinate parallelDist package is missing. dist() function of stats is used.
				Please install the package which is defined in "Suggests", if other distances than available in dist() or faster distance computation is necessary'
		)
		DataDists = as.matrix(dist(DataOrDistance, method = method))
	}else{
		DataDists = as.matrix(parallelDist::parDist(DataOrDistance, method = method,...))
		AnzVar = ncol(DataDists)
		AnzData = nrow(DataDists)
	}
  }# end if(isSymmetric(DataOrDists))

  if (is.null(LC[1]))
    LC = setGridSize(DataDists)
  else{
	if(!is.vector(LC)) stop('LC has to be a vector')
    if(is.list(LC)) stop('LC has to be a vector not a list')
	if(!is.numeric(LC)) stop('LC has to be numeric')
  }
  DBAnzahl = AnzData
  if (is.null(Cls))
    Cls = rep(1, DBAnzahl)
  
  Lines = LC[1]
  Columns = LC[2]
  if (Lines %% 2 != 0) {
    stop('Lines has to be even')
  }
  if (Columns %% 4 != 0) {
    stop('Columns has to be dividable by four')
  }
  #achtung unter umstaenden darf auch Lines==Columns nicht gelten, muss ich noch pruefen
  # in setPolarGrid
  if (Lines > Columns) {
    stop('Number of Columns has to be higher than number of Lines')
  }
  if (Lines * Columns < DBAnzahl) {
    stop('Map is too small for the Number of DataBots, please chose higher Lines and Columns')
  }
  Rmax = Lines / 2
  
  #reldiffp=function(x,y){
  #  if(x+y==0) return(0)
  #  return(signif((y-x)/(0.5*(x+y)),2)*100)
  #}
  
  # Initialisierung von ben?tigten InputVariablen
  ################################################################################################
  # Mein Algorithmus
  ################################################################################################
  ListeDerPositionsSchablonen = setPolarGrid(
    Lines = Lines,
    Columns = Columns,
    QuadOrHexa = QuadOrHexa,
    PlotIt = F
  )                     #Binary matrix of all free DataBotsPositions(=0), Positions of DataBots will be defined by 1
  GridRadii = ListeDerPositionsSchablonen$GridRadii                     #Radii Matrix of all possible Positions of DataBots in Grid
  GridAngle = ListeDerPositionsSchablonen$GridAngle                   #Angle Matrix of all possible Positions of DataBots in Grid
  AllallowedDBPosR0 = ListeDerPositionsSchablonen$AllallowedDBPosR0         #Radius-Matrix in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  AllallowedDBPosPhi0 = ListeDerPositionsSchablonen$AllallowedDBPosPhi0#Angle-Matrix in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  ## Geht schief wenn LC zu klein
  init = sample(x = Lines * Columns,
                replace = F,
                size = DBAnzahl)
  binary = matrix(0, nrow = Lines, ncol = Columns)
  dballind = which(binary == 0, arr.ind = T)
  DataBotsPos = dballind[init, ]
  #richtig cooler Trick
  AllDataBotsPos = DataBotsPos[, 1] + 1i * DataBotsPos[, 2]
  Nullpunkt = which(AllallowedDBPosR0 == 0, arr.ind = T)
  
  if (QuadOrHexa) {
    #Bei quadratischen Gitter ist minimaler Radius groesser 3 noetig
    Rmin = setRmin(AllallowedDBPosR0, Lines, Columns, DBAnzahl, p = 0.05)
  } else{
    # !QuadOrHexa
    #Rmin=1
    Rmin = setRmin(AllallowedDBPosR0, Lines, Columns, DBAnzahl, p = 0.05)
  }# end if QuadOrHexa
  #pp=seq(1,100,4)/100 #FuerAnzahl gleichzeitig springender DataBots
  m = (0.5 - 0.05) / (Rmax - Rmin)
  b = 0.5 - m * Rmax
  pp = m * 1:Rmax + b
  #pp=(1:Rmax)/Rmax-Rmin/Rmax+1/DBAnzahl
  #pp=pmax(0.05,pp)
  # pp=pmin(0.5,pp)
  # pp=pmax((1:Rmax)/Rmax-Rmin/Rmax,0.15)
  jumpthreshold = 0
  stressverlauf = c()
  if (!Silent)
    message('Operator: Preparing.')
  
  fokussiertlaufind = 1
  rvec = seq(from = Rmax, by = -1, to = Rmin)
  stress = Inf
  alpha = Rmax * 0.01
  
  RadiusPositionsschablone = ListeDerPositionsSchablonen$AllallowedDBPosR0
  IndPossibleDBPosR = findPossiblePositionsCsingle(RadiusPositionsschablone, Rmax, alpha, Lines)
  
  eppocheradiusreduziert = c()
  OutputDistance = rDistanceToroidCsingle(
    Re(AllDataBotsPos),
    Im(AllDataBotsPos),
    ListeDerPositionsSchablonen$AllallowedDBPosR0,
    Lines,
    Columns,
    as.vector(Nullpunkt)
  )
  Nachbahrschaftsfunktion = 1 - OutputDistance ^ 2 / (pi * Rmax ^ 2)
  Nachbahrschaftsfunktion[Nachbahrschaftsfunktion < 0] = 0
  N = sum(Nachbahrschaftsfunktion)
  StressConstAditiv = sum(Nachbahrschaftsfunktion * DataDists) / N

  dummy=0;
  numberOfSteps=length(rvec)
	if(Silent){
		progress = txtProgressBar(min = dummy, max = numberOfSteps+1, style = 3)
	}else{
	  ProzentualeZeitfolge=round(sort((rvec-Rmin)/(Rmax-Rmin),decreasing=F)*100,0)
	  ProzentualeZeitfolge[numberOfSteps]=99
	}
  if (!Silent)
    message('Operator: Starting algorithm')
  #----------------------------------------------------------------------------#
  for (Radius in rvec) {
		dummy=dummy+1
    if (!Silent){  
      message(paste0('Operator: ', ProzentualeZeitfolge[dummy],'% calculated.'))
      #message(paste0('Operator: Current focus: ', Radius))
			#setTxtProgressBar
    }else{
			 progressm = setTxtProgressBar(progress, dummy)
		}
      
    
    Eppoche = 1
    Jumping = TRUE
    limit = ceiling(1 / pp[Radius]) # Ab welcher Eppoche wird Abbruchbedingung geprueft
    steigungsverlaufind = 20#limit #wieviele zurueckliegende eppochen werden maximal geprueft
    if (PlotIt) {
      bmu = getCartesianCoordinates(AllDataBotsPos,
                                    GridRadius = GridRadii,
                                    GridAngle,
                                    QuadOrHexa = QuadOrHexa)
      string = paste0('Radius ', Radius, ', Eppoche ', Eppoche)
      plotSwarm(bmu, Cls, main = string)
    } #end if
    #Zeitfaktor: Je naeher DataBots springen, desto schneller riechen sie erneut entspricht, weniger DBs springen pro Eppoche
    nBots = round(pp[Radius] * DBAnzahl)
    #nBots=round(0.05*DBAnzahl)#s. Fast and reliable ESOM learning
    List = PswarmCurrentRadiusC2botsPositive(AllDataBotsPos,
                                             Radius,
                                             DataDists,
                                             IndPossibleDBPosR,
                                             RadiusPositionsschablone,
                                             pp,
                                             Nullpunkt,
                                             Lines,
                                             Columns,
                                             nBots,
                                             limit,
                                             steigungsverlaufind,
                                             StressConstAditiv,
                                             Debug)
    AllDataBotsPos = List$AllDataBotsPos
    stressverlauf = c(stressverlauf, List$stressverlauf)
    eppocheradiusreduziert = c(eppocheradiusreduziert, List$fokussiertlaufind)
    if (!Silent) {
      message(paste0('Operator: ', tail(eppocheradiusreduziert, 1), '.iteration'))
      message(
        paste0(
          'Operator: weak Nash equilibrium found. Paypoff maximized with ',
          signif(RelativeDifference(List$stressverlauf[1],tail(List$stressverlauf, 1)),2),
          ' %'
        )
      )
    }
  } #end for rvec
  #----------------------------------------------------------------------------#
  
  
  if (!Silent){  
    message(paste0('Operator: 100 % calculated.'))
  }else{
    progressm = setTxtProgressBar(progress, dummy+1)
    close(progress) 
  }

  bmu = getCartesianCoordinates(AllDataBotsPos,
                                GridRadius = GridRadii,
                                GridAngle,
                                QuadOrHexa = QuadOrHexa)
  #Possible Minor Rounding error
  if(max(bmu[,1])>Lines)
    Lines=Lines+1
  if(max(bmu[,2])>Columns)
    Columns=Columns+1
  
  return(list(
    ProjectedPoints = bmu,
    LC              = c(Lines, Columns),
    Control         = list(stressverlauf          = stressverlauf,
                           eppocheradiusreduziert = eppocheradiusreduziert,
                           LetzteEppocheStress    = stress
    )
  ))
}
