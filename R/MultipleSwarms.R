MultipleSwarms=function(InputDistances,NumberOfSwarms=12,LogicalProcessors=4,PlotIt=F,Cls=NULL,Silent=T,Debug=FALSE){
# bmus=pswarmCpp(DataOrDists,PlotIt=T,Cls)
# Laesst multiple Pswarm Schwarmalgorithmen parallel ueber einen Datensatz laufen
# anhand spieltheorie wird pro eppoche das beste nash equilibirum ausgewaehlt, siehe diss
# polar oriatated swarm
#
# INPUT
# Dists[1:n,1:d]            distance matrix, the  matrix has to be symmetric, 
#                           if its not symetric data matrix is assumed and Dists=DistanceMatrix(data); 
#                           data matrix is array of data: n cases in rows, d variables in columns, matrix is not symmetric
#
# OPTIONAL
# NumberOfSwarms            Number of pararell Pswarms
# LogicalProcessors         Numeber of LogicalProcessors. How many games should be simoultansly computated
#
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
  
  LC=c(NULL,NULL)

  requireNamespace('parallel')
  
  if (missing(InputDistances)) {
    stop('Distances are Missing.')
  } else{
    if (is.list(InputDistances)) {
      stop('InputDistances is a list! It has to be a matrix.')
    }
    if (!is.matrix(InputDistances)) {
      {
        warning('InputDistances is not a matrix. Trying to circumvent...')
        if (is.data.frame(InputDistances))
          InputDistances = data.matrix(InputDistances)
        else
          InputDistances = as.matrix(InputDistances)
      }
    }
  }
  if (!Silent)
    print('Operator: Setting options')
  
  QuadOrHexa = F #hesagonales Gitter => dichteste Kugelpackung
  # if(PlotIt)
  #   tryCatch({requireNamespace("plotrix")},error=function(ex) {})#zum plotten
  #
  if (isSymmetric(InputDistances)) {
    DataDists = InputDistances
    AnzVar = ncol(InputDistances)
    AnzData = nrow(InputDistances)
  } else{
    #!isSymmetric
    warning('Distances are not in a symmetric matrix, Datamatrix is assumed and dist() ist called')
    
    DataDists = as.matrix(dist(InputDistances, method = "euclidean", diag =
                                 TRUE))
    AnzVar = ncol(DataDists)
    AnzData = nrow(DataDists)
  }# end if(isSymmetric(DataOrDists))
  
  if (is.null(LC[1]))
    LC = setGridSize(InputDistances)
  else{
    if (!is.vector(LC))
      stop('LC has to be a vector')
    if (is.list(LC))
      stop('LC has to be a vector not a list')
    if (!is.numeric(LC))
      stop('LC has to be numeric')
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
  
  reldiffp = function(x, y) {
    if (x + y == 0)
      return(0)
    return(signif((y - x) / (0.5 * (x + y)), 2) * 100)
  }
  
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
  DataBotsPos = dballind[init,]
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
    print('Operator: Starting algorithm')
  
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
  #library(parallel)
 
  if(NumberOfSwarms<LogicalProcessors)
    numWorkers = NumberOfSwarms #min=8, max=48 bei 4-8gb ram und 4 prozessoren
  else
    numWorkers = LogicalProcessors 
  Workers = parallel::makeCluster(numWorkers, type = "PSOCK")
  
  for (Radius in rvec) {
    if (!Silent)
      print(paste0('Operator: Current focus: ', Radius))
    
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
    Listen = parallel::parLapply(
      cl = Workers,
      1:NumberOfSwarms,
      FUN = function(i,
                     AllDataBotsPos,
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
      return(
        PswarmCurrentRadiusC2botsPositive(
          AllDataBotsPos,
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
          Debug
        )
      ),
      AllDataBotsPos,
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
      Debug
    )
    # Das Kriterium zur auswahl ist mir noch unklar

    #ResultsStress=sapply(1:NumberOfSwarms,function(i,xx) return(length(xx[[i]]$stressverlauf)),Listen)
    ResultsStress=sapply(1:NumberOfSwarms,function(i,xx) return(reldiffp(xx[[i]]$stressverlauf[1], tail(xx[[i]]$stressverlauf, 1))),Listen)

    ind = which.max(ResultsStress)
    #Listen=list(List1,List2,List3)
    List = Listen[[ind]]
    AllDataBotsPos = List$AllDataBotsPos
    stressverlauf = c(stressverlauf, List$stressverlauf)
    eppocheradiusreduziert = c(eppocheradiusreduziert, List$fokussiertlaufind)
    if (!Silent) {
      print(paste0(
        'Operator: ',
        tail(eppocheradiusreduziert, 1),
        '.iteration'
      ))
      print(
        paste0(
          'Operator: weak Nash equilibrium found. Paypoff maximized with ',
          reldiffp(List$stressverlauf[1], tail(List$stressverlauf, 1)),
          ' %'
        )
      )
    }
  } #end for rvec
  parallel::stopCluster(Workers)
  bmu = getCartesianCoordinates(AllDataBotsPos,
                                GridRadius = GridRadii,
                                GridAngle,
                                QuadOrHexa = QuadOrHexa)
  
  return(list(
    ProjectedPoints = bmu,
    LC = c(Lines, Columns),
    Control = list(
      stressverlauf = stressverlauf,
      eppocheradiusreduziert = eppocheradiusreduziert,
      LetzteEppocheStress = stress
    )
  ))
}
