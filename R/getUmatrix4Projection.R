getUmatrix4Projection=function(Data,ProjectedPoints,PlotIt=TRUE,Cls=NULL,toroid=T,Tiled=F,ComputeInR=F){
#V=getUmatrix4Projection(ProjectedPoints,Data)
#V=getUmatrix4Projection(ProjectedPoints,Data,TRUE,Cls,TRUE)
# Generalisierte U-Matrix fuer Projektionsverfahren
# Erklaerung der Funktion in
#D:\Subversion\lehre\Vorlesungen\KnowledgeDiscovery\01Text\04ProjektionenUndVisualisierung\81DatenbionischeProjektionen\05Schwaerme\DataBionicSwarm.pptx
#Folie 45-51
# INPUT
# Data[1:n,1:d]                          array of data: n cases in rows, d variables in columns
# ProjectedPoints[1:n,OutputDimension]   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
                                         # see Projections/R/..
# OPTIONAL
# PlotIt                                 bool, defaut=FALSE, if =TRUE: U-Marix of every current Position of Databots will be shown
# toroid
## ComputeInR                                  =T: Rcode, =F Cpp Code

# Output
# Umatrix[Lines,Columns                     umatrix (see ReadUMX() dbt.DataIO)
# EsomNeurons[Lines,Columns,weights]       3-dimensional numeric array (wide format), not wts (long format)
# Bmu[1:n,OutputDimension]                 GridConverted Projected Points information converted by convertProjectionProjectedPoints() to predefined Grid by Lines and Columns
# gplotres                                Ausgabe von ggplot
# unbesetztePositionen                    Umatrix[unbesetztePositionen] =NA
# author: MT 06/2015
#1.Editor; MT 12/2015

warning('depricated, calls GeneralizedUmatrix() function')
requireNamespace('GeneralizedUmatrix')
return(GeneralizedUmatrix::GeneralizedUmatrix(Data=Data,ProjectedPoints=ProjectedPoints,PlotIt=PlotIt,Cls=Cls,Toroid=toroid,Tiled=Tiled,ComputeInR=ComputeInR))
}
