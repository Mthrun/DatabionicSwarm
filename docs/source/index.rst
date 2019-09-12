.. DatabionicSwarm documentation master file, created by
   sphinx-quickstart on Wed Sep 11 20:29:54 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*****************************************
Short Intro to the Databionic Swarm (DBS)
*****************************************

*Michael C. Thrun*

*27 Jan 2019*

***************
1 Introduction
***************

DBS is a flexible and robust clustering framework that consists of three independent modules. The first module is the parameter-free projection method Pswarm, which exploits the concepts of self-organization and emergence, game theory, swarm intelligence and symmetry considerations. The second module is a parameter-free high-dimensional data visualization technique, which generates projected points on a topographic map with hypsometric colors, called the generalized U-matrix. The third module is a clustering method with no sensitive parameters. The clustering can be verified by the visualization and vice versa. The term DBS refers to the method as a whole. For further details, see Databionic swarm in [Thrun, 2018], chapter 8.

***********************************
2 First Example: Automatic approach
***********************************

Here one example is presented using the automatic approach without any user interaction with shiny. Further automatic examples and a comparison to 26 common clustering algorithms is provided in http://www.deepbionics.org/Projects/ClusteringAlgorithms.html. If you want to verify your clustering result externally, you can use Heatmap or SilhouettePlot of the CRAN package DataVisualizations.

2.1 First Module: Projection of high-dimensional Data
#####################################################

First generate a 2d projection, the DistanceMatrix has to be defined by the user.
::

	library(DatabionicSwarm)

::
	
	## Package 'DatabionicSwarm' version 1.1.1.
	
	## Type 'citation('DatabionicSwarm')' for citing this R package in publications.

::

	data('Hepta')
	
	InputDistances=as.matrix(dist(Hepta$Data))
	
	projection=Pswarm(InputDistances)
