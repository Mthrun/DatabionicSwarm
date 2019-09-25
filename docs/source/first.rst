**First Example: Automatic approach**
=====================================

Currently, this Documentation is under development. For now please see to the vignette in order to learn how the algorithm can be applied.

Here one example is presented using the automatic approach without any user interaction with shiny. Further automatic examples and a comparison to 26 common clustering algorithms is provided in http://www.deepbionics.org/Projects/ClusteringAlgorithms.html. If you want to verify your clustering result externally, you can use **Heatmap** or **SilhouettePlot** of the CRAN package ``DataVisualizations``.

First Module: Projection of high-dimensional Data
---------------------------------------------------


First generate a 2d projection, the DistanceMatrix has to be defined by the user.

.. code-block:: R

	library(DatabionicSwarm)

::
	
	## Package 'DatabionicSwarm' version 1.1.1.
	## Type 'citation('DatabionicSwarm')' for citing this R package in publications.

.. code-block:: R

	data('Hepta')
	InputDistances=as.matrix(dist(Hepta$Data))
	projection=Pswarm(InputDistances)

Second Module: Generalized Umatrix
-----------------------------------

Here the Generalized Umatrix is calculated using a simplified emergent self-organizing map algorithm. Then, the visualization of Generalized Umatrix is done by a 3D landscape called ``topographic map`` with hypsometric tints. Seven valleys are shown resulting in seven main clusters. The resulting visualization will be toroidal meaning that the left borders cyclically connects to the right border (and bottom to top). It means there are no **real** borders in this visualizations. Instead, the visualization is **continuous**. This can be visualized using the ``Tiled=TRUE`` option of ``plotTopographicMap``.

Note, that the **NoLevels** option is only set to load this vignette faster and should normally not be set manually. It describes the number contour lines placed relative to the hypsometric tints. All visualizations here are small and a low dpi is set in knitr in order to load the vignette faster.

.. code-block:: R

	library(DatabionicSwarm)
	library(GeneralizedUmatrix)
	visualization=GeneratePswarmVisualization(Data = Hepta$Data,projection$ProjectedPoints,projection$LC)
	GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)

::

	## Loading required namespace: matrixStats

.. code-block:: R

	rgl::rgl.close()#please ignore, indicates that this plot should not be saved in Rmarkdown
	
Third Module: Automatic Clustering
--------------------------------------

The number of clusters can be derived from ``dendrogram`` ``(PlotIt=TRUE)`` or the visualization. Therefore we choose the seven valleys as the number of clusters. The function DBSclustering has one parameter to be set. Typically, the default setting **StructureType = TRUE** works fine. However, for density-based structures sometimes ``StructureType = FALSE`` of the function **DBSclustering** yields better results. Please verify with the visualization or the Dendrogram. For the Dendrogram choose ``PlotIt=TRUE`` in the function **DBSclustering**. In the case of **BestmatchingUnits**, the parameter **LC** defines the size of the grid with Lines and columns where the position (0,0) lies in the left upper corner. In the case of **ProjectedPoints**, the point (0,0) lies in the left bottom corner. The transformation is normally done automatically. However, sometimes the user wishes to skip the visualization and use projected points directly. Then **LC** can be changed accordingly to LC[c(2,1)]. Seldom, there could be a rounding error leading to an error catch. In such a case try LC+1.

.. code-block:: R

	library(DatabionicSwarm)
	library(GeneralizedUmatrix)
	Cls=DBSclustering(k=7, Hepta$Data, visualization$Bestmatches, visualization$LC,PlotIt=FALSE)

::

	## Loading required namespace: parallelDist
	
::

	## 
	##      PLEASE NOTE:  The components "delsgs" and "summary" of the
	##  object returned by deldir() are now DATA FRAMES rather than
	##  matrices (as they were prior to release 0.0-18).
	##  See help("deldir").
	##  
	##      PLEASE NOTE: The process that deldir() uses for determining
	##  duplicated points has changed from that used in version
	##  0.0-9 of this package (and previously). See help("deldir").
	
.. code-block:: R

	GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls,NoLevels=10)
	
.. raw:: html
   :file: TopographicMap.html
   
