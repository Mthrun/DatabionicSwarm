
**Second Example: Interactive Approach**
========================================

In this example, we show how to improve an automatic clustering accordingly to the topographic map using a shiny interface. The examples are not run in Rmarkdown, because CRAN wants to check the results regularly and this example is above their time limit. You may see http://www.deepbionics.org/Rpackages.html for the complete examples.

First Module: Projection of High-dimensional Data
--------------------------------------------------

First, we generate a 2d projection with instant visualization of annealing steps ``(PlotIt=TRUE)``. This shows the non-linear process of concentrating on global structures first and later on local structures. Such an approach enables to entangle non-linear high-dimensional structures. If the user does not define the DistanceMatrix, it is automatically set to Euclidean because the data itself can be the input for ``Pswarm``.

.. code-block:: R

	rgl::rgl.close()#please ignore, closes last plot of rgl to not show it multiple times
	library(DatabionicSwarm)
	data('Lsun3D')
	projection=Pswarm(Lsun3D$Data,Cls=Lsun3D$Cls,PlotIt=T,Silent=T)

Second Module: Generalized Umatrix
----------------------------------

If **Non-Euclidean** Distances are used, Please Use SammonsMapping from the ``ProjectionBasedClustering`` package with the correct OutputDimension to generate a new data matrix from the distances (see SheppardDiagram of ``DataVisualization`` Package or KruskalStress). Here the Generalized Umatrix is calculated using a simplified emergent self-organizing map algorithm. Then the topographic map is visualized based on the information of the Generalized Umatrix.

.. code-block:: R

	library(DatabionicSwarm)
	library(GeneralizedUmatrix)
	visualization=GeneratePswarmVisualization(Data = Lsun3D$Data,projection$ProjectedPoints,projection$LC)

	GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)
	rgl::rgl.close()#please ignore, indicate that this plot should not be saved in Rmarkdown

Third Module: Interactive Clustering
------------------------------------

The number of clusters can be derived from dendrogram ``(PlotIt=TRUE)`` or the visualization. In this example, outliers should be marked manually in the visualization after the process of automatic clustering. Therefore we choose the three central valleys as the number of clusters. Often, it helps to generate first the shape of an island out of the continuous topographic map because then you already have the most prominent mountains marked as the borders of the visualizations. Then you can improve the clustering by redefining valleys interactively or marking outliers lying in vulcanos. It is strongly suggested to verify such a clustering externally, e.g. **Heatmap** or some unsupervised index.

.. code-block:: R

	library(DatabionicSwarm)
	library(GeneralizedUmatrix)
	Cls=DBSclustering(k=3, Lsun3D$Data, visualization$Bestmatches, visualization$LC,PlotIt=FALSE)
	GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)

Generating the Shape of an Island out of the Topograpahic Map
-------------------------------------------------------------

To generate the 3D landscape in the shape of an island from the toroidal topographic map visualization you may cut your island interactively around high mountain ranges. Currently, I am unable to show the output in R markdown :-( If you know how to resolve the Rmarkdown issue, please mail me: info@deepbionics.org

.. code-block:: R

	library(DatabionicSwarm)
	library(ProjectionBasedClustering)
	library(GeneralizedUmatrix)

	Imx = ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(visualization$Umatrix,visualization$Bestmatches,Cls)

	GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches, Cls=Cls,Imx = Imx)
	
Manually Improving the Clustering Using the Topograpahic Map
------------------------------------------------------------

In this example, the four outliers can be marked manually with mouse clicks using the shiny interface. Currently, I am unable to show the output in R markdown :-( Please try it out yourself:

.. code-block:: R

	library(ProjectionBasedClustering)
	Cls2=ProjectionBasedClustering::interactiveClustering(visualization$Umatrix, visualization$Bestmatches, Cls)
	
