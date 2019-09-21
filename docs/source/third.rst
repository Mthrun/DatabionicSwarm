
**Quality Measures of Projection and Clustering**
=================================================

Delaunay Classification Error
-----------------------------

Using insights of graph theory, the Delaunay classification error calculates for each projected point the direct neighborhood based on the **Delaunay graph**. Every direct Connection is weighted with the high-dimensional distance between the two corresponding data points and sorted per neighborhood by these weights. In the next step all sorted projected points points of the direct neighborhood of each projected points aquire new weights according to the harmonic function. Then, the prior classification is used to check which points do not belong to these direct neighborhoods of projected points. The weights of these points are summed up. A lower value indicates a better two-dimensional projection of the high-dimensional Input space. A higher value indicates a worse two-dimensional projection of the high-dimensional Input space.

.. code-block:: R

	DelaunayClassificationError(Lsun3D$Data,projection$ProjectedPoints,Lsun3D$Cls)
	
You can also compare various projections method to a common baseline together:

.. code-block:: R

	DCEpswarm=DelaunayClassificationError(Lsun3D$Data,projection$ProjectedPoints,Lsun3D$Cls)$DCE
	baselineproj=ProjectionBasedClustering::NeRV(Lsun3D$Data)
	DCEpswarm=DelaunayClassificationError(Lsun3D$Data,baselineproj,Lsun3D$Cls)$DCE
	RelativeDifference(DCEpswarm,baselineproj)

This has the advantage of an clear range of [−2,2]. Further Details can be read in the conference presentation attached to [Thrun/Ultsch,2018] on ResearchGate.

Clustering Accuracy
-------------------

The accuracy is defined as follows:``Accuracy=No.oftruepositives/No.ofcases`` The number of **true positives** is the number of labeled data points for which the label defined by a prior classification is identical to the label defined after the clustering process. The best of all permutation of labels of the clustering algorithm regarding the accuracy is chosen, because the labels are arbitrarily defined by any algorithm. See details in conference presentation attached to [Thrun et al.,2018] on ResearchGate.

.. code-block:: R
	
	ClusteringAccuracy(Lsun3D$Cls,Cls)

