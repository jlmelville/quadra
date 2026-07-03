# Metric overview

This package provides two ways to evaluate the performance of an
embedding. The most generic is to consider “neighborhood preservation”:

1.  Find the *k*-nearest neighbors of a point in the input space.
2.  Do the same in the output space.
3.  Count the overlap.

This only requires calculating the Euclidean distance matrix in the
input and output spaces, without any other assumptions about the type of
dimensionality reduction carried out. This package provides some quality
measures based around this concept.

Alternatively, if some sort of labelling is applied to the points, each
point can be treated as the target in a retrieval procedure:

1.  Rank all the other points by distance to the target point.
2.  See how highly in the ranked list the points with the same label are
    found.
3.  Construct a Receiver Operating Characteristic (ROC) curve, or
    something like it (e.g. Precision-Recall curve)
4.  Calculate the Area Under the Curve (AUC).
5.  Average over all points.

This only needs the output distance matrix, but requires the sort of
labeling usually reserved for data intended for supervised
classification. Quadra can also provide some help with this, but
requires the [PRROC package](https://cran.r-project.org/package=PRROC)
to be installed.

`quadra` also contains random triplet and random pair distance methods
for global structure preservation. See the global preservation article
for those metrics.
