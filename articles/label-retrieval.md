# Label retrieval

If some sort of labelling is applied to the points, each point can be
treated as the target in a retrieval procedure:

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

Rows with undefined AUC values, such as rows that cannot define both
positive and negative examples, are excluded from overall and per-label
averages. If no rows remain for an average, that average is `NA_real_`.

## Further Reading

Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
Information retrieval perspective to nonlinear dimensionality reduction
for data visualization. *Journal of Machine Learning Research*,
*11*(Feb), 451-490.
[http://www.jmlr.org/papers/v11/venna10a.html](http://www.jmlr.org/papers/v11/venna10a.md)

### Precision-Recall AUC and Receiver Operating Characteristic Area Under the Curve

Davis, J., & Goadrich, M. (2006, June). The relationship between
Precision-Recall and ROC curves. In *Proceedings of the 23rd
international conference on Machine learning* (pp. 233-240). ACM.
<http://pages.cs.wisc.edu/~jdavis/davisgoadrichcamera2.pdf>

Keilwagen, J., Grosse, I., & Grau, J. (2014). Area under
precision-recall curves for weighted and unweighted data. *PloS One*,
*9*(3), e92209. <https://dx.doi.org/10.1371/journal.pone.0092209>
