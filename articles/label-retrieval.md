# Label retrieval

Use this guide when known labels define what should be retrieved near
each observation. Choose ROC AUC for rank discrimination across relevant
and irrelevant candidates. Prefer PR AUC when relevant observations are
rare and precision among retrieved observations matters.

[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
and
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
treat each observation as a query, rank the remaining observations by
distance, and regard matching labels as relevant. These are supervised
retrieval summaries. Both functions require the [PRROC
package](https://cran.r-project.org/package=PRROC).

Perfect retrieval has AUC 1 for both measures. Random ROC AUC is 0.5.
The random PR AUC baseline is the prevalence of relevant observations
among the remaining candidates, so it changes with label imbalance.

``` r

library(quadra)

iris_x <- as.matrix(iris[, -5])
pca_iris <- stats::prcomp(iris_x, retx = TRUE, rank. = 2)$x
dout <- as.matrix(stats::dist(pca_iris))
```

``` r

roc <- roc_auc(dout, iris$Species)
pr <- pr_auc(dout, iris$Species)

auc_summary <- rbind(
  ROC = c(overall = roc$av_auc, unlist(roc$label_av)),
  PR = c(overall = pr$av_auc, unlist(pr$label_av))
)
round(auc_summary, 3)
```

    ##     overall setosa versicolor virginica
    ## ROC   0.930      1      0.902     0.889
    ## PR    0.876      1      0.822     0.808

Every Iris class has 50 observations, so each query has 49 relevant
candidates among the other 149 and a random PR baseline of about 0.329.
The overall Iris scores are above their random baselines. Setosa is
retrieved perfectly in this PCA embedding, while the lower versicolor
and virginica scores identify labels for which nearby observations are
less consistently from the same class. The difference between ROC and PR
values is not itself a winner-takes-all comparison: the measures use
different baselines and emphasize different retrieval behavior.

## Inspecting Query-Level Results

By default, each function returns `list(av_auc, label_av)`. `av_auc`
weights each defined query equally; `label_av` contains per-label means
in first-appearance order. Set `ret_extra = TRUE` to append the
row-aligned `query_auc` vector and its defined-query count.

``` r

roc_detail <- roc_auc(dout, iris$Species, ret_extra = TRUE)
head(round(roc_detail$query_auc, 3))
```

    ## [1] 1 1 1 1 1 1

``` r

roc_detail$n_defined
```

    ## [1] 150

All 150 Iris queries are defined. In other data, queries without both
relevant and irrelevant candidates are `NA_real_` and excluded from the
overall and per-label means. If no query is defined, the overall average
is `NA_real_`.

Both functions require a finite square distance matrix and one
nonmissing label per observation. The query itself is excluded from its
ranking. Labels without a defined query remain `NA_real_` in `label_av`.

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
