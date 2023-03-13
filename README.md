# ClusterGVis <img src="man/clusterGVis.png" align="right" height="200" />

<!-- badges: start -->

For better clusteing and visualizing the time series gene expression data of RNA-SEQ with ***fuzzy c-means*** algorithm from e1071 package or ***Kmeans*** from ComplexHeatmap package, the function **blockwiseModules** from [**WGCNA**](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) package outputs also can be accepted and can be used to produce graphs for visualization. Here supply [**ClusterGVis**](https://github.com/junjunlab/ClusterGVis) package to cluster and visualize time-serie gene expression data in a more concise and elegant way with one-step operation. You can also do enrichment analysis for each clusters with using **clusterProfiler** in [**enrichCluster**](https://github.com/junjunlab/ClusterGVis) function. ClusterGVis allows you to create ***publication-quality figures***.

Thanks for the contributions for [**clusterProfiler**](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [**Mfuzz**](https://www.bioconductor.org/packages/release/bioc/html/Mfuzz.html) and [**ComplexHeatmap**](https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html)!

<!-- badges: end -->

## Installation

You can install the development version of transPlotR like so:

``` r
# Note: please update your ComplexHeatmap to the latest version!
# install.packages("devtools")
devtools::install_github("junjunlab/ClusterGVis")
```

## Citation

> Jun Z (2022). *ClusterGVis: One-step to Cluster and Visualize Gene Expression Matrix.*  https://github.com/junjunlab/ClusterGVis, https://github.com/junjunlab/ClusterGVis/wiki/document

## Figures

![1669187565352](https://user-images.githubusercontent.com/64965509/203490013-b2b33188-4d16-4991-b87e-acee2479e643.png)

## Document

> - https://github.com/junjunlab/ClusterGVis/wiki/document
> - https://github.com/junjunlab/ClusterGVis/wiki/version-0.0.2

## Related blogs

> - [**ClusterGVis 对基因表达时间序列聚类和可视化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507094&idx=1&sn=7c2872e4e7d92f0f16831f9e3b13f6ca&chksm=c184e6e7f6f36ff10ec1e41b1e45e90ffe8f0918878a6045fe0471c77729ea6af5d7e14beb5b&token=503374955&lang=zh_CN#rd)
> - [**ClusterGVis 的问题解答及优化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507124&idx=1&sn=bea21af4c86246715aed0219d4478aea&chksm=c184e6c5f6f36fd3a41222b014dd35ceeba8f983258fc36c287eab2188d4cbe3956f1e041d63&token=503374955&lang=zh_CN#rd)
> - [**ClusterGVis 对指定 cluster 添加注释**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507173&idx=1&sn=8c384e0e8678d0b20086b31a3bc1fa70&chksm=c184e694f6f36f82489e1e514d68d3ad5e80577616d197c9f8784176601ac342b682ed02b9f9&token=503374955&lang=zh_CN#rd)
> - [**ClusterGVis 继续查缺补漏**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507228&idx=1&sn=472e8fc2a17041b94043ac79e2018903&chksm=c184e76df6f36e7b3a6e0691a4140ebf030d323b9128cc3ab4c05c97d4000523dc71982bbecb&token=503374955&lang=zh_CN#rd)
> - [**ClusterGVis 添加多个样本注释**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507410&idx=1&sn=c33809620a13392f420a9bc1160400ac&chksm=c184e7a3f6f36eb5286ca59cba9cc0f81bb2e1d8e6f769588faf5f324c945b8cfe8f50708780&token=133699415&lang=zh_CN#rd)
> - [**bugs 报告和修复**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507735&idx=1&sn=d8236c12a07beecc5d6c181b196a9a78&chksm=c184e566f6f36c7072f382be27259127b4fa9c0b1228c891f5cfc35869861b3d9b8f6e9b0824&token=139164705&lang=zh_CN#rd)
> - [**ClusterGVis 对接单细胞啦**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247508319&idx=1&sn=6fe9cc5ea16468de8f8b15f00cdf8d22&chksm=c1849b2ef6f31238be4e1e414179470ea3084a9f7c62a9e334bbedc4f3605ffe8d924cba50e6&token=1432898004&lang=zh_CN#rd)
