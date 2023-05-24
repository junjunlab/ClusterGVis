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

> Jun Zhang (2022). *ClusterGVis: One-step to Cluster and Visualize Gene Expression Matrix.*  https://github.com/junjunlab/ClusterGVis

## Figures

![image](https://user-images.githubusercontent.com/64965509/226291391-0a3b3d5f-f3ef-499e-9815-bf1abe9442e9.png)

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
> - [**听说你还想添加 KEGG 注释?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247508362&idx=1&sn=68c90f6b6cf328f220c2926e1ff68df6&chksm=c1849bfbf6f312eddc8839706dc81984fb6da270dec3505b351b69f775e49e92ea585b3c5d1c&token=1432898004&lang=zh_CN#rd)
> - [**ClusterGVis 对接 monocle2 拟时序热图**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247509140&idx=1&sn=3f46ed8760be054b173a60642d8fe608&chksm=c1849ee5f6f317f36c7cb7e496181c1d9b431cb0c7a674f8cd31e7960b8f0830c4e8e0490b7d&token=46732270&lang=zh_CN#rd)
