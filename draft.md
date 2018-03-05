# Does Lake and Stream Connectivity Control Phosphorus Retention in Lakes?

## Abstract

Lake water residence time and depth are known to be strong predictors of phosphorus (P) retention. However, there is substantial variation in P retention among lakes with the same depth and residence time. One potential explanatory factor for this remaining variation is that connectivity among lakes as well as connectivity between lakes and nearby streams influences either lake P trapping or the type of stream delivered P. Therefore, we examined the extent to which connectivity among lakes versus connectivity between lakes and streams contributes to differences in P retention among lakes. Specifically, we evaluated the effect of lake-based and stream-based connectivity metrics on P retention using a hierarchical parameterization of the Vollenweider equation. We compared Vollenweider’s k in lakes with high and low stream versus lake connectivity. We found that variation in k is more strongly associated with lake connectivity metrics compared to stream connectivity metrics. This result suggests that lake-associated processes, which likely control P trapping, play a larger role in determining P retention than stream-associated processes which may affect variation in the type of delivered P.

## Introduction

A comprehensive understanding of P cycling at broad scales is necessary to manage the risks of eutrophication from excess nutrient loading (citation). One of the critical steps in the P cycle occurs when P inputs to lakes are subject to processing as a result of uptake, burial, and resuspension (citation). The extent of such processing varying amoung lakes and is known as _P retention_. The individual components of in-lake P processing have been the subject of intense study at the individual lake level. **For example, citation A discovered this while citation B showed that.** However, process-based study of these details is not feasible for hundreds of lakes at the regional scale. 

In the present study, we examined the extent to which connectivity among lakes versus connectivity between lakes and streams contributes to differences in P retention among lakes.

## Methods

### Connectivity metrics

We calculated several stream connectivity metrics including average link length, stream order ratio, and the closest distance to an upstream lake. Both average link length and stream order ration are approximations of stream network complexity (citation, fractal stream book). It is calculated by dissolving stream network nodes that do not occur at a "fork" and dividing the total stream length in the upstream watershed by the total number of stream reaches. Stream order ratio is calculated by dividing the number of first-order streams in the upstream watershed of a given lake by the total number of higher order streams. For this analysis we used Strahler stream order (citation). A high stream order ratio indicates a high complexity of streams in the watershed. P inputs to a lake with a high stream order ratio have likely undergone less in-stream processing as high order streams are known to have long/short "spiral distances"(citation). Closest distance to an upstream lake may covary with the type of P inputs. It is calculated on a path-distance basis (as the fish swims) rather than on a straight-line basis (as the crow flies). 

In addition to stream connectivity metrics we calculated several lake connectivity metrics including total upstream lake area, number of upstream lakes, and the presence of an upstream lake. 

We examined the influence of many of our connectivity metrics on P cycling both at the scale of an invidual lake watershed and at the scale of the watershed of the entire upstream lake network. 

### Vollenweider

First, we used the Vollenweider model to establish the overall relationship between P retention and water residence time among lakes included in our dataset. 

### Software

## Results

Lakes with both high lake connectivity and high stream connecitvity were most often found in the northern portions of the study area.

![Population-level fit the Vollenweider P retention model where the dark circles are empirically measured values of P retention and residence time, while the light circles and vertical bars are posterior means and 95% credible intervals respectively.](/home/jose/Documents/Science/Dissertation/Analysis/public/dissertation-analysis_files/figure-html/global_vollenweider-1.png)

**Table of conditional inference tree splits**

![Partitioned (hierarchical) fits of the Vollenweider P retention model where the lines and shaded polygons represent the posterior mean and 95% credible interval respectively. Darker lines and polygons represent the left side of the parition which is usually smaller in magnitude.](/home/jose/Documents/Science/Dissertation/Analysis/public/dissertation-analysis_files/figure-html/model_on_partitions-2.png)

![Comparisons of the posterior distribution of Vollenweider's k among lake sub-populations partioned according to Table 1 where the subscript 1 denotes a subpopulation on the left side of the partition which is usually smaller in magnitude.](/home/jose/Documents/Science/Dissertation/Analysis/public/dissertation-analysis_files/figure-html/model_on_partitions-1.png)

**Conceptual model**

**Lake examples?**

## Discussion

Stream connectivity is really stream _complexity_. 