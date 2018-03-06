<!-- markdown - embed externally stored table
What am I calling Vollenweider's k? -->

# Does Lake and Stream Connectivity Control Phosphorus Retention in Lakes?

## Abstract

Lake water residence time and depth are known to be strong predictors of phosphorus (P) retention. However, there is substantial variation in P retention among lakes with the same depth and residence time. One potential explanatory factor for this remaining variation is that connectivity among lakes as well as connectivity between lakes and nearby streams influences either lake P trapping or the type of stream delivered P. Therefore, we examined the extent to which connectivity among lakes versus connectivity between lakes and streams contributes to differences in P retention among lakes. Specifically, we evaluated the effect of lake-based and stream-based connectivity metrics on P retention using a hierarchical parameterization of the Vollenweider equation. We compared Vollenweider’s k in lakes with high and low stream versus lake connectivity. We found that variation in k is more strongly associated with lake connectivity metrics compared to stream connectivity metrics. This result suggests that lake-associated processes, which likely control P trapping, play a larger role in determining P retention than stream-associated processes which may affect variation in the type of delivered P.

## Introduction

### Importance of P retention

A comprehensive understanding of P cycling at broad scales is necessary to predict nutrient concentrations among many different lake types and to better manage the risks of eutrophication from excess nutrient loading (Smith 2003). These predictions can ultimately be used to implement management actions aimed at avoiding phenomena such as algal blooms and fish kills which can negatively affect human health and degrade ecosystem services. One of the primary challenges of predicting nutrient concentrations is that lake nutrient concentrations are not simply a function of total inputs to the lake, but rather a function of complex interactions between lake characteristics and biogeochemical cycles. For example, lake phosphorus concentration is a function of biological response to phosphorus loading yet this response is not a simple function of load but also depends on lake characteristics such as depth and biogeochemical cycling with iron (Sondergaard 2003, Wagner et al. 2011). Such interactions may obscure direct relationships between nutrient concentration and total nutrient inputs.

Although both nitrogen (N) and phosphorus loading contribute to eutrophication risk, studies on lake nutrient retention commonly focus solely on phosphorus (P). This focus on P retention is a result of strong evidence that P rather than N is the primary limiting nutrient in lakes (Schindler 1977). Although later studies found that N and N:P ratios strongly control primary production in some lakes, P concentration is strongly correlated with phytoplankton primary production in most lakes (Downing et al. 2001, Elser et al. 2007).

### Connectivity 

One of the critical steps in the P cycle occurs when P inputs to lakes are subject to processing as a result of uptake, burial, and resuspension (citation). The extent of such processing varies among lakes and is known as _P retention_. 

The individual components of in-lake P processing have been the subject of intense study at the individual lake level. **For example, citation A discovered this while citation B showed that.** However, process-based study of these details is not feasible for hundreds of lakes at the regional scale. 

### 

In the present study, we examined the extent to which connectivity among lakes versus connectivity between lakes and streams contributes to differences in P retention among lakes.

## Methods

### P Retention Models

Lake P retention is typically modelled as a non-linear function of water residence time using a set of equations known as the Vollenweider equations (Vollenweider 1975). Although there are several variants of these equations, we use a 2-parameter model that has been shown to have good performance in multiple cross-sectional studies (Brett and Benjamin 2007, Cheng et al. 2010). The equation takes the form:

\begin{equation}
R_{p} = 1 - \frac{1}{1 + k\tau^{x}}
\end{equation}

where $k$ and $x$ are adjustable (fitted) parameter coefficients and $\tau$ is water residence time. The value of $k$ has been variously interpreted as either a bulk loss rate or as an approximation of sedimentation rate (citation). We first used the non-hierarchical Eq. 1 to model the overall relationship between P retention and water residence time for all lakes in the study dataset. Next, we explored hierarchical versions of Eq. 1 with respect to $k$ and specific sub-populations of the overall lakes dataset $k_{j}$:

\begin{equation}
R_{p} = 1 - \frac{1}{1 + k_{j}\tau^{x}}
\end{equation}

\begin{equation}
k_{j} = a_{i}
\end{equation}

We assigned lake membership in each sub-population $a_{i}$ based on a variety of lake characteristics and connectivity metrics where the splitting criteria was determined from the results of fitting conditional inference trees. The conditional inference tree (ctree) technique is a machine learning method where binary splits of the independent variables are recursively repeated to find the split that maximizes association with the depedendent variable (citation). In the present study, our dependent variable was P retention while our independent variable was a single lake characteristics or connectivity metric (Table 1). The advantage of the ctree technique over the more widely used classification-regression tree (CART) technique is that tree growth stopping rules are pre-specified. As a result, this avoids some of the subjectivity associated with post-hoc tree pruning.

We fit both non-hierarchical and hierarchical models using the non-linear extension to the `brms` package to interface with the Stan statistical program (Buekner, Stan Development Team 2017). In both models, we set an informative prior on $k$ and $x$ of N(1.3, 0.1) and N(0.45, 0.1) respectively. These priors were based on the confidence intervals presented in Brett and Benjamin (2007) and qualitatively matched those used by Cheng et al. (2010). We used the default settings of `brms` to generate posterior estimates which specify running four chains of 2,000 iterations each with no thinning and random intial values between -2 and 2.

### Connectivity metrics

We calculated several stream connectivity metrics including average link length, stream order ratio, and the closest distance to an upstream lake. Both average link length and stream order ratio are approximations of stream network complexity (citation, fractal stream book). Low connectivity lakes often have highly complex stream networks with a short average link length and a high stream order ratio. The first step in calculating average link length is to dissolve stream network nodes that do not occur at a fork (junction). Next, we divided the total stream length in the upstream watershed by the total number of stream reaches. We calculated (Strahler) stream order ratio by dividing the number of first-order streams in the upstream watershed of a given lake by the total number of higher order streams (citation). Finally, we calculated closest distance to an upstream lake on a path-distance basis (as the fish swims) rather than on a straight-line basis (as the crow flies).

In addition to stream connectivity metrics, we calculated several lake connectivity metrics including total upstream lake area, number of upstream lakes, and presence/abscence of an upstream lake. For the purposes of calculating lake connectivity metrics, we defined a lake as any waterbody with an area greater than 4 ha (0.04 km2). In general terms, a low connectivity lake is one with low total upstream lake area and little to none upstream lakes.

We calculated both lake and stream connectivity metrics at multiple scales (Figure 1). First, we calculated connectivity metrics at the scale of individual lake watersheds (iws scale). We defined an individual lake watershed as the area draining to a particular lake exclusive of any upstream areas that drain into a lake greater than or equal to 10 ha (0.1 km2). We also calculated connectivity metrics at the scale of entire upstream lake networks (nws scale). Network watersheds are defined as the area draining into any part the upstream network regardless of the presence or abscence of upstream lakes (Figure 1).

![Diagram of the individual lake watershed and network watershed scales. In this example, the lake-watershed for lake 3 includes the lake-watershed of lake 2 because lake 2 is smaller than 10 ha.](figures/iws_nws.png)

We calculated all stream connectivity metrics using the `streamnet` and `nhdR` packages (citations). The algorithms in the `streamnet` package use the `sf` R package as well as the `v.net` and `v.stream.order` modules (Jasiewicz and Metz 2011) included in GRASS GIS (GRASS Development Team 2017). We calculated lake connectivity metrics using data from the LAGOSNE dataset (citation). Both sets of connectivity metrics were calculated using the National Hydrography Dataset (NHD) as a primary input. We used data on lake P retention and residence time from the National Eutrophication Survey (citations). All processed data and code are available at [DOI]. 

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

A high stream order ratio indicates a high complexity of streams in the watershed. P inputs to a lake with a high stream order ratio have likely undergone less in-stream processing as high order streams are known to have long/short "spiral distances"(citation).