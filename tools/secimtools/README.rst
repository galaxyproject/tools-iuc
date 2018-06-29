SECIMTools 1.0.0 Galaxy Wrappers
================================

Synopsis
========

SECIMTools project aims to develop a suite of tools for processing of
metabolomics data, which can be run in a standalone mode or via Galaxy Genomics
Framework.


Motivation
==========

The SECIMTools a set of python tools that are available both as standalone and
wrapped for use in Galaxy. The suite includes a comprehensive set of quality
control metrics (retention time window evaluation and various peak evaluation
tools), visualization techniques (hierarchical cluster heatmap, principal
component analysis, linear discriminant analysis, modular modularity
clustering), basic statistical analysis methods (partial least squares -
discriminant analysis, analysis of variance), advanced classification methods
(random forest, support vector machines), and advanced variable selection tools
(least absolute shrinkage and selection operator LASSO and Elastic Net).

Automated Installation
======================

Galaxy should be able to automatically install the dependencies if the conda
resolver is used. A secimtools pypi package is also available.


Manual Installation
===================

For those not using Galaxy's automated installation from the ToolShed, put
the XML and Python files into the ``tools/secimtools`` folder and add the XML
files to your ``tool_conf.xml`` as usual. For example, use::

  <section id="secimtools" name="SECIM Tools 1.0.0)">
    <tool file="secimtools/anova_fixed.xml"/>
    <tool file="secimtools/bland_altman_plot.xml"/>
    <tool file="secimtools/blank_feature_filtering_flags.xml"/>
    <tool file="secimtools/coefficient_variation_flags.xml"/>
    <tool file="secimtools/compare_flags.xml"/>
    <tool file="secimtools/compound_identification.xml"/>
    <tool file="secimtools/data_normalization_and_rescaling.xml"/>
    <tool file="secimtools/distribution_features.xml"/>
    <tool file="secimtools/distribution_samples.xml"/>
    <tool file="secimtools/drop_flags.xml"/>
    <tool file="secimtools/hierarchical_clustering_heatmap.xml"/>
    <tool file="secimtools/imputation.xml"/>
    <tool file="secimtools/kruskal_wallis.xml"/>
    <tool file="secimtools/lasso_enet_var_select.xml"/>
    <tool file="secimtools/linear_discriminant_analysis.xml"/>
    <tool file="secimtools/log_and_glog_transformation.xml"/>
    <tool file="secimtools/macros.xml"/>
    <tool file="secimtools/magnitude_difference_flags.xml"/>
    <tool file="secimtools/mahalanobis_distance.xml"/>
    <tool file="secimtools/merge_flags.xml"/>
    <tool file="secimtools/modify_design_file.xml"/>
    <tool file="secimtools/modulated_modularity_clustering.xml"/>
    <tool file="secimtools/multiple_testing_adjustment.xml"/>
    <tool file="secimtools/mzrt_match.xml"/>
    <tool file="secimtools/partial_least_squares.xml"/>
    <tool file="secimtools/principal_component_analysis.xml"/>
    <tool file="secimtools/random_forest.xml"/>
    <tool file="secimtools/remove_selected_features_samples.xml"/>
    <tool file="secimtools/retention_time_flags.xml"/>
    <tool file="secimtools/run_order_regression.xml"/>
    <tool file="secimtools/scatter_plot_2D.xml"/>
    <tool file="secimtools/scatter_plot_3D.xml"/>
    <tool file="secimtools/standardized_euclidean_distance.xml"/>
    <tool file="secimtools/subset_data.xml"/>
    <tool file="secimtools/summarize_flags.xml"/>
    <tool file="secimtools/svm_classifier.xml"/>
    <tool file="secimtools/threshold_based_flags.xml"/>
    <tool file="secimtools/tool_conf.xml"/>
    <tool file="secimtools/ttest_single_group.xml"/>
    <tool file="secimtools/ttest.xml"/>
  </section>

You will also need to install the ``secimtools`` python package within the
Galaxy Environment.


History
=======

======= ============================================================== 
Version Changes
------- -------------------------------------------------------------- 
v1.0.0  - First set of wrappers published to the Galaxy ToolShed
======= ============================================================== 


Bug Reports
===========

You can open support requests at https://github.com/secimTools/SECIMTools/issues


License (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
