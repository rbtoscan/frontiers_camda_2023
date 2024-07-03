This repository contains all scripts and processed ouputs used for analysis and interpretation for the CAMDA 2023 challenge and the Frontiers Bioinformatics manuscript submission.

# Antimicrobial Resistance in Diverse Urban Microbiomes: Uncovering Patterns and Predictive Markers

### Authors
**Rodolfo Brizola Toscan<sup>1,\*</sup>, Wojciech Lesiński<sup>2</sup>, Piotr Stomma<sup>2</sup>, Balakrishnan Subramanian<sup>3</sup>, Paweł Łabaj<sup>2</sup> and Witold R. Rudnicki<sup>2,3</sup>**

1. Małopolska Centre of Biotechnology, Jagiellonian University, Cracow, Poland  
2. Faculty of Computer Science, University of Białystok, Białystok, Poland  
3. Computational Center, University of Białystok, Białystok, Poland

**Correspondence**:  
Rodolfo Brizola Toscan  
rodolfo.toscan@doctoral.uj.edu.pl

# Abstract
Antimicrobial resistance (AMR) poses a significant global health threat, exacerbated by urbanization and anthropogenic activities. This study investigates the distribution and dynamics of AMR within urban microbiomes from six major U.S. cities using metagenomic data provided by the CAMDA 2023 challenge. We employed a range of analytical tools to investigate sample resistome, virome, and mobile genetic elements (MGEs) across these urban environments. Our results demonstrate that AMR++ and Bowtie outperform other tools in detecting diverse and abundant AMR genes, with binarization of data enhancing classification performance. The analysis revealed that a portion of resistome markers are closely associated with MGEs, and their removal drastically impacts the resistome profile and the accuracy of resistome modeling. These findings highlight the importance of preserving key MGEs in resistome studies to maintain the integrity and predictive power of AMR profiling models. This study underscores the heterogeneous nature of AMR in urban settings and the critical role of MGEs, providing valuable insights for future research and public health strategies.

![mermaid-diagram-2024-07-02-215545](https://github.com/rbtoscan/frontiers_camda_2023/assets/87976680/bcae5296-0a9c-4cc3-8b64-d20cb2090bf2)
The diagram provides a high-level view of the workflow used in this study, illustrating the stages from data preparation through to prediction. The process includes quality control of sequencing data, various profiling methods for resistome, virome, and mobile genetic elements, followed by integration and analysis using feature selection and machine learning techniques.
