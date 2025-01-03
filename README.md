# Pathway analysis
The Pathway Response Gene Sets (PRGS) was composed of 15 signaling pathways, each with some responsive genes. These genes were consistently deregulated across experiments and specific to the perturbed pathway.

The overall process flow of the PRGS was as follows:

1, The pathway response gene set was pre-treated. (preproccess.R)

2, 15 groups of pathway response genes were obtained. (Ppdd.R)

3, For each pathway, a Gaussian transformation was performed on each dataset. (Gaussian_transformation.R)

4, we combined the 15 groups of pathway response genes into a pathway-responsive gene set (PRGS). (PRGS.R)

