# Weakly Supervised Ordinal Support Vector Machine (WSO-SVM)
R implementation of the WSO-SVM algorithm. The main algorithm is in ` WSO-SVM Run.R`.

WSO-SVM is a data-inclusive ML model, which was specifically designed for a unique dataset of image-localized biopsies with spatially matched multiparametric MRI. It integrates all sources of data to train a robust model: (a) Labeled/biopsy samples with precise labels class 1 and 2. (b) Unlabeled samples from tumoral region, class 1 or 2 (unknown). (c) Normal brain samples from normal brain, assigned to class 3. And these three classes obey the intrinsic order.

Simply, it leverages a combination of data sources including precisely labeled samples and imprecisely labeled samples for three ordinal classes.

For the general weakly supervised ordinal learning, please refer to the model 'Hybrid Ordinal' Learner with the code via link: https://github.com/lujiawang13/Hybrid-Ordinal-Learner


We appreciate it if you would please cite our original paper if you found the package useful for your work.
The original paper: (will add soon)

