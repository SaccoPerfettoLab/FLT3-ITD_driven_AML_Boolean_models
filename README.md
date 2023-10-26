# Unveiling the signaling network of FLT3-ITD AML improves drug sensitivity prediction

## Abstract

Currently, the identification of patient-specific therapies in cancer is mainly informed by personalized genomic analysis. In the setting of acute myeloid leukemia (AML), patient-drug treatment matching fails in a subset of patients harboring atypical internal tandem duplications (ITDs) in the tyrosine kinase domain of the FLT3 gene. To address this unmet medical need, here we develop a systems-based strategy that integrates multiparametric analysis of crucial signaling pathways, patient-specific genomic and transcriptomic data with a prior-knowledge signaling network using a Boolean-based formalism. By this approach, we derive personalized predictive models describing the signaling landscape of AML FLT3-ITD positive cell lines and patients. These models enable us to derive mechanistic insight into drug resistance mechanisms and suggest novel opportunities for combinatorial treatments. Interestingly, our analysis reveals that the JNK kinase pathway plays a crucial role in the tyrosine kinase inhibitor response of FLT3-ITD cells through cell cycle regulation. Finally, our work shows that patient-specific logic models have the potential to inform precision medicine approaches.

Read the full paper on [eLife journal](https://doi.org/10.7554/eLife.90532.1)

## Repository

The repository has two main directories:

-   *Boolean model and simulations:* it contains the code to (i) generate and visualize the cell-derived FLT3-ITD Boolean models, (ii) perform simulations of both cells and patients. The code is divided in:
A) [Model generation](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/Boolean%20model%20and%20simulations/A.Model_generation.html)
B) [Model visualization](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/Boolean%20model%20and%20simulations/B.Model_visualization.html)
C) [In silico validation](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/Boolean%20model%20and%20simulations/C.In_silico_validation.html)
D) [Combinatory treatment inference](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/Boolean%20model%20and%20simulations/D.Combinatory_treatment_inference.html)
E) [Patients simulation on cell-derived Boolean models](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/Boolean%20model%20and%20simulations/E.Patients_simulations.html)

-   *FLT3 ITD patients data preprocessing:* it contains the code to process FLT3 ITD patients data.
1. [FLT3-ITD localization in 14 patients AML cohort](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/FLT3%20ITD%20patients%20data%20preprocessing/1.FLT3_ITD_annotation.html)
2. [Annotation of functional impact of other mutations](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/FLT3%20ITD%20patients%20data%20preprocessing/2.Mutation_annotation.html)
3. [RNA-seq processing](https://raw.githack.com/SaccoPerfettoLab/FLT3-ITD_driven_AML_Boolean_models/main/FLT3%20ITD%20patients%20data%20preprocessing/3.-RNA-seq-processing.html)
