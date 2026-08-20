# In_silico_breeding
Genomic mating code with input genotype (VCF) and phenotype file of a founder population and conduct genomic mating to determine best crossing plans. 
Genotype must be phased and imputed (BEAGLE). 
Phenotype may ideally be BLUP.

Using PopVar and AlphaSimR.

Import Genotype (vcf).R
- Code for importing phased genotype into R ready matrix format.
- Saves as RDS file for future use.

PopVar Genomic Mating.R
- Code for prediction of genetic variance in bi-parental populations without simulation of progeny genotype.
- Provides crossing simulation of all unique crosses in founder population. 
- Relatively slow.

(V1) Genomic Selection Model Training.R
- Code for training genomic selection model based on genotype and phenotype.
- Saves model and marker effects.
- Includes 5 fold CV for model comparison.

(V2.1) Progeny Simulation Genomic Mating.R
- Code for simulation of progeny genotype and prediction of phenotype. 
- Provides crossing simulation of all unique crosses in founder population. 
- Requires saved genomic selection model. 

(V2.2) Progeny Simulation Crossing (Planned)
- Code for simulation of progeny genotype and prediction of phenotype. 
- Provides crossing simulation of SPECIFIED crossing plan. 
- Requires saved genomic selection model. 
