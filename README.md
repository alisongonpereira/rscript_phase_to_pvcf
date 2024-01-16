This R script is designed to facilitate genetic data analysis by processing VCF files in conjunction with output from the PHASE software. It takes as input a VCF file (Variant Call Format) that only contains the GT field and an output file from PHASE software (out_pairs). The script performs validation checks on these inputs to ensure correct formatting and readability.

Key Features:

   - Input Processing: The script requires a VCF file, which can be prepared using the bcftools command:
    
    bcftools annotate --remove 'FORMAT/[any],FORMAT/[additional],FORMAT/[fields]' multiple_fields.vcf -o gt_only.vcf

   - Phased VCF Generation: The main goal is to generate a phased VCF file from the provided inputs.
   - Quality Control Reporting: Reports samples with imputation and phasing probabilities under a user-defined value, set to 0.8 by default.

This how what GPT4 tried to illustrate the pipeline. I got a bit confused by it, not gonna lie.

![image](https://github.com/alisongonpereira/rscript_phase_to_pvcf/assets/101121606/395415fc-4950-4791-b019-8f87992d86cb)


Reference:

    Stephens, M., and Donnelly, P. (2003). "A comparison of Bayesian methods for haplotype reconstruction from population genotype data." American Journal of Human Genetics, 73:1162-1169.
    Stephens, M., Smith, N., and Donnelly, P. (2001). "A new statistical method for haplotype reconstruction from population data." American Journal of Human Genetics, 68, 978--989.
    Stephens, M., and Scheet, P. (2005). "Accounting for Decay of Linkage Disequilibrium in Haplotype Inference and Missing-Data Imputation." American Journal of Human Genetics, 76:449-462.
    Li, N., and Stephens, M. (2003). "Modelling Linkage Disequilibrium, and identifying recombination hotspots using SNP data." Genetics, 165:2213-2233.
    Crawford et al. (2004). "Evidence for substantial fine-scale variation in recombination rates across the human genome." Nature Genetics, 36: 700-706​

​​​.
