# Initiation of smoking and susceptibility to type 2 diabetes: a mendelian randomization study

Christopher S Thom, Zhuoran Ding, Benjamin F Voight

Corresponding author: Benjamin F Voight

If you have questions or would like to discuss any of this: thomc@email.chop.edu


# Summary statistics from the following GWAS were used as part of this study
 We used data from individuals of European ancestry only. All data sets were analyzed in genome build hg19. 

Smoking traits

Liu, M. et al. Association studies of up to 1.2 million individuals yield new insights into the genetic etiology of tobacco and alcohol use. Nature Genetics 51, 237–244 (2019)



Type 2 diabetes

Klarin, D. et al. Genetics of blood lipids among ~300,000 multi-ethnic participants of the Million Veteran Program. Nat. Genet. 50, 1514–1523 (2018)



Coronary artery disease 

van der Harst, P. & Verweij, N. Identification of 64 Novel Genetic Loci Provides an Expanded View on the Genetic Architecture of Coronary Artery Disease. Circ. Res. 122, 433–443 (2018) 


Body mass index

Yengo L, Sidorenko J, Kemper KE, Zheng Z, Wood AR, Weedon MN, Frayling TM, Hirschhorn J, Yang J, Visscher PM, GIANT Consortium. (2018). Meta-analysis of genome-wide association studies for height and body mass index in ~700,000 individuals of European ancestry. Hum Mol Genet. 2018 Oct 15;27(20):3641-3649. doi: 10.1093/hmg/ddy271.



# Scripts and Data set explanation - PDFs for forest plots and scatter plots, xls files for tables, are included here:

   MR-smk_int.T2D.R - the MR analyses of smoking initiation on T2D and CAD
   
   MR-smk_ces.T2D.R  - the MR analyses of smoking cessation on T2D (and CAD)
   
   MR-cig_day.T2D.R - the MR analyses of cigarettes per day on T2D (and CAD)
   
   summary_stats_forestplots.R - for generating all the forestplots
   
   instrument_summary_tables.R - for pulling out all the instruments (See SupplementalTables.xlsx)
    
   Table1.xls - Example calculations for converting effect sizes (beta) and standard errors into 'change in type 2 diabetes odds ratio per 2-fold change in smoking initiation exposure risk' 
   
   SmkT2D.LDSC.Scripts - command line scripts and explanation for deriving LDSC results
   
   LDSCresults - full results of all trait combinations; subset are discussed in paper
   
   
