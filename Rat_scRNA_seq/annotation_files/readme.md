panglao: mm_genelist.xlsx

### genelist on "Brain Cell Type Specific Gene Expression and Co-expression Network Architectures"
### based on 41598_2018_27293_MOESM3_ESM.xlsx
### https://ngdc.cncb.ac.cn/celltaxonomy/celltype   as reference
```
library(Hmisc)
Marker <- list(Neutrophil=capitalize(tolower(c('Stfa2','Stfa2l1','Fcnb'))),
               Micro=capitalize(tolower(c('CCL4','CCL3','CSF1R','CX3CR1','P2RY12'))),
               Endo=capitalize(tolower(c('APOLD1','FLT1','RGS5','PTPRB','TM4SF1','CD34'))),
               Astro=capitalize(tolower(c('AQP4','GJA1','GJB6','SLC4A4','SLC1A2'))),
               Neu=capitalize(tolower(c('RELN','CNR1','GAD2','SYT1','GAD1'))),
               Oli=capitalize(tolower(c('PLP1','MOBP','CLDN11','MBP','UGT8','Olig2','Olig1'))),
               Opcs=capitalize(tolower(c('PDGFRA','TNR','PCDH15','SHC4','VCAN'))),
               Bcell=capitalize(tolower(c('cd19','Ms4a1','Cd79b','Cd79a'))),
               Tcell=capitalize(tolower(c('CD3E','CD3D','cd8a','Cd4','CD3G'))),
               NK=capitalize(tolower(c('Klrb1','Klra1','Cma1'))),
               Pericyte=capitalize(tolower(c('Acta2','Tbx18','Slc38a11','Postn','Itga7'))),
               Mastercell=capitalize(tolower(c('Fcer1a'))),
               Macroph=capitalize(tolower(c('Fcna','Ccl12','Ms4a7'))),
               Pericyte=capitalize(tolower(c('FCGR3A','MS4A7','Slc38a11','Postn','Itga7'))),
)
```
