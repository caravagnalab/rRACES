#### CNA

Driver_cna_PCAWG_1_1_ <- read_csv("Driver_cna_PCAWG-1 (1).csv") %>%
  dplyr::arrange(chr, from, to) %>% select(!c('...1')) %>% select(chr, from, to, segment_id) %>%
  mutate(from=as.integer(from), to=as.integer(to))

write.table(Driver_cna_PCAWG_1_1_, file= 'CNA_Driver.csv', 
            col.names=F, sep='\t',row.names = FALSE)

# conversion through https://genome.ucsc.edu/cgi-bin/hgLiftOver

converted = read_tsv("DriverCNA_converted_coordinates.csv")[,1:4]
colnames(converted) = c('chr', 'from_hg38', 'to_hg38', 'segment_id')
converted$dupl = duplicated(converted$segment_id)
converted = converted %>% filter(dupl == F)

old_cna_drivers <- read_csv("Driver_cna_PCAWG-1 (1).csv") %>%
  dplyr::arrange(chr, from, to) %>% select(!c('...1')) %>%
  rename(from_hg19=from, to_hg19=to)

double_coords = left_join(converted, old_cna_drivers) %>% select(!dupl) %>%
  rowwise() %>% mutate(type=strsplit(type, '-')[[1]][1])

write.csv(double_coords, file = 'CNA_Drivers_PCAWG.csv')

### SNVs

drivers_PCAWG <- read_csv("drivers_PCAWG.csv") %>%
  dplyr::arrange(chr, from, to) %>% select(!c('...1'))


# exon_variant

find_substitution = function(dataset,n){
  row = dataset[n,]
  mut=strsplit(row$aa_mutation, ':')[[1]] 
  consequence = strsplit(row$consequence_type, ':')[[1]][1:length(mut)]
  mut_df = data.frame(
    'transcript_affected'= row$transcript_affected,
    'driver_gene'= row$driver_label,
    'substitutions'= mut,
    'consequence'= consequence
      ) %>% filter(mut!='NA')
  return(mut_df)
}
mt=lapply(1:nrow(drivers_PCAWG), find_substitution, dataset=drivers_PCAWG)
mt=Reduce(rbind, mt) %>% filter(!is.na(substitutions))
annotated_drivers= left_join(mt,drivers_PCAWG)

cons_to_keep= c('exon_variant','frameshift_variant',
                'stop_gained','missense_variant','disruptive_inframe_deletion',
                'inframe_deletion','stop_lost')
drivers = annotated_drivers %>% filter(consequence %in% cons_to_keep) %>% 
  mutate(driver_label= paste0(driver_gene, ' ', substitutions)) %>%
  select("chr","from","to","ref","alt","ttype",driver_label, driver_gene)

drivers$dupl = duplicated(drivers$driver_label)
drivers = drivers %>% filter(dupl == F)
drivers = drivers %>% select(!dupl) %>% 
  mutate(mut_type= 
           ifelse(ref %in% c('A','C', 'T', 'G') & alt  %in% c('A','C', 'T', 'G') , 'SNV', 'indel')
         )

write.csv(drivers, file = 'SNVs_and_indel_Driver.csv')

drivers_coord = drivers %>% mutate(driver_id= paste0(chr,':',from,':',to,':',driver_label)) %>%
  select(chr,from,to, driver_id)
write.table(drivers_coord, file= 'SNVs_and_indel_Driver_coords.csv', 
            col.names=F, sep='\t',row.names = FALSE)

# conversion through https://genome.ucsc.edu/cgi-bin/hgLiftOver

new_drivers_coords= read_tsv('SNV_and_indel_drivers_new_coords.csv') 
colnames(new_drivers_coords)= c('chr', 'from_hg38','to_hg38', 'driver_label', 'sub') 
new_drivers_coords=new_drivers_coords %>% mutate(driver_label= paste0(driver_label,' ', sub))

drivers = drivers %>% mutate(driver=driver_label) %>% 
  mutate(driver_label=paste0(chr,':',from,':',to,':',driver_label))
drivers = left_join(drivers, new_drivers_coords)

drivers = drivers %>%
  rename('from_hg19'='from','to_hg19'='to') %>% select(!c(driver_label,sub))

gene_symbols <- c("TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "BRAF", "PTEN", "PIK3CA", "ALK", "ERBB2",
                  "CDKN2A", "MYC", "NOTCH1", "PTPN11", "VHL", "RET", "KIT", "JAK2", "FLT3", "TERT",
                  "APC", "RB1", "NF1", "HRAS", "SMAD4", "PTCH1", "CTNNB1", "NF2", "FGFR2", "MET",
                  "MAP2K1", "MAP2K2", "MAX", "IDH1", "IDH2", "FGFR3", "NRAS", "CBL", "STK11", "MSH2",
                  "MSH6", "EPCAM", "RAD51", "CHEK2", "FLT3", "DNMT3A", "TET2", "IDH1", "IDH2", "KEAP1")

drivers = drivers %>% filter(driver_gene %in% gene_symbols)
write.csv(drivers, file='SNVs_indels_Drivers_PCAWG.csv')

