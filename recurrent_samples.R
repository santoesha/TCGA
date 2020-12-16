#Distribute primary+recurrent tumors 
characteristics <- read_tsv(file='../TCGA/gdc_sample_sheet.2020-07-07.tsv')
recurrent <- characteristics %>% dplyr::filter(`Sample Type` == "Recurrent Tumor")
recurrent <- intersect(recurrent$`Case ID`,characteristics$`Case ID`)

write_delim(as.data.frame(recurrent[1:4]),"set1_recurrent_samples.txt")
write_delim(as.data.frame(recurrent[5:9]),"set2_recurrent_samples.txt")
write_delim(as.data.frame(recurrent[10:13]),"set3_recurrent_samples.txt")

