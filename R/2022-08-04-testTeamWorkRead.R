#library("dplyr")
library("readr")

#reading data
data_dbcan<-"../inst/extdata/test_data/Prueba1_Bin_154_2_1.fna.faaoverview.txt"
dbcan_df<-read_delim(data_dbcan, delim="\t", col_names = F)%>%
#  drop_na(.data$X12) %>%
#  select(.data$X1, .data$X15) %>%
#  drop_na() %>%
#  distinct() %>%
#  separate_rows(.data$X15, sep="\\|") %>%
#  filter(str_detect(.data$X15, "KEGG"))