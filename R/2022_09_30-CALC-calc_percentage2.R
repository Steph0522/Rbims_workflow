#' @title Calculate the percentage.
#' @description  Calculate the percentage of KO in certain pathway.
#' @usage calc_percentage(tibble_ko, y_axis, data_experiment=NULL)
#' @param tibble_ko a tibble object from mapping_ko.
#' @param y_axis a character, indicating the pathway to analyze.
#' @param data_experiment optional. a data frame object 
#' containing metadata information.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @examples
#' calc_percentage(ko_bin_mapp, Pathway)    
#' @export
calc_percentage2<-function(tibble_ko,
                          y_axis,
                          data_experiment=NULL){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
  # Select data -------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway", 
                    "KO", "dbCAN", "domain_name", "Pfam")
  
  # Transform from wide to long -------------------------------------------####
  numeric_table<- tibble_ko %>%
    group_by({{y_axis_enquo}}) %>% 
    summarise_if(is.numeric, sum) %>% 
    mutate_at(c(1), ~replace(., is.na(.), "Unassigned"))%>% 
    column_to_rownames(y_axis_label) %>% 
    dplyr::select_at(vars(matches("bin"))) 
  
percentage_table<-as.data.frame(t(t(numeric_table)/colSums(numeric_table)))*100
Table_with_percentage<- percentage_table %>% 
  rownames_to_column(var = y_axis_label) %>% 
  pivot_longer(cols = -!!y_axis_enquo,
               values_to = "Percentage", 
               names_to = "Bin_name")
  
  # Join data experiment --------------------------------------------------####
  if(is.null(data_experiment) == F){
    Table_with_percentage<-Table_with_percentage %>%
      left_join(data_experiment, by="Bin_name")
  }
  
  return(Table_with_percentage)
}