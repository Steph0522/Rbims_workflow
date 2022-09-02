#' @title Extract abundance profile of InterProScan output.
#' @description Reads a table object created with InterProScan and generates
#' a profile table of abundance with the hits of the KEGG, PFAM or 
#' INTERPRO databases. The output of KEGG database can be used within 
#' mapping_ko.
#' @usage read_interpro(data_interpro, 
#' database=c("KEGG", "Pfam", "INTERPRO",
#' "TIGRFAM", "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
#' "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", 
#' "MobiDBLite","Hamap", "Gene3D", "Coils", "CDD"), profile = TRUE)
#' @param data_interpro a table, output of InterProScan on tsv format.
#' InterProScan should have been run with -pa option to be able to use the 
#' KEGG option, in the database argument.
#' @param database a character indicating for which database do you want to
#' get the abundance profile. Valid options are "KEGG", "PFAM" or "INTERPRO".
#' @param profile a logical value indicating if you want to print a profile 
#' or not. This option is valid for "PFAM" and "INTERPRO" database. 
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import tibble dplyr stringr tidyr janitor rlang
#' @examples
#' \dontrun{
#' read_interpro(data_interpro="inst/extdata/Interpro_test.tsv", database="INTERPRO", 
#' profile = F)
#' }
#' @export
read_picrust2<-function(data_picrust, 
                        profile=TRUE, write=FALSE,
                        database= c("KO", "EC", "pathway")){
# Extract functions or gene----------------------------------------------####
    table_picrust<-suppressWarnings(
      suppressMessages(read_delim(data_picrust,
                                  delim="\t", 
                                  col_names = T)  )) 
# Choosing database----------------------------------------------####
    
if(database=="KO"){
          colnames(table_picrust)[1] <- "KO"}
if(database=="EC"){
      colnames(table_picrust)[1] <- "EC"}   
if(database=="pathway"){
     colnames(table_picrust)[1] <- "pathway"}   
    
# Profile or not ----------------------------------------------####
    
 if(isFALSE(profile)){
      picrust<-table_picrust%>% 
        pivot_longer(-1,names_to= "Bin_name", 
                     values_to="Abundance")
      
    } else{
      picrust<-table_picrust
    }

# Write data or not --------------------------------------------------------------####

if(isTRUE(write)){
  write_tsv(picrust, paste0("picrust_output_", format(Sys.time(), "%b_%d_%X"), ".tsv"))
}
else{
  return(picrust)
}

# Return ----------------------------------------------------------------####
return(interpro)
}


