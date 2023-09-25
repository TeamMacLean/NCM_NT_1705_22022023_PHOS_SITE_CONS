
make_base_df <- function(sites=NA, lifestyles=NA) {
## Makes the basic dataframe from the raw data  
  if (is.na(sites)) {
    sites <-  here::here("results/merged_test_sites.csv")
  }
  
  if (is.na(lifestyles)) {
    lifestyles <- here::here("lib", "all_proteomes_lifestyles.csv")
  }
  
  df <- readr::read_csv(
    sites,
    col_types = readr::cols(
      phos_site_id = readr::col_character(),
      mo_protein_id = readr::col_character(),
      species_compared = readr::col_character(),
      has_ortholog = readr::col_logical(),
      best_hit_ortholog = readr::col_character(),
      has_hit_in_ortholog = readr::col_logical(),
      best_hit_ortholog_score = readr::col_double(),
      best_hit_ortholog_identity = readr::col_character(),
      num_p_sites_in_mo = readr::col_integer(),
      num_sites_matched_in_ortholog = readr::col_integer(),
      sites_found_in_ortho = readr::col_integer(),
      mo_peptide = readr::col_character(),
      mo_seq_in_hit = readr::col_character(),
      orth_seq_in_hit = readr::col_character(),
      annotated_seq = readr::col_character(),
      mod_pattern = readr::col_character()
    )
    
  ) |>
    dplyr::transmute(
      phos_site_id,
      mo_protein_id,
      species_compared,
      has_ortholog,
      best_hit_ortholog,
      num_p_sites_in_mo,
      num_sites_matched_in_ortholog,
      prop_found = 100 * num_sites_matched_in_ortholog / num_p_sites_in_mo
    )
  
  trophism_info <-
    readr::read_csv(lifestyles)  |> 
    dplyr::rename("Saprotroph" = "Saprophyte") |>
    dplyr::select(
      name,
      tag,
      Saprotroph,
      Endophyte,
      Symbiont,
      Commensal,
      Biotroph,
      Hemibiotroph,
      Necrotroph
    ) |>
    dplyr::filter(tag != "Mo8") |>
    tidyr::pivot_longer(
      cols = c(
        Saprotroph,
        Endophyte,
        Symbiont,
        Commensal,
        Biotroph,
        Hemibiotroph,
        Necrotroph
      ),
      names_to = "trophism",
      values_to = "trophism_val"
    ) |>
    dplyr::group_by(trophism, tag) |>
    dplyr::filter(trophism_val == TRUE) |>
    dplyr::summarise(trophism = dplyr::first(trophism),
              name = name) |>
    dplyr::arrange(tag)
  
  appressorium_info <-
    readr::read_csv(lifestyles)  |>
    dplyr::rename(
      "Hyaline" = "Hyaline Appressorium",
      "Compound" = "compound appressorium",
      "Melanized" = "Melanized Appressorium",
      "None" = "No appressorium"
    ) |>
    dplyr::select(name, tag, "Compound", "Hyaline", "Melanized", "None") |>
    dplyr::filter(tag != "Mo8") |>
    tidyr::pivot_longer(
      cols = c("Compound", "Hyaline", "Melanized", "None"),
      names_to = "appressorium",
      values_to = "appressorium_val"
    ) |>
    dplyr::group_by(appressorium, tag) |>
    dplyr::filter(appressorium_val == TRUE) |>
    dplyr::summarise(appressorium = dplyr::first(appressorium),
              name = name) |>
    dplyr::arrange(tag)
  
  pmk1_info <-
    readr::read_csv(lifestyles)  |>
    dplyr::rename("pmk1_path" = "Pmk1 pathogenicity") |>
    dplyr::filter(tag != "Mo8") |>
    dplyr::mutate("pmk1_path" = dplyr::if_else(pmk1_path, "PMK1 Required", "PMK1 Not Required")) |>
    dplyr::select(name, tag, "pmk1_path")
  
  
  df <-
    dplyr::left_join(df, appressorium_info, by = c("species_compared" = "tag")) |> 
    dplyr::distinct()
  
  return(df)
}


convert_ids <- function(df) {
  # convert site data from cluster in format "MGG_01258T0 [68-87]-2xPhospho [S7(100); T14(100)]" 
  # to site data in phos-cons format MGG_00111T0-S13-S16
  
  # Note that the phos sites without positions 
  # "MGG_04185T0 [225-250]-1xPhospho [S/T]" 
  # come out as MGG_04185T0-NA so are not mappable to the IDs in the phos-cons data
  # which are like MGG_00111T0-S13-S16
  df |> 
    dplyr::mutate(gene = substr(id, 1, 11),
           offset = stringr::str_extract(id, "\\[\\d+-\\d+\\]")
    ) |> 
    tidyr::separate_wider_delim(offset, names = c("start", "stop"), delim = "-") |> 
    dplyr::mutate(
      start = stringr::str_remove(start,"\\["),
      start = as.integer(start),
      stop = stringr::str_remove(stop, "\\]"),
      stop = as.integer(stop),
      sites = stringr::str_extract(id, "\\[[STY]\\d+.*\\]"),
      sites = stringr::str_remove_all(sites, "\\(\\d+\\)"),
      sites = stringr::str_remove_all(sites, "\\(\\d+\\.\\d+\\)"),
      sites = stringr::str_remove_all(sites, "\\["),
      sites = stringr::str_remove_all(sites, "\\]"),
      sites = stringr::str_extract_all(sites, "[STY]\\d+"),
      site_letters = purrr::map(sites, function(x){stringr::str_extract(x, "[STY]")}),
      site_numbers = purrr::map(sites, function(x){as.integer(stringr::str_extract(x, "\\d+"))}),
      proper_starts = purrr::map2(site_numbers, start, function(x,y){ (x + y)-1 }),
      site_and_start = purrr::map2(proper_starts, site_letters, function(x,y){paste0(y,x)}),
      sites = purrr::map_chr(site_and_start, function(x){paste0(x, collapse="-")}),
      new_id = paste(gene, sites, sep="-")
    )  |> 
    dplyr::select(id, gene, new_id)
}


get_cluster_number <- function(cluster_info=NA) {
  #relies on cluster numbers calculated in analysis/0001_heatmap_development.Rmd
  #these are the kmeans cluster numbers.
  if (is.na(cluster_info)) {
    cluster_info = here::here("analysis/cluster_gene_info.csv")
  }
  readr::read_csv(cluster_info) |> 
    dplyr::select(phos_site_id, cluster_number) |> 
    dplyr::distinct() |> 
    dplyr::filter(! is.na(cluster_number))
}



## given two phlyo objects, returns the clades (nodes) that are different between them

get_nodes_different <- function(x, y){
  key1 <- ape::makeNodeLabel(x, "md5sum")$node.label
  key2 <- ape::makeNodeLabel(y, "md5sum")$node.label
  mk12 <- match(key1, key2)
  mk21 <- match(key2, key1)
  
  list( 
        x =  which(is.na(mk12)) + ape::Ntip(x),
        y = which(is.na(mk21)) + ape::Ntip(y)
  )
}