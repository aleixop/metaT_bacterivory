library(tidyverse)
library(vroom)
library(viridis)
library(patchwork)
library(ggtext)

select <- dplyr::select

theme_set(theme_minimal(base_size = 12))

scientific_10 <-
  function(x) {   parse(text=gsub("e\\+*", "0^", scales::scientific_format()(x))) }

split_header <- function(x,pos){
  map_chr(.x = x,
          .f = ~ str_split(string = .x,
                           pattern = '_')[[1]][pos])
}

expand_annotations <- function(annotation, group, simplify = TRUE){
  
  possible_groups <- 
    colnames(annotation)
  
  if (group == 'COG_category'){
    sep <- ''
  } else if (group %in% possible_groups){
    sep <- ','
  } else {
    stop(paste0('Group has to be one of the following: ', paste0(possible_groups, collapse = ' ')))
  }
  
  expanded_annotation <- 
    annotation %>% 
    mutate({{group}} := str_split(string = .data[[group]],pattern = sep)) %>% 
    unnest({{group}})
  
  if (simplify == TRUE){
    expanded_annotation <- 
      expanded_annotation %>% 
      dplyr::select(1,all_of({{group}}))
  }
  
  if (group == 'eggNOG_OGs'){
    expanded_annotation <- 
      expanded_annotation %>% 
      mutate(eggNOG_OGs = str_remove(eggNOG_OGs, pattern = '@.*$')) %>% 
      unique()
  }
  
  if (group == 'KEGG_Pathway'){
    expanded_annotation <- 
      expanded_annotation %>% 
      mutate(KEGG_Pathway = str_remove_all(KEGG_Pathway, 'ko|map')) %>% 
      unique()      
  }
  
  return(expanded_annotation)
  
}

collapse_counts_by_annotation <- function(annotation = annotations_min2TPM,
                                          counts = quantification_raw_min2TPM, 
                                          group,
                                          normalized = F,
                                          simplify = T,
                                          method = 'none'){
  
  expanded_annotation <- 
    expand_annotations(annotation = annotation,
                       group = group,
                       simplify = simplify)
  
  if (sum(grepl('[Ss]ample',colnames(counts)))){ # if '[Ss]ample' in colnames means the data is in tidy format
    counts_tidy <- 
      counts
  } else {
    message('Transforming to tidy format...')
    counts_tidy <- 
      counts %>% 
      pivot_longer(names_to = 'Sample',
                   values_to = 'Abundance',
                   -1)
  }
  
  counts_collapsed <- 
    counts_tidy %>% 
    left_join(expanded_annotation) %>% 
    group_by(.data[[group]], Sample) %>% 
    summarise(Abundance = sum(Abundance))
  
  if (!method %in% c('DESeq2','edgeR')){
    return(counts_collapsed)
  } else {
    counts_wide <- 
      counts_collapsed %>% 
      mutate({{group}} := case_when(.data[[group]] == '-' ~ 'not_annotated',
                                    is.na(.data[[group]]) ~ 'not_annotated',
                                    TRUE ~ .data[[group]])) %>% 
      group_by(.data[[group]], Sample) %>% 
      summarise(Abundance = sum(Abundance)) %>% 
      pivot_wider(names_from = Sample,
                  values_from = Abundance) %>% 
      column_to_rownames({{group}})
  }
  
  if (method == 'DESeq2'){
    require(DESeq2)
    dummy_metadata <- 
      data.frame(row.names = colnames(counts_wide),
                 foo = rep(1,ncol(counts_wide)))
    obj <-
      DESeqDataSetFromMatrix(countData = round(counts_wide),
                             colData = dummy_metadata,
                             design = ~ 1)
    if (normalized == TRUE){
      rld <- 
        rlog(obj)
      
      norm_counts <-
        rld %>% 
        assay()
      
      colnames(norm_counts) <- 
        colnames(rld)
    } else{
      return(obj)
    }
  } else if (method == 'edgeR'){
    require(edgeR)
    obj <- 
      DGEList(counts = counts_wide)
    if (normalized == TRUE){
      norm_counts <- 
        obj %>%
        calcNormFactors(method = 'TMM') %>% 
        cpm(normalized = T)
    } else {
      return(obj)
    }
  }
  return(norm_counts)
} 

hclust_order_from_tidy <- function(df, value_col, sample_col, id_col, clust.method = 'ward.D2', dist.method = 'euclidean', data_only = F){
  
  clust_out <- 
    df %>% 
    dplyr::select(value_col, sample_col, id_col) %>% 
    pivot_wider(names_from = sample_col,
                values_from = value_col,
                id_cols = id_col,
                values_fill = 0) %>% 
    column_to_rownames(var = id_col) %>% 
    dist(method = dist.method) %>% 
    hclust(method = clust.method)
  
  order <- 
    clust_out$labels[clust_out$order]
  
  if (data_only){
    return(clust_out)
  } else {
    return(order)
  }
} 

read_annotation_files <- function(files, colnames = T){
  
  exp_names <- 
    str_remove(files, '.*/') %>% 
    str_remove('_.*')
  
  df <- 
    map_df(.x = exp_names,
           .f = ~ vroom(files[str_detect(files,.x)],
                        col_names = colnames) %>% 
             mutate(Experiment = str_to_title(.x)))
  
  return(df)
}

read_quantification_files <- function(files){
  
  exp_names <- 
    str_remove(files, '.*/') %>% 
    str_remove('_.*')
  
  df <- 
    map(.x = exp_names,
        .f = ~ vroom(files[str_detect(files,.x)]) %>%
          pivot_longer(names_to = 'Sample',
                       values_to = 'Abundance',
                       -1) %>%
          mutate(Experiment = str_to_title(.x))) %>% 
    bind_rows()
  
  return(df)
}

read_files_by_experiment <- function(files, colnames = T){
  
  exp_names <- 
    str_remove(files, '.*/') %>% 
    str_remove('_.*') 
  
  df_list <-
    exp_names %>% 
    set_names(
      map(.x = .,
          .f = ~ vroom(files[str_detect(files,.x)], col_names = colnames)),
    .)
  
  names(df_list) <- str_to_title(exp_names)
  
  return(df_list)
  
}

read_fasta_by_experiment <- function(files, type){
  
  require(Biostrings)
  
  exp_names <- 
    str_remove(files, '.*/') %>% 
    str_remove('_.*')  
  
  if (type == 'nucl'){
  
  fasta_list <-
    exp_names %>% 
    set_names(
      map(.x = .,
          .f = ~ readDNAStringSet(files[str_detect(files,.x)])),
      .)
  
  } else if (type == 'aa') {
    
    fasta_list <-
      exp_names %>% 
      set_names(
        map(.x = .,
            .f = ~ readAAStringSet(files[str_detect(files,.x)])),
        .)
    
  } else {
    stop("Type must be one of the following: 'aa','nucl'")
  }
  
  names(fasta_list) <- str_to_title(exp_names)
  
  return(fasta_list)
  
}

read_rds_by_experiment <- function(files){
  
  exp_names <- 
    str_remove(files, '.*/') %>% 
    str_remove('_.*')  
  
  df_list <-
    exp_names %>% 
    set_names(
      map(.x = .,
          .f = ~ readRDS(files[str_detect(files,.x)])),
      .)
  
  names(df_list) <- str_to_title(exp_names)
  
  return(df_list)
  
}

vegan_formatter <- function(df, row_names, sample_col = 'Sample', abun_col = 'Abundance', fill = 0){
  
  df_wide <- 
    df %>% 
    dplyr::select(sample_col, abun_col, row_names) %>% 
    pivot_wider(names_from = sample_col,
                values_from = abun_col, 
                values_fill = fill) %>% 
    column_to_rownames(row_names) %>% 
    t()
  
  return(df_wide)
  
}

nmds_from_tidy_df <- function(df, row_names, sample_col = 'Sample', abun_col = 'Abundance', fill = 0, stress = F, try = 20){
  
  require(vegan)
  set.seed(1)
  
  df_wide <- 
    vegan_formatter(df, sample_col = sample_col, abun_col = abun_col, fill = fill, row_names = row_names)
  
  nmds_result <- 
    metaMDS(comm = df_wide,
            distance = 'bray',
            autotransform = F,
            try = try)
  
  stress_value <- 
    round(nmds_result$stress, 3)
    
  nmds_df <- 
    scores(nmds_result)$sites %>% 
    as_tibble(rownames = sample_col)
  
  if (stress == T){
    return(stress_value)
  } else {
    return(nmds_df)
  }
}

pcoa_from_tidy_df <- function(df, row_names, sample_col = 'Sample', abun_col = 'Abundance', fill = 0, axis_perc = F){

  require(ape)
  set.seed(1)
  
  df_wide <- 
    vegan_formatter(df, sample_col = sample_col, abun_col = abun_col, fill = fill, row_names = row_names)
  
  dist <- 
    vegdist(df_wide, method = 'bray')
  
  pcoa <- 
    pcoa(dist)
  
  axis <- 
    c(round(100*pcoa$values[1,2],1),
      round(100*pcoa$values[2,2],1))
  
  cat(paste0('First axis explains ',axis[1],'% and second axis ',axis[2],'%.\n'))
  
  pcoa_df <- 
    as_tibble(pcoa$vectors[,1:2], rownames = sample_col)
  
  if (axis_perc == T){
    return(axis)
  } else {
    return(pcoa_df)
  }
  
}

  
pca_from_tidy_df <- function(df, row_names, sample_col = 'Sample', abun_col = 'Abundance', fill = 0, scale = F, axis_perc = F){
  
  require(tidyverse)
  
  df_wide <- 
    vegan_formatter(df, sample_col = sample_col, abun_col = abun_col, fill = fill, row_names = row_names)
  
  require(ggfortify)
  
  pca <- 
    prcomp(df_wide, scale. = scale)
  
  axis <- 
    c(round(100*summary(pca)$importance[2,1],1),
      round(100*summary(pca)$importance[2,2],1))
  
  cat(paste0('First axis explains ',axis[1],'% and second axis ',axis[2],'%.\n'))
  
  pca_df <-
    pca %>% 
    autoplot(label = T) %>% 
    ggplot_build() %>% 
    .$data %>% 
    .[[2]] %>% 
    dplyr::select({{sample_col}} := label, PC1 = x, PC2 = y) %>% 
    as_tibble()

  if (axis_perc == T){
    return(axis)
  } else {
    return(pca_df)
  }  
}

filter_by_explained_abundance <- function(df, abun_col, names_col, explained_abun){
  
  names_to_keep <- 
    df %>%
    group_by(across(names_col)) %>% 
    summarise(!!abun_col := sum(.data[[abun_col]])) %>% 
    arrange(desc(.data[[abun_col]])) %>% 
    mutate(cumulative = cumsum(.data[[abun_col]]/sum(.data[[abun_col]]))) %>% 
    filter(lag(cumulative, default = 0) <= explained_abun) %>% 
    pull({{names_col}}) %>% 
    unique()
  
  filt_df <- 
    df %>% 
    filter(.data[[names_col]] %in% names_to_keep) 
  
  return(filt_df)
}

factorize_by_abundance <- function(df, to_factor, abun_col, method){
  
  if (method == 'total'){
    
    order <- 
      df %>% 
      group_by(across(to_factor)) %>% 
      filter(.data[[abun_col]] == max(.data[[abun_col]])) %>% 
      arrange(desc(.data[[abun_col]]))
    
  } else if (method == 'max'){
    
    order <- 
      df %>% 
      group_by(across(to_factor)) %>% 
      summarise(total = sum(.data[[abun_col]])) %>% 
      arrange(-total)

  } else {
    stop("Method must be one of the following: 'total', 'max'")
  }
  
  df_factorized <- 
    df %>% 
    ungroup() %>% 
    mutate(across(to_factor, ~ factor(.x, levels = unique(pull(order,.x)))))
  
  return(df_factorized)
}

heatmap_plotter <- function(df, value_col, sample_col, id_col, label = F, 
                            clust.method = 'ward.D2', dist.method = 'euclidean', 
                            log = F, dendrogram = T, facet = F, pseudocount = F){
  
  require(patchwork)
  require(ggtree)
  
  if (log){
    if (pseudocount){
      df <- 
        df %>% 
        ungroup() %>% 
        mutate({{value_col}} := log10(.data[[value_col]] + min(.data[[value_col]][.data[[value_col]] > 0])))
    } else {
      df <- 
        df %>% 
        mutate({{value_col}} := log10(.data[[value_col]]))
    }
  }
  
  clust <- 
    hclust_order_from_tidy(df = df, 
                           value_col = value_col,
                           sample_col = sample_col,
                           id_col = id_col,
                           clust.method = clust.method,
                           dist.method = dist.method,
                           data_only = T)
  
  order <- 
    clust$labels[clust$order]
  
  df_ordered <- 
    df %>% 
    mutate({{id_col}} := factor(.data[[id_col]], levels = order))
  
  p_heatmap <- 
    df_ordered %>% 
    ggplot(aes_string(x = sample_col, y = id_col)) +
    geom_tile(aes_string(fill = value_col)) +
    scale_fill_viridis_c(option = 'plasma', name = NULL) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = NULL,
         x = NULL)
  
  if (facet != F){
    
    p_heatmap <- 
      p_heatmap +
      facet_grid(~ .data[[facet]], 
                 scales = 'free_x',
                 space = 'free')
  
  }
  
  if (dendrogram){
    
    p_dendro <- 
      ggtree(clust, ladderize = F)
      # geom_tiplab(align = T) +
      # scale_x_continuous(expand=expansion(0.2))
    
    labs_dendro <- 
      p_dendro$data %>% 
      arrange(y) %>% 
      filter(!is.na(label)) %>% 
      pull(label)
    
    if (!identical(labs_dendro, order)){
      stop('Labels in dendrogram and heatmap do not match.')
    }
    
    p_heatmap <-
      p_dendro + p_heatmap + plot_layout(widths = c(0.3,0.7))
  }
  
  if (label != F){
    
    labels <- 
      df_ordered %>% 
      dplyr::select(id_col, label) %>% 
      unique()
    
    p_labels <- 
      labels %>% 
      ggplot(aes_string(y = id_col, x = 0)) +
      geom_text(aes_string(label = label), hjust = 0) +
      theme_void() +
      xlim(0,2)
    
    return(p_heatmap | p_labels)
  } else {
    return(p_heatmap)
  }
}

calculate_fold_change <- function(df, sample1, sample2, sample_col = 'Sample', abun_col = 'Abundance'){
  
  df_FC <- 
    df %>% 
    filter(.data[[sample_col]] %in% c(sample1, sample2)) %>% 
    pivot_wider(names_from = sample_col,
                values_from = abun_col,
                values_fill = 0) %>% 
    mutate(FC = .data[[sample2]] / .data[[sample1]],
           log2FC = log2(.data[[sample2]]) - log2(.data[[sample1]]))
  
  return(df_FC)
}

N50 <- function(lengths){
  len.sorted <- rev(sort(lengths))
  N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1]
  return(N50)
}

list_pivot_longer <- function(table, names_to = 'Sample', values_to = 'Abundance', excluded_cols = 1){
  
  res <- 
    names(table) %>% 
    set_names(
      map(.x = .,
          .f = ~ table[[.x]] %>% 
            pivot_longer(names_to = 'Sample',
                         values_to = 'Abundance',
                         -all_of(excluded_cols))),
      .)
  
  return(res)
}

prepare_annotations_cp <- function(test_genes, annotation, group, annotations_aux, description_col, output){
  
  annotation_filt <- 
    annotation %>% 
    expand_annotations(group = group) %>% 
    filter(.data[[group]] != '-', !is.na(.data[[group]])) %>% 
    left_join(annotations_aux)
  
  genes <- # keep only annotated genes in genes to test
    annotation_filt %>% 
    filter(Name %in% test_genes) %>% 
    pull(Name)
  
  term2gene <- 
    annotation_filt %>% 
    select(all_of(group), Name)
  
  term2name <- 
    annotation_filt %>% 
    select(all_of(c(group,description_col)))
  
  if (output == 'genes'){
    return(genes)
  } else if (output == 'annotation'){
    return(annotation_filt)
  } else if (output == 'term2gene'){
    return(term2gene)
  } else if (output == 'term2name'){
    return(term2name)
  }
}

enrichment_analysis <- function(test_genes, annotation, group, annotations_aux, description_col){
  
  require(clusterProfiler)
  
  input_objects <-
    c('genes','term2gene','term2name') %>% 
    set_names(
      map(.x = .,
          .f = ~ prepare_annotations_cp(
            test_genes = test_genes,
            annotation = annotation,
            group = group,
            annotations_aux = annotations_aux,
            description_col = description_col,
            output = .x
          )
      )
      ,.)
  
  enrichment <- 
    enricher(gene = input_objects$genes,
             TERM2GENE = input_objects$term2gene,
             TERM2NAME = input_objects$term2name) %>% 
    as_tibble()
  
}

go_enricher <- function(ontology, annotation, de_genes){
  
  require(topGO)  
  
  genes <- 
    tibble(gene = names(annotation),
           de = factor(as.integer(gene %in% de_genes))) %>% 
    deframe()
  
  topgo_object <-
    new("topGOdata",
        ontology = ontology,
        allGenes = genes,
        annot = annFUN.gene2GO,
        gene2GO = annotation)  
  
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  
  test_object <- 
    topgo_object %>% 
    getSigGroups(test.stat)
  
  GenTable(topgo_object,
           classic = test_object,
           topNodes = 1000,
           orderBy = 'weight') %>%
    filter(classic <= 0.05) %>% 
    as_tibble()
}

gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
clr = function(x, base= exp(1)){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

import_clust_results <- function(data_dir = 'data/clust/', group, obj){
  
  names_dict <- 
    c(pfam = 'PFAM_id',
      ko = 'KEGG_ko',
      og = 'eggNOG_OGs')
  
  results_dir <- 
    list.files(paste0(data_dir,group),
               pattern = 'Results',
               full.names = T)
  
  
  
  results_file <- 
    list.files(results_dir,
               pattern = 'Clusters_Objects.tsv',
               full.names = T)
  
  clust_results <- 
    read_tsv(results_file) %>% 
    filter(across(1) != 'Genes') %>% 
    pivot_longer(names_to = 'Cluster',
                 values_to = names_dict[[group]],
                 everything()) %>% 
    filter(!is.na(.data[[names_dict[[group]]]])) %>% 
    mutate(Cluster = str_remove(Cluster, ' \\(.*'))
  
  if (obj == 'clusters'){
    
    return(clust_results)
    
  } else if (obj == 'table'){
    
    results_files <- 
      list.files(paste0(results_dir,'/Processed_Data'),
                 full.names = T)
    
    tables <- 
      read_files_by_experiment(files = results_files)
    
    tables_tidy <- 
      names(tables) %>% 
      set_names(
        map(.x = .,
            .f = ~ tables[[.x]] %>%
              dplyr::rename(!!names_dict[[group]] := 'Genes') %>% 
              pivot_longer(names_to = 'Sample',
                           values_to = 'Abundance',
                           -1) %>% 
              mutate(Sample = factor(Sample, levels = metadata$Sample)) %>% 
              left_join(clust_results)),
        .)  
    
    return(tables_tidy)
  } else {
    stop("Object must be one of: 'clusters', 'table'")
  }
}

create_hulls <- function(df, group, x, y){
  
  df %>% 
    group_by(.data[[group]]) %>% 
    dplyr::slice(chull(.data[[x]],.data[[y]]))
  
}

kegg_brite_finder <- function(ko){
  
  require(KEGGREST)
  
  brites <- 
    keggGet(ko)[[1]]$BRITE
  
  cats_values <- 
    set_names(c('X',LETTERS[1:6]),0:6)
  
  tab <- 
    brites %>% 
    as_tibble() %>% 
    mutate(cat = cats_values[as.character(str_count(str_remove(value,'[^ ].*')))],
           Text = trimws(value)) %>% 
    select(-value)
  
  cats <- unique(tab$cat)
  
  cats_num <- 
    set_names(1:length(cats), cats)
  
  for (i in cats){
    
    tab <- 
      tab %>% 
      mutate({{i}} := case_when(cats_num[[i]] == cats_num[cat] ~ Text,
                                cats_num[[i]] > cats_num[cat] ~ '-',
                                TRUE ~ NA_character_))
    
  }
  
  res <- 
    tab %>% 
    fill({{cats}}) %>% 
    mutate(KEGG_ko = paste0('ko:',trimws(str_extract(Text, pattern = 'K.[0-9]*? '))),
           BRITE = str_remove(str_remove(X,'.*\\[BR:'),'\\]')) %>% 
    filter(str_detect(Text,'^K[0-9]')) %>% 
    select(BRITE,KEGG_ko, everything(), -cat, -Text, -X)
  
  return(res)
  
}

subset_fasta <- function(names, fasta, out_file = F, add_suffix = F){
  
  require(Biostrings)
  
  names(fasta) <- 
    str_remove(names(fasta), ' .*')
  
  fasta_filt <- 
    fasta[names(fasta) %in% names]
  
  if (add_suffix != F & length(fasta_filt) > 0){
    names(fasta_filt) <- 
      paste0(names(fasta_filt),'_',add_suffix)
  }
  
  if (out_file != F){
    writeXStringSet(fasta_filt, filepath = out_file)
  }
  
  return(fasta_filt)
}

experiment_bind_rows <- function(df){
  
  experiments <- c('Mar18','Jul17','Sep20','Nov18')
  
  df %>% 
    bind_rows(.id = 'Experiment') %>% 
    mutate(Experiment = factor(Experiment, levels = experiments))
  
  
}

is.outlier <- function(x){
  
  q1 <- quantile(x, 0.25, na.rm = T)
  q3 <- quantile(x, 0.75, na.rm = T)
  iqr <- q3 - q1
  
  !between(x =x,
           left = q1-1.5*iqr,
           right = q3+1.5*iqr)
  
}

simple_lca <- function(df, tax_cols, id_col = 'Name'){
  
  # tax cols should be from lowest to highest (Domain -> Species)
  
  df_grouped <- 
    df %>% 
    group_by(across(all_of(id_col)))    
  
  df_unique_hits <- 
    df_grouped %>% 
    filter(.data[[id_col]] %in% .data[[id_col]][n() == 1])
  
  df_multiple_hits <- 
    df_grouped %>% 
    filter(.data[[id_col]] %in% .data[[id_col]][n() > 1])   
  
  df_filtered <- 
    df_multiple_hits
  
  for (col in tax_cols){
    
    df_filtered <- 
      df_filtered %>%
      mutate({{col}} := case_when(n_distinct(.data[[col]]) > 1 ~ NA_character_,
                                  TRUE ~ .data[[col]]))
    
  }
  
  final_df <- 
    df_filtered %>% 
    unique() %>% 
    bind_rows(df_unique_hits)
  
  return(final_df)
  
}

read_blast <- function(df){
  
  read_tsv(df, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  
}

experiment_set_names <- function(){
  
  experiments <- c('Jul17','Mar18','Nov18','Sep20')
  
  experiments %>% 
  set_names()
  
}

percentage <- function(x){100*x/sum(x)}

ggtext_color <- function(text, color){
  paste0("<span style = 'color:",color,";'>",text,"</span>")
}

quantification_maker <- function(clusters){
  
  clusts_quant_dfs <- 
    experiment_set_names() %>% 
    map(~ quantification_tpm_dfs[[.x]] %>% 
          left_join(annotation_dfs[[.x]] %>% select(Name, KEGG_ko)) %>% 
          expand_annotations('KEGG_ko', simplify = F) %>% 
          filter(KEGG_ko %in% clusters$KEGG_ko) %>% 
          left_join(clusters) %>% 
          select(-KEGG_ko) %>% 
          unique())
  
  return(clusts_quant_dfs)
}

quantification_samples_category <- function(clusts_quants_dfs, samples = growth_samples){
  
  category_sample_quant_dfs <- 
    experiment_set_names() %>% 
    map(~ clusts_quants_dfs[[.x]] %>% 
          filter(Sample %in% samples) %>% 
          group_by(Sample, group) %>% 
          summarise(Abundance = sum(Abundance)))      
  
  return(category_sample_quant_dfs)
}


fold_change_category <- function(clusts_quant_dfs){
  
  fc_category <- 
    experiment_set_names() %>% 
    map(~ clusts_quant_dfs[[.x]] %>%
          group_by(Sample, group) %>% 
          summarise(Abundance = sum(Abundance)) %>% 
          left_join(metadata) %>% 
          group_by(group, State) %>% 
          summarise(Abundance = mean(Abundance)) %>%
          pivot_wider(names_from = State, values_from = Abundance, values_fill = 0) %>%
          mutate(fc_growth_lag = growth/lag,
                 fc_growth_decline = ifelse(.x == 'Jul17',NA_integer_,growth/decline)))
  
  return(fc_category)
}

add_fc_hk <- function(fc_category){
  
  fc_category_with_ref <- 
    experiment_set_names() %>% 
    map(.x = .,
        .f = ~ fc_category[[.x]] %>% 
          pivot_longer(names_to = 'State', values_to = 'mean_tpm', -c(group, contains('_'))) %>% 
          select(group, State, everything()) %>% 
          mutate(fc_growth_lag_ref = fc_hk_ref[[.x]][['fc_growth_lag']],
                 fc_growth_decline_ref = fc_hk_ref[[.x]][['fc_growth_decline']],
                 up_lag_growth = fc_growth_lag > fc_growth_lag_ref) %>% 
          pivot_wider(names_from = State, values_from = mean_tpm) %>% 
          select(group, !contains('_'), everything()))  
  
  return(fc_category_with_ref)
}

plot_fc <- function(fc_category_with_ref, plot = T, remove_exp_labels = F, title = '', category_sample_quant_dfs){
  
  order_main_category <- 
    fc_category_with_ref %>% 
    experiment_bind_rows() %>% 
    summarise(de_sum = sum(up_lag_growth),
              growth = sum(growth)) %>% 
    arrange(-de_sum, -growth)
  
  p_fc_main_category_1 <- 
    fc_category_with_ref %>% 
    experiment_bind_rows() %>% 
    mutate(group = factor(group, levels = order_main_category$group),
           de = if_else(up_lag_growth == T, 'Yes','No'),
           de = factor(de, levels = c('Yes','No'))) %>% 
    ggplot(aes(x = Experiment, y = fct_rev(group))) +
    geom_point(aes(size = fc_growth_lag, color = de)) +
    labs(x = NULL, y = NULL, subtitle = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 9),
          legend.position = 'bottom',
          panel.grid.minor.x = element_blank()) +
    scale_size_area(name = 'FC', limits = c(0,10)) +
    scale_color_manual(name = 'FC > FC(HK)', 
                       values = viridis(option = 'mako', n = 9)[c(4,7)]) +
    guides(size = guide_legend(order = 1),
           color = guide_legend(override.aes = list(size=3)))
  
  p_fc_main_category_2 <-
    category_sample_quant_dfs %>% 
    experiment_bind_rows() %>% 
    mutate(group = factor(group, levels = order_main_category$group)) %>% 
    ggplot(aes(x = Abundance, y = fct_rev(group))) +
    geom_jitter(width = 0, height = 0.05, alpha = 0.5) +
    labs(x = 'TPM', y = NULL) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_log10()
  
  if (remove_exp_labels == F){
    p_fc_main_category <-
      p_fc_main_category_1 + p_fc_main_category_2 +
      plot_layout(widths = c(0.4,0.6))
  } else {
    p_fc_main_category <-
      (p_fc_main_category_1+theme(axis.text.x = element_blank())) + 
      p_fc_main_category_2 +
      plot_layout(widths = c(0.4,0.6))
  }
  
  if (plot == 'fc'){
    return(p_fc_main_category_1)
  } else if (plot == 'tpm'){
    return(p_fc_main_category_2)
  } else {
    return(p_fc_main_category)
  }
}

quant_filterer <- function(list_df){
  
  experiment_set_names() %>% 
    map(~ list_df[[.x]] %>% 
          mutate(total = rowSums(.[-1])) %>% 
          filter(total > 0) %>% 
          select(Name, metadata$Sample[metadata$Sample %in% colnames(.)]))
  
}

species_collapser <- function(list_df){
  
  experiment_set_names() %>% 
    map(~ list_df[[.x]] %>% 
          mutate(EukProt_ID = str_remove(Name, '_.*')) %>% 
          group_by(EukProt_ID) %>% 
          summarise(across(where(is.double), ~ sum(.x))))
  
}

plot_comparison_categories_trophic_mode <- function(cat, ncol, nrow){
  
  set.seed(3)
  
  quantification_sp_tmm_categories_sp25 %>% 
    experiment_bind_rows() %>% 
    filter(category %in% cat) %>% 
    mutate(group = factor(group, levels = order_groups)) %>%
    arrange(group) %>% 
    mutate(group2 = str_remove(group, 'H\\+-transporting ')) %>% 
    mutate(group2 = factor(group2, levels = unique(.$group2))) %>% 
    left_join(trophic_modes_final_species) %>% 
    mutate(trophic_mode = factor(trophic_mode, levels = c('Heterotroph','Mixotroph','Phototroph'))) %>% 
    left_join(metadata) %>% 
    ggplot(aes(x = trophic_mode, y = perc)) +
    geom_jitter(aes(color = trophic_mode),width = 0.2, height = 0) +
    geom_violin(aes(fill = trophic_mode), alpha = 0.7) +
    labs(y = '',
         x = '',
         subtitle = cat) +
    scale_color_manual(values = trophic_colors, name = '') +
    scale_fill_manual(values = trophic_colors, name = '') +
    facet_wrap(~ group2, scales = 'free_y',
               labeller = labeller(group2 = label_wrap_gen(25)),
               nrow = nrow,
               ncol = ncol) +
    theme(legend.position = 'bottom',
          axis.text.x = element_blank())
}

plot_comparison_categories_trophic_mode_state <- function(cat, wrap = 25){
  
  set.seed(3)
  
  quantification_sp_tmm_categories_sp25 %>% 
    experiment_bind_rows() %>% 
    filter(category %in% cat) %>% 
    mutate(group = factor(group, levels = order_groups)) %>%
    arrange(group) %>% 
    mutate(group2 = str_remove(group, 'H\\+-transporting ')) %>% 
    mutate(group2 = factor(group2, levels = unique(.$group2))) %>% 
    left_join(trophic_modes_final_species) %>% 
    mutate(trophic_mode = factor(trophic_mode, levels = c('Heterotroph','Mixotroph','Phototroph'))) %>% 
    left_join(metadata) %>% 
    ggplot(aes(x = trophic_mode, y = perc)) +
    geom_jitter(aes(color = trophic_mode),width = 0.2, height = 0) +
    geom_violin(aes(fill = trophic_mode), alpha = 0.7) +
    labs(y = '',
         x = '',
         subtitle = cat) +
    scale_color_manual(values = trophic_colors, name = '') +
    scale_fill_manual(values = trophic_colors, name = '') +
    ggh4x::facet_grid2(State ~ group2, 
                       scales = 'free_y', 
                       independent = 'y',
                       axes = 'all',
                       labeller = label_wrap_gen(wrap)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_blank())
}

calculate_perc_categories_groups <- function(all_df, sp_df){
  
  experiment_set_names() %>% 
    map(~ all_df[[.x]] %>% 
          filter(Name %in% annotation_categories$Name) %>% 
          left_join(annotation_categories) %>% 
          mutate(Sample = factor(Sample, levels = metadata$Sample)) %>% 
          group_by(EukProt_ID, Sample, group, category) %>% 
          summarise(Abundance = sum(Abundance)) %>% 
          left_join(sp_df[[.x]] %>% dplyr::rename('Total_abundance' = Abundance)) %>% 
          mutate(perc = 100*Abundance/Total_abundance) %>% 
          left_join(eukprot_tax %>% select(EukProt_ID, Name_to_Use)))  
  
}

overlap_facet_plotter <- function(df, title){ 
  ggplot(data = df,
         aes(x = KO_Name_b, y = fct_rev(KO_Name_a))) +
    geom_tile(aes(fill = perc)) +
    scale_fill_viridis_c(option = 'mako', name = '% shared', limits = c(0,100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = NULL, y = NULL) +
    labs(subtitle = title) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 25))
}
