library(R6)

QProMS <- R6::R6Class(
  classname = "QProMS",
  public = list(
    ####################
    # Input parameters #
    data = NULL,
    input_type = NULL,
    intensity_type = NULL,
    expdesign = NULL,
    color_palette = NULL,
    # parameters for data wrangling #
    filtered_data = NULL,
    valid_val_filter = NULL,
    valid_val_thr = NULL,
    pep_filter = NULL,
    pep_thr = NULL,
    rev = NULL,
    cont = NULL,
    oibs = NULL,
    ###############################
    # parameters for normalization #
    normalized_data = NULL,
    norm_methods = NULL,
    is_norm = FALSE,
    vsn_norm_run_once = FALSE,
    ############################
    # parameters for imputation #
    imputed_data = NULL,
    imp_methods = NULL,
    is_mixed = NULL,
    is_imp = FALSE,
    imp_run_once = FALSE,
    #################
    # parameters For Statistics #
    tested_condition = NULL,
    univariate = NULL,
    clusters_def = NULL,
    clusters_number = NULL,
    stat_table = NULL,
    fold_change = 1,
    p_adj_method = NULL,
    alpha_ttest = NULL,
    anova_table = NULL,
    cluster_table = NULL,
    alpha_anova = NULL,
    # Functions #
    loading_data = function(input_path, input_type){
      
      self$data <- data.table::fread(input = input_path) %>%
        tibble::as_tibble(.name_repair = janitor::make_clean_names)
      
      self$input_type <- input_type
    },
    define_colors = function(){
      n_of_color <- max(self$expdesign %>% dplyr::count(replicate) %>% dplyr::pull(n))
      self$color_palette <- head(viridis::viridis(n = n_of_color + 1, direction = 1), -1)
    },
    make_expdesign = function(start_with = "lfq_intensity_"){
      ## qui mettere tutti gli if in base all'intensity type
      
      self$expdesign <- self$data %>%
        dplyr::select(gene_names, dplyr::starts_with(start_with)) %>%
        tidyr::pivot_longer(!gene_names, names_to = "key", values_to = "intensity") %>%
        dplyr::distinct(key) %>%
        dplyr::mutate(label = stringr::str_remove(key, start_with)) %>%
        dplyr::mutate(condition = stringr::str_remove(label, "_[^_]*$")) %>%
        dplyr::mutate(replicate = stringr::str_remove(label, ".*_"))
      
      self$define_colors()
      
      if(self$input_type == "max_quant"){
        self$pg_preprocessing()
      }
    },
    define_tests = function(){
      conditions <-
        dplyr::distinct(self$expdesign, condition) %>% pull(condition)
      
      tests <-
        tidyr::expand_grid(cond1 = conditions, cond2 = conditions) %>%
        dplyr::filter(cond1 != cond2) %>%
        dplyr::mutate(test = paste0(cond1, "_vs_", cond2)) %>%
        dplyr::pull(test)
      
      return(tests)
    },
    pg_preprocessing = function(){
      ########################################################################
      #### This function prepare the proteing groups in the QProMS format ####
      #### and remove duplicates.                                         ####
      ########################################################################
      
      ### this firts part remove duplicate and missing gene names
      ### in proteinGroups.txt input
      
      ## Indentify all duplicate gene names 
      ## and add after __ the protein iD
      
      data <- self$data
      expdesign <- self$expdesign
      
      list_unique_gene_names <- data %>%
        dplyr::select(protein_i_ds, gene_names, id) %>%
        dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*")) %>%
        ## every protein gorups now have only 1 gene name associated to it
        dplyr::rename(unique_gene_names = gene_names) %>%
        janitor::get_dupes(unique_gene_names) %>%
        dplyr::mutate(unique_gene_names = dplyr::case_when(
          unique_gene_names != "" ~ paste0(unique_gene_names, "__",
                                           stringr::str_extract(protein_i_ds, "[^;]*")),
          TRUE ~ stringr::str_extract(protein_i_ds, "[^;]*"))) %>%
        dplyr::select(unique_gene_names, id)
      
      ## update data that now don't have dupe or missing spot
      data_unique <- dplyr::left_join(data, list_unique_gene_names, by = "id") %>%
        dplyr::mutate(gene_names = dplyr::case_when(
          unique_gene_names != "" ~ unique_gene_names, 
          TRUE ~ gene_names)) %>%
        dplyr::select(-unique_gene_names) %>%
        dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*"))
      
      ### this second part standardize the data in the right format
      
      data_standardized <- data_unique %>% 
        dplyr::select(
          gene_names,
          dplyr::all_of(expdesign$key),
          peptides,
          razor_unique_peptides,
          unique_peptides,
          reverse,
          potential_contaminant,
          only_identified_by_site
        ) %>%
        tidyr::pivot_longer(
          !c(gene_names,
             peptides,
             razor_unique_peptides,
             unique_peptides,
             reverse,
             potential_contaminant,
             only_identified_by_site),
          names_to = "key",
          values_to = "raw_intensity"
        ) %>%
        dplyr::inner_join(., expdesign, by = "key") %>%
        dplyr::mutate(raw_intensity = log2(raw_intensity)) %>%
        dplyr::mutate(raw_intensity = dplyr::na_if(raw_intensity, -Inf)) %>%
        dplyr::mutate(bin_intensity = dplyr::if_else(is.na(raw_intensity), 0, 1)) %>%
        dplyr::select(-key)
      
      self$data <- data_standardized
    },
    data_wrangling = function(valid_val_filter = "alog", valid_val_thr = 0.75, 
                              pep_filter = "peptides", pep_thr = 2, 
                              rev = TRUE, cont = TRUE, oibs = TRUE, rescue_cont = NULL) {
      
      ##############################################################
      #### this function is divided in 2 steps:                 ####
      #### the first apply filer specific to maxquant input.    ####
      #### the second part filer the data base on valid values. ####
      ##############################################################
      
      # store inputs for summary table
      self$valid_val_filter <- valid_val_filter
      self$valid_val_thr <- valid_val_thr
      self$pep_filter <- pep_filter
      self$pep_thr <- pep_thr
      self$rev <- rev
      self$cont <- cont
      self$oibs <- oibs
      
      
      
      # setup object parameters
      self$vsn_norm_run_once <- FALSE
      self$imp_run_once <- FALSE
      data <- self$data
      
      if (self$input_type == "max_quant"){
        ### pep filter puo essere:
        ## c("peptides", "unique", "razor")
        
        data_wrang <- data %>%
          dplyr::mutate(potential_contaminant = dplyr::case_when(
            gene_names %in% rescue_cont ~ "", TRUE ~ potential_contaminant)) %>%
          ## remove reverse, potentialcontaminant and oibs from data base on user input
          {if(rev)dplyr::filter(., !reverse == "+") else .} %>%
          {if(cont)dplyr::filter(., !potential_contaminant == "+") else .} %>%
          {if(oibs)dplyr::filter(., !only_identified_by_site == "+") else .} %>%
          ## filter on peptides:
          {if(pep_filter == "peptides"){dplyr::filter(., peptides >= pep_thr)}
            else if (pep_filter == "unique") {dplyr::filter(., unique_peptides >= pep_thr)}
            else {dplyr::filter(., razor_unique_peptides >= pep_thr)}}
      }else{
        data_wrang <- data
      }
      
      ## different type of strategy for filter missing data:
      ## c("alog", "each_grp", "total") alog -> at least one group
      
      filtered_data <- data_wrang %>%
        {if(valid_val_filter == "total")dplyr::group_by(., gene_names)
          else dplyr::group_by(., gene_names, condition)} %>%
        dplyr::mutate(miss_val = dplyr::n() - sum(bin_intensity)) %>%
        dplyr::mutate(n_size = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(gene_names) %>%
        ## rage compreso tra 0 e 100% espresso in valori tra 0 e 1
        {if(valid_val_filter == "alog") dplyr::filter(., any(miss_val <= round(n_size * (1 - valid_val_thr), 0)))
          else dplyr::filter(., all(miss_val <= round(n_size * (1 - valid_val_thr), 0)))} %>%
        dplyr::ungroup() %>%
        dplyr::select(gene_names, label, condition, replicate, bin_intensity, raw_intensity) %>% 
        dplyr::rename(intensity = raw_intensity)
      
      self$filtered_data <- filtered_data
      
    },
    plot_protein_counts = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      expdes <- self$expdesign
      
      p <- data %>%
        dplyr::group_by(label) %>%
        dplyr::summarise(counts = sum(bin_intensity)) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(., expdes, by = "label") %>%
        dplyr::mutate(replicate = as.factor(replicate)) %>%
        dplyr::group_by(condition) %>%
        echarts4r::e_charts(replicate, renderer = "svg") %>%
        echarts4r::e_bar(counts) %>%
        echarts4r::e_x_axis(name = "Replicates") %>%
        echarts4r::e_y_axis(name = "Counts") %>%
        echarts4r::e_tooltip(trigger = "item") %>%
        echarts4r::e_color(self$color_palette) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_y_axis(
          name = "Counts",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>%
        echarts4r::e_x_axis(
          name = "Replicate",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 14,
            lineHeight = 60
          )
        ) %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore", "dataView"))
      
      return(p)
    },
    plot_total_counts = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      
      counts <- data %>% 
        dplyr::count(label) %>% 
        dplyr::slice(1L) %>%
        dplyr::pull(n)
      
      
      echarts4r::e_chart(renderer = "svg") %>%
        echarts4r::e_gauge(
          as.numeric(counts),
          "Proteins",
          min = 0,
          max = 10000,
          progress = list(show = TRUE)
        ) %>%
        echarts4r::e_title(text = "Total number of proteins", subtext = "Identify across all dataset") %>%
        echarts4r::e_color("#440154")
    },
    plot_protein_coverage = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      
      p <- data %>%
        dplyr::group_by(gene_names) %>%
        dplyr::summarise(counts = sum(bin_intensity)) %>%
        dplyr::ungroup() %>%
        dplyr::select(counts) %>%
        table() %>%
        tibble::as_tibble() %>%
        dplyr::rename(occurrence = n) %>%
        echarts4r::e_charts(counts, renderer = "svg") %>%
        echarts4r::e_bar(occurrence) %>%
        echarts4r::e_y_axis(name = "Counts") %>%
        echarts4r::e_tooltip(trigger = "item") %>%
        echarts4r::e_color(self$color_palette) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_y_axis(
          name = "Counts",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore", "dataView"))
      
      return(p)
    },
    plot_cv = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      
      p <- data %>% 
        dplyr::group_by(gene_names, condition) %>% 
        dplyr::summarise(
          mean = mean(intensity, na.rm = TRUE),
          sd = sd(intensity, na.rm = TRUE),
          CV = round(sd / mean, 3)
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(condition) %>% 
        echarts4r::e_chart() %>% 
        echarts4r::e_boxplot(
          CV,
          colorBy = "data",
          outliers = FALSE,
          itemStyle = list(color = "#DADADA", borderWidth = 2)
        ) %>%  
        echarts4r::e_tooltip(trigger = "axis") %>% 
        echarts4r::e_title(text = "Sample CV", subtext = "Coefficient of variation") %>% 
        echarts4r::e_color(self$color_palette)
      
      return(p)
    },
    plot_missing_data = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      
      p <- data %>%
        dplyr::group_by(label) %>%
        dplyr::mutate(bin_intensity = dplyr::if_else(bin_intensity == 1, "Valid", "Missing")) %>%
        dplyr::count(bin_intensity) %>%
        tidyr::pivot_wider(id_cols = label, names_from = bin_intensity, values_from = n) %>%
        {if(ncol(.) == 2) dplyr::mutate(., Missing = 0)else . } %>%
        dplyr::ungroup() %>%
        dplyr::mutate(total = Valid + Missing) %>%
        dplyr::mutate(perc_present = paste0(round(Valid*100/total, 1), "%")) %>%
        dplyr::mutate(perc_missing = paste0(round(Missing*100/total, 1), "%")) %>%
        echarts4r::e_charts(label, renderer = "svg") %>%
        echarts4r::e_bar(Valid, stack = "grp", bind = perc_present) %>%
        echarts4r::e_bar(Missing, stack = "grp", bind = perc_missing) %>%
        echarts4r::e_x_axis(name = "Samples", axisLabel = list(interval = 0, rotate = 45)) %>%
        echarts4r::e_y_axis(name = "Counts") %>%
        echarts4r::e_tooltip(trigger = "item") %>%
        echarts4r::e_color(c("#21918c", "#440154")) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_y_axis(
          name = "Counts",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore", "dataView"))
      
      return(p)
    },
    plot_missval_distribution = function(){
      
      # controllare che ci sia il self$filtered_data
      data <- self$filtered_data
      
      p <- data %>%
        dplyr::group_by(gene_names) %>%
        dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                         missval = any(is.na(intensity))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(missing_value = dplyr::if_else(missval, "Missing", "Valid")) %>%
        dplyr::mutate(missing_value = factor(missing_value, levels = c("Valid", "Missing"))) %>%
        dplyr::group_by(missing_value) %>%
        echarts4r::e_charts(renderer = "svg") %>%
        echarts4r::e_density(
          mean,
          smooth = TRUE,
          areaStyle = list(opacity = 0),
          symbol = "none"
        ) %>%
        echarts4r::e_y_axis(
          name = "Densiry",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>%  
        echarts4r::e_x_axis(
          name = "log2 Intensity",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>%
        echarts4r::e_color(c("#21918c", "#440154")) %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore"))
      
      return(p)
    },
    normalization = function(norm_methods = "None"){
      
      # store inputs for summary table
      self$norm_methods <- norm_methods
      
      if(is.null(self$filtered_data)){
        print("error")
      }
      
      self$imp_run_once <- FALSE
      data <- self$filtered_data
      
      if(norm_methods == "None"){
        self$is_norm <- FALSE
      }else{
        if(!self$vsn_norm_run_once){
          
          self$vsn_norm_run_once <- TRUE
          
          ## convert tibble data into a matrix
          raw_matrix <- data %>%
            tidyr::pivot_wider(id_cols = gene_names,
                               names_from = label,
                               values_from = intensity) %>%
            tibble::column_to_rownames("gene_names") %>%
            as.matrix()
          
          ## Variance stabilization transformation on matrix
          vsn_fit <- vsn::vsn2(2 ^ raw_matrix, verbose = FALSE)
          norm_matrix <- vsn::predict(vsn_fit, 2 ^ raw_matrix)
          
          ## return a table with QProMS object format
          normalized_data <- norm_matrix %>%
            tibble::as_tibble(rownames = "gene_names") %>%
            tidyr::pivot_longer(cols = !gene_names,
                                names_to = "label",
                                values_to = "norm_intensity") %>%
            dplyr::full_join(data, .x, by = c("gene_names", "label")) %>%
            dplyr::mutate(intensity = norm_intensity) %>% 
            dplyr::select(-norm_intensity) %>% 
            dplyr::relocate(intensity, .after = last_col())
          
          self$normalized_data <- normalized_data
        }
        self$is_norm <- TRUE
      }
      
    },
    plot_distribution = function(){
      
      if(self$is_norm){
        data <- self$normalized_data
      }else{
        data <- self$filtered_data
      }
      
      p <- data %>%
        dplyr::mutate(intensity = round(intensity, 2)) %>%
        dplyr::group_by(condition, label) %>%
        echarts4r::e_charts(renderer = "svg") %>%
        echarts4r::e_boxplot(
          intensity,
          colorBy = "data",
          layout = 'horizontal',
          outliers = FALSE,
          itemStyle = list(borderWidth = 2)
        ) %>%
        echarts4r::e_tooltip(trigger = "item") %>%
        echarts4r::e_color(self$color_palette) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_y_axis(
          name = "Intensity",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 16,
            lineHeight = 60
          )
        ) %>% 
        echarts4r::e_toolbox_feature(feature = "saveAsImage")
      
      return(p)
    },
    plot_correlation = function(){
      # define if use normalize or row intensity
      if(self$is_norm){
        data <- self$normalized_data
      }else{
        data <- self$filtered_data
      }
      
      mat <- data %>%
        dplyr::select(gene_names, label, intensity) %>%
        tidyr::pivot_wider(names_from = label, values_from = intensity) %>%
        dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), .fns = ~ !is.na(.x))) %>%
        tibble::column_to_rownames("gene_names") %>%
        cor() %>% 
        round(digits = 2)
      
      p <- mat %>% 
        echarts4r::e_charts(renderer = "svg") %>%
        echarts4r::e_correlations(order = "hclust", visual_map = FALSE) %>%
        echarts4r::e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
        echarts4r::e_y_axis(axisLabel = list(interval = 0, rotate = 0), position = "right") %>%
        echarts4r::e_tooltip(trigger = "item", formatter = htmlwidgets::JS("
          function(params){
          return('X: ' + params.value[0] + '<br />Y: ' + params.value[1] + '<br />Value: ' + params.value[2])
          }")) %>%
        echarts4r::e_title("Correlation matrix", subtext = "Pearson correlation") %>%
        echarts4r::e_visual_map(
          min = min(mat),
          max = 1,
          bottom = 150,
          inRange = list(color = c("#440154", "#31688e", "#35b779"))
        ) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage"))
      
      return(p)
    },
    imputation = function(imp_methods = "mixed", shift = 1.8, scale = 0.3) {
      
      # store inputs for summary table
      self$imp_methods <- imp_methods
      
      # Define if use normalize or row intensity
      if(self$is_norm){
        data <- self$normalized_data
      }else{
        data <- self$filtered_data
      }
      
      if(imp_methods == "mixed"){
        self$is_mixed <- TRUE
      }else{
        self$is_mixed <- FALSE
      }
      
      if(imp_methods == "mixed" | imp_methods == "perseus"){
        self$is_imp <- TRUE
        if(!self$imp_run_once){
          
          self$imp_run_once <- TRUE
          
          if(self$is_mixed){
            data_mixed <- data %>%
              dplyr::group_by(gene_names, condition) %>%
              dplyr::mutate(for_mean_imp = dplyr::if_else((sum(bin_intensity) / dplyr::n()) >= 0.75, TRUE, FALSE)) %>%
              dplyr::mutate(mean_grp = mean(intensity, na.rm = TRUE)) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(imp_intensity = dplyr::case_when(
                bin_intensity == 0 & for_mean_imp ~ mean_grp,
                TRUE ~ as.numeric(intensity))) %>%
              dplyr::mutate(intensity = imp_intensity) %>% 
              dplyr::select(-c(for_mean_imp, mean_grp, imp_intensity))
            
            data <- data_mixed
          }
          ## this funcion perform classical Perseus imputation
          ## sice use random nomral distibution i will set a set.seed()
          set.seed(11)
          
          imputed_data <- data %>%
            dplyr::group_by(label) %>%
            # Define statistic to generate the random distribution relative to sample
            dplyr::mutate(
              mean = mean(intensity, na.rm = TRUE),
              sd = sd(intensity, na.rm = TRUE),
              n = sum(!is.na(intensity)),
              total = nrow(data) - n
            ) %>%
            dplyr::ungroup() %>%
            # Impute missing values by random draws from a distribution
            # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
            dplyr::mutate(imp_intensity = dplyr::case_when(
              is.na(intensity) ~ rnorm(total, mean = (mean - shift * sd), sd = sd * scale),
              TRUE ~ intensity
            )) %>%
            dplyr::mutate(intensity = imp_intensity) %>%
            dplyr::select(-c(mean, sd, n, total, imp_intensity))
          
          self$imputed_data <- imputed_data
          
        }
        
      }else{
        self$is_imp <- FALSE
      }
      
      
    },
    plot_imputation = function(){
      
      if(!self$is_norm & !self$is_imp){
        data <- self$filtered_data
      }else if(self$is_norm & !self$is_imp){
        data <- self$normalized_data
      }else{
        data <- self$imputed_data
      }
      
      p <- data %>%
        dplyr::group_by(condition) %>%
        echarts4r::e_charts(renderer = "svg") %>%
        echarts4r::e_density(
          intensity,
          smooth = TRUE,
          areaStyle = list(opacity = 0),
          symbol = "none"
        ) %>%
        echarts4r::e_color(self$color_palette) %>%
        echarts4r::e_theme("QProMS_theme") %>% 
        echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore"))
      
      return(p)
    }, 
    plot_pca = function(view_3d = FALSE){
      
      # verificare che ci sia
      data <- self$imputed_data
      
      ## generate a matrix from imputed intensiy
      mat <- data %>%
        dplyr::select(gene_names, label, intensity) %>%
        tidyr::pivot_wider(id_cols = "gene_names",
                           names_from = "label",
                           values_from = "intensity") %>%
        tibble::column_to_rownames("gene_names") %>%
        as.matrix()
      
      ## perform PCA
      
      pca <- prcomp(t(mat), center = TRUE, scale = TRUE) 
      
      ## calculate persentage of each PC
      pca_var <- pca$sdev^2
      pca_var_perc <- round(pca_var/sum(pca_var)*100, 1)
      
      ## create a data.frame for the first 3 PC
      pca_table <- data.frame(
        label = rownames(pca$x),
        x = pca$x[, 1],
        y = pca$x[, 2],
        z = pca$x[, 3]
      ) %>% 
        dplyr::left_join(self$expdesign, by = "label") 
      
      ## generate plot
      if(!view_3d){
        p <- pca_table %>%
          dplyr::group_by(condition) %>%
          echarts4r::e_chart(x, renderer = "svg") %>%
          echarts4r::e_scatter(y, symbol_size = c(10, 10), bind = replicate) %>%
          echarts4r::e_tooltip(
            trigger = "item",
            formatter = htmlwidgets::JS("
        function(params){
          return('Rep: ' + params.name);
        }
      ")
          ) %>%
          echarts4r::e_title(text = "PCA", subtext = "Principal component analysis") %>%
          echarts4r::e_x_axis(
            name = paste0("PC1 - ", pca_var_perc[1], " %"),
            nameLocation = "center",
            nameTextStyle = list(
              fontWeight = "bold",
              fontSize = 15,
              lineHeight = 50
            )
          ) %>%
          echarts4r::e_y_axis(
            name = paste0("PC2 - ", pca_var_perc[2], " %"),
            nameLocation = "center",
            nameTextStyle = list(
              fontWeight = "bold",
              fontSize = 15,
              lineHeight = 50
            )
          ) %>% 
          echarts4r::e_color(self$color_palette) %>% 
          echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore", "dataView"))
      }else{
        p <- pca_table %>%
          dplyr::group_by(condition) %>%
          echarts4r::e_chart(x) %>%
          echarts4r::e_scatter_3d(y, z, symbol_size = c(10, 10), bind = replicate) %>%
          echarts4r::e_tooltip(
            trigger = "item",
            formatter = htmlwidgets::JS("
        function(params){
          return('Rep: ' + params.name);
        }
      ")
          ) %>%
          echarts4r::e_title(text = "PCA", subtext = "Principal component analysis") %>%
          echarts4r::e_x_axis_3d(
            name = paste0("PC1 - ", pca_var_perc[1], " %"),
            nameLocation = "center",
            nameTextStyle = list(
              fontWeight = "bold",
              fontSize = 15,
              lineHeight = 50
            )
          ) %>%
          echarts4r::e_y_axis_3d(
            name = paste0("PC2 - ", pca_var_perc[2], " %"),
            nameLocation = "center",
            nameTextStyle = list(
              fontWeight = "bold",
              fontSize = 15,
              lineHeight = 50
            )
          ) %>%
          echarts4r::e_z_axis_3d(
            name = paste0("PC3 - ", pca_var_perc[3], " %"),
            nameLocation = "center",
            nameTextStyle = list(
              fontWeight = "bold",
              fontSize = 15,
              lineHeight = 50
            )
          ) %>%
          echarts4r::e_legend() %>% 
          echarts4r::e_color(self$color_palette) %>% 
          echarts4r::e_toolbox_feature(feature = c("saveAsImage", "restore"))
      }
      
      return(p)
    },
    tidy_vector = function(data, name) {
      
      tidy_vec <- data %>%
        tibble::as_tibble(rownames = NA) %>%
        tibble::rownames_to_column(var = "gene_names") %>%
        dplyr::rename(!!name := value)
      
      return(tidy_vec)
    },
    stat_t_test_single = function(data, test, fc, alpha, p_adj_method, ...){
      
      cond_1 <- stringr::str_split(test, "_vs_")[[1]][1]
      cond_2 <- stringr::str_split(test, "_vs_")[[1]][2]
      
      self$univariate <- TRUE
      self$fold_change <- fc
      
      mat <- data %>%
        dplyr::filter(condition == cond_1 | condition == cond_2) %>%
        dplyr::mutate(label_test = paste(condition, replicate, sep = "_")) %>%
        tidyr::pivot_wider(id_cols = "gene_names",
                           names_from = "label_test",
                           values_from = "intensity") %>%
        tibble::column_to_rownames("gene_names") %>%
        dplyr::relocate(dplyr::contains(cond_2), .after = dplyr::last_col()) %>%
        na.omit() %>% 
        as.matrix()
      
      a <- grep(cond_1, colnames(mat))
      b <- grep(cond_2, colnames(mat))
      
      p_values_vec <- apply(mat, 1, function(x) t.test(x[a], x[b], ...)$p.value)
      
      p_values <- p_values_vec %>%
        as_tibble(rownames = NA) %>%
        tibble::rownames_to_column(var = "gene_names") %>%
        dplyr::rename(p_val = value)
      
      fold_change <- apply(mat, 1, function(x) mean(x[a]) - mean(x[b])) %>% #metterlo in log2?
        as_tibble(rownames = NA) %>%
        tibble::rownames_to_column(var = "gene_names") %>%
        dplyr::rename(fold_change = value)
      
      p_ajusted <- p.adjust(p_values_vec, method = p_adj_method) %>% 
        as_tibble(rownames = NA) %>%
        tibble::rownames_to_column(var = "gene_names") %>%
        dplyr::rename(p_adj = value)
      
      stat_data <- fold_change %>% 
        dplyr::full_join(., p_values, by = "gene_names") %>% 
        dplyr::full_join(., p_ajusted, by = "gene_names") %>% 
        dplyr::mutate(significant = dplyr::if_else(abs(fold_change) >= fc & p_adj <= alpha, TRUE, FALSE)) %>% 
        dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_significant") := significant) %>% 
        dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_p_val") := p_val) %>% 
        dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_fold_change") := fold_change) %>% 
        dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_p_adj") := p_adj)
      
      return(stat_data)
    },
    stat_t_test = function(data, test, fc = 1, alpha = 0.05, p_adj_method = "BH") {
      
      if(length(test) == 1){
        stat_table_single <-
          self$stat_t_test_single(
            data = data,
            test = test,
            fc = fc,
            alpha = alpha,
            p_adj_method = p_adj_method
          )
        
        complete_stat_table <- 
          data %>%
          tidyr::pivot_wider(id_cols = "gene_names",
                             names_from = "label",
                             values_from = "intensity") %>%
          dplyr::left_join(stat_table_single, by = "gene_names")
        
      } else {
        stat_table_map <-
          purrr::map(
            .x = test,
            .f = ~ self$stat_t_test_single(
              data = data,
              test = .x,
              fc = fc,
              alpha = alpha,
              p_adj_method = p_adj_method
            )
          ) %>%
          purrr::reduce(dplyr::full_join, by = "gene_names")
        
        complete_stat_table <- 
          data %>%
          tidyr::pivot_wider(id_cols = "gene_names",
                             names_from = "label",
                             values_from = "intensity") %>%
          dplyr::left_join(stat_table_map, by = "gene_names")
      }
      
      
      self$stat_table <- complete_stat_table
      self$tested_condition <- test
      
    },
    stat_anova = function(data, alpha = 0.05){
      
      self$alpha_anova <- alpha
      self$univariate <- FALSE
      
      p_values_vec <- data %>%
        split(.$gene_names) %>%
        purrr::map_dbl(~ summary(stats::aov(intensity ~ condition, .x))[[1]][["Pr(>F)"]][[1]]) %>%
        data.frame() %>%
        tibble::rownames_to_column(var = "gene_names") %>%
        tibble::deframe()
      
      p_values <- p_values_vec %>%
        self$tidy_vector(name = "p_val")
      
      p_ajusted <- p.adjust(p_values_vec, method = "BH") %>% 
        self$tidy_vector(name = "p_adj")
      
      stat_data <- data %>% 
        tidyr::pivot_wider(id_cols = "gene_names", names_from = "label", values_from = "intensity") %>% 
        dplyr::full_join(., p_values, by = "gene_names") %>% 
        dplyr::full_join(., p_ajusted, by = "gene_names") %>% 
        dplyr::mutate(significant = dplyr::if_else(p_adj <= alpha, TRUE, FALSE))
      
      self$anova_table <- stat_data
    },
    e_arrange_list = function(list, rows = NULL, cols = NULL, width = "xs", title = NULL) {
      plots <- list
      
      if (is.null(rows)) {
        rows <- length(plots)
      }
      
      if (is.null(cols)) {
        cols <- 1
      }
      
      w <- "-xs"
      if (!isTRUE(getOption("knitr.in.progress"))) {
        w <- ""
      }
      
      x <- 0
      tg <- htmltools::tagList()
      for (i in 1:rows) {
        r <- htmltools::div(class = "row")
        
        for (j in 1:cols) {
          x <- x + 1
          cl <- paste0("col", w, "-", 12 / cols)
          if (x <= length(plots)) {
            c <- htmltools::div(class = cl, plots[[x]])
          } else {
            c <- htmltools::div(class = cl)
          }
          r <- htmltools::tagAppendChild(r, c)
        }
        tg <- htmltools::tagAppendChild(tg, r)
      }
      
      if (!isTRUE(getOption("knitr.in.progress"))) {
        htmltools::browsable(
          htmltools::div(
            class = "container-fluid",
            htmltools::tags$head(
              htmltools::tags$link(
                rel = "stylesheet",
                href = "https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css",
                integrity = "sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO",
                crossorigin = "anonymous"
              )
            ),
            htmltools::h3(title),
            tg
          )
        )
      } else {
        if (!is.null(title)) {
          htmltools::div(title, tg)
        } else {
          tg
        }
      }
    },
    plot_volcano_single = function(test, highlights_names, text_color, bg_color) {
      
      table <- self$stat_table %>% 
        dplyr::select(gene_names, starts_with(test)) %>% 
        rename_at(vars(matches(test)), ~ str_remove(., paste0(test, "_")))
      
      min_thr <- table %>% 
        dplyr::filter(significant) %>% 
        pull(p_val) %>% 
        max()
        
      
      left_line <-
        tibble(p_val = c(-log10(min_thr),-log10(min_thr), max(-log10(table$p_val))),
               fold_change = c(min(table$fold_change),-self$fold_change,-self$fold_change))
      
      right_line <-
        tibble(p_val = c(max(-log10(table$p_val)),-log10(min_thr),-log10(min_thr)),
               fold_change = c(self$fold_change, max(table$fold_change), self$fold_change))
      
      p <- table %>%
        dplyr::mutate(color = dplyr::case_when(fold_change > 0 & significant ~ "#cf4446", 
                                               fold_change < 0 & significant ~ "#0d0887",
                                               TRUE ~ "#e9ecef")) %>%
        dplyr::group_by(color) %>%
        dplyr::mutate(fold_change = round(fold_change, 2)) %>%
        dplyr::mutate(p_val = -log10(p_val)) %>%
        dplyr::mutate(p_val = round(p_val, 3)) %>%
        echarts4r::e_chart(fold_change, renderer = "svg") %>%
        echarts4r::e_scatter(p_val, legend = FALSE, bind = gene_names) %>%
        echarts4r::e_tooltip(
          formatter = htmlwidgets::JS(
            "
      function(params){
        return('<strong>' + params.name +
                '</strong><br />FC: ' + params.value[0] +
                '<br />p.val: ' + params.value[1])
                }
    "
          )
        ) %>%
        echarts4r::e_add_nested("itemStyle", color) %>%
        echarts4r::e_data(left_line, fold_change) %>%
        echarts4r::e_line(
          p_val,
          legend = FALSE,
          color = "#440154",
          symbol = "none",
          lineStyle = list(type = "dashed", width = .8)
        ) %>%
        echarts4r::e_data(right_line, fold_change) %>%
        echarts4r::e_line(
          p_val,
          legend = FALSE,
          color = "#440154",
          symbol = "none",
          lineStyle = list(type = "dashed", width = .8)
        ) %>%
        echarts4r::e_toolbox() %>%
        echarts4r::e_toolbox_feature(feature = "dataZoom") %>%
        echarts4r::e_toolbox_feature(feature = "saveAsImage") %>%
        echarts4r::e_title(text = test,
                           left = "center") %>%
        echarts4r::e_x_axis(
          name = "Fold_change",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 15,
            lineHeight = 50
          )
        ) %>%
        echarts4r::e_y_axis(
          name = "-log(Pvalue)",
          nameLocation = "center",
          nameTextStyle = list(
            fontWeight = "bold",
            fontSize = 15,
            lineHeight = 50
          )
        )
      
      if (!is.null(highlights_names)) {
        for (name in highlights_names) {
          highlights_name <- table %>%
            dplyr::filter(gene_names == name) %>%
            dplyr::mutate(p_val = -log10(p_val)) %>%
            dplyr::select(xAxis = fold_change,
                          yAxis = p_val,
                          value = gene_names) %>% as.list()
          
          p <- p %>%
            echarts4r::e_mark_point(
              data = highlights_name,
              silent = TRUE,
              label = list(color = text_color),
              itemStyle = list(color = bg_color)
            )
        }
      }
      
      return(p)
    },
    plot_volcano = function(test, highlights_names = NULL, text_color = "#440154", bg_color = "#fde725"){
      if(length(test)==1){
        p <- self$plot_volcano_single(
          test = test, 
          highlights_names = highlights_names,
          text_color = text_color, 
          bg_color = bg_color
        )
      }else{
        volcanos <- purrr::map(
          .x = test,
          .f = ~ self$plot_volcano_single(
            test = .x,
            highlights_names = highlights_names,
            text_color = text_color, 
            bg_color = bg_color
          )
        )
        p <- self$e_arrange_list(volcanos, cols = length(volcanos), rows = 1)
      }
      
      return(p)
    },
    result_matrix = function(){
      # questa dovrebbe essere una fuzione privata
      mat <- self$anova_table %>% 
        dplyr::filter(significant) %>% 
        dplyr::select(-c(p_val, p_adj, significant)) %>% 
        tibble::column_to_rownames("gene_names") %>% 
        as.matrix() 
      
      mat_scaled = t(apply(mat, 1, scale))
      colnames(mat_scaled) <- colnames(mat)
      
      return(mat_scaled)
    },
    plot_heatmap = function(n_cluster = 0, clustering_distance = "euclidean", clustering_method = "average"){
      
      expdesign <- self$expdesign
      mat <- self$result_matrix()
      
      condition_anno <- setNames(as.character(expdesign$condition),
                                 as.character(expdesign$label))
      
      condition_col <- setNames(self$color_palette,
                                pull(distinct(expdesign, condition)))
      
      hm = Heatmap(
        matrix = mat,
        name = "Scaled \nIntensity",
        clustering_distance_rows = clustering_distance,
        # c("euclidean", "pearson", "spearman", "kendall")
        clustering_method_rows = clustering_method,
        # c("single", "complete", "average", "median", "centroid")
        row_km = n_cluster,
        show_row_dend = TRUE,
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
          cond = condition_anno,
          col = list(cond = condition_col),
          show_legend = FALSE
        )
      ) 
      
      if(n_cluster == 0){
        self$clusters_def <- FALSE
      }else{
        self$clusters_def <- TRUE
        self$define_cluster(hm)
      }
      
      return(draw(hm))
      
      
    },
    define_cluster = function(heatmap){
      
      mat <- self$result_matrix()
      
      cluster_order <- ComplexHeatmap::row_order(heatmap)
      
      
      for (i in 1:length(cluster_order)) {
        if (i == 1) {
          clu <- row.names(mat[cluster_order[[1]],])
          out <-
            cbind(clu, paste("cluster", names(cluster_order)[i], sep = ""))
          colnames(out) <- c("gene_names", "cluster")
        } else {
          clu <- row.names(mat[cluster_order[[i]],])
          clu <-
            cbind(clu, paste("cluster", names(cluster_order)[i], sep = ""))
          out <- rbind(out, clu)
        }
      }
      
      anova_significant <- self$anova_table %>% 
        dplyr::select(c(gene_names, p_val, p_adj, significant))
      
      cluster_mat <- mat %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("gene_names") %>% 
        dplyr::left_join(as.data.frame(out), by = "gene_names") %>% 
        dplyr::left_join(anova_significant, by = "gene_names")
      
      self$cluster_table <- cluster_mat
      self$clusters_number <- cluster_mat %>% 
        dplyr::distinct(cluster) %>% 
        dplyr::pull()
      
    },
    summary_result_table = function(){
      # controllare che tutti i parametri non siano nulli
      
      gt_table <- tibble::tibble(
        Parameters = c(
          "Total Proteins",
          "Condition",
          "Replicate",
          "Valid value strategy",
          "Valid value threshold",
          "Peptides column",
          "Peptides threshold",
          "Reverse",
          "Contaminant",
          "Only Identify by Site" 
          # "Significant by t.test",
          # "Threshold p adjusted",
          # "Fold change",
          # "Significant by ANOVA",
          # "Threshold p adjusted"
        ),
        
        Value = as.character(c(
          nrow(self$stat_table),
          table(self$expdesign$replicate)[1],
          table(self$expdesign$condition)[1],
          self$valid_val_filter,
          paste0(self$valid_val_thr*100, " %"),
          self$pep_filter,
          paste0("≥ ", self$pep_thr),
          if(self$rev){"Removed"}else{"Not removed"},
          if(self$cont){"Removed"}else{"Not removed"},
          if(self$oibs){"Removed"}else{"Not removed"}
          # self$stat_table %>% dplyr::filter(significant) %>% nrow(),
          # paste0("≤ ", self$alpha_ttest),
          # paste0("±", self$fold_change),
          # self$anova_table %>% dplyr::filter(significant) %>% nrow(),
          # paste0("≤ ", self$alpha_anova)
        )),
        
        Description = c(
          "Total number of proteins quantified in whole experiment.", 
          "Number of condition defined by experimental design.", 
          "Number of replicates for each condition defined by experimental design.",
          "Strategy used for filter missing data: at least in one group (alog), in each group (each_grp) or total.", 
          "Trheshold used to filter missing data based on the strategy selected.",
          "Column present in proteinGroups.txt that is selected to apply the filter on minimun number of peptides.", 
          "Trheshold used to filter peptides colum seleted. Only protein with plus or equal numbers are reteined.",
          "Column present in proteinGroups.txt that define proteins derived from the reversed part of the decoy database. These should be removed.", 
          "Column present in proteinGroups.txt that define proteins found to be a commonly occurring contaminant. These should be removed.",
          "Column present in proteinGroups.txt that define proteins identified only by a modification site. These should be removed."
          # paste0("Number of significant proteins in ", self$tested_condition),
          # "Trheshold used to filter for significant in t.test analysis",
          # "Fold change trheshold used during the analysis",
          # "Number of significant proteins by ANOVA test", 
          # "Trheshold used to filter for significant in ANOVA analysis"
        )
      ) %>% gt::gt(rowname_col = "Parameters") %>% 
        gt::tab_header(
          title = gt::md("**Summary results**")
        ) %>% 
        gt::tab_row_group(
          label = gt::md("**Data wranglilng parameters**"),
          rows = 4:10
        ) %>% 
        # gt::tab_row_group(
        #   label = gt::md("**Statistics**"),
        #   rows = 11:15
        # ) %>% 
        gt::tab_row_group(
          label = gt::md("**Experimental design**"),
          rows = 1:3
        ) %>% 
        gt::tab_footnote(
          footnote = gt::md("From *MaxQuant* output."),
          locations = cells_body(columns = Description, rows = stringr::str_detect(Description, "proteinGroups.txt"))
        )
      
      return(gt_table)
    },
    gsea_go_single = function(test, org_db, ontology, pvalue_cutoff, max_row){
      
      table <- self$stat_table %>% 
        dplyr::select(gene_names, starts_with(test)) %>% 
        dplyr::rename_at(dplyr::vars(matches(test)), ~ stringr::str_remove(., paste0(test, "_")))
      
      gsea_vec <- table %>%
        dplyr::arrange(-fold_change) %>%
        dplyr::select(gene_names, fold_change) %>%
        tibble::deframe()
      
      gsea_result <- clusterProfiler::gseGO(
        geneList     = gsea_vec,
        OrgDb        = org_db,
        ont          = ontology,
        keyType      = 'SYMBOL',
        minGSSize    = 100,
        maxGSSize    = 500,
        pvalueCutoff = pvalue_cutoff,
        verbose      = FALSE
      )
      
      row <- nrow(gsea_result@result)
      cutoff <- 0.7
      simplifed <- FALSE

      while(row >= max_row){
        simplifed <- TRUE
        gsea_result_simp <- clusterProfiler::simplify(gsea_result, cutoff = cutoff)
        cutoff <- cutoff - 0.1
        row <- nrow(gsea_result_simp@result)
      }

      if(simplifed){
        gsea_result <- gsea_result_simp
      }
      
      # cut_off <- seq(from = 0.1, to = 1, by = 0.1)
      # 
      # gsea_result <- purrr::map(
      #   .x = cut_off,
      #   .f = ~ clusterProfiler::simplify(gsea_result, cutoff = .x)
      # ) %>% 
      #   purrr::keep(~ nrow(.x) <= max_row) %>% 
      #   purrr::pluck(last)
      # 
      gsea_result <- gsea_result@result %>%
        tibble::as_tibble() %>%
        dplyr::mutate(identifier = test)
      
      return(gsea_result)
    },
    gsea_go = function(org_db, ontology = "BP", pvalue_cutoff = 0.05, max_row = 20){
      
      gsea_go_results <-  
        purrr::map(
          .x = self$tested_condition,
          .f = ~ self$gsea_go_single(
            test = .x,
            org_db = org_db,
            ontology = ontology,
            pvalue_cutoff = pvalue_cutoff,
            max_row = max_row
          )
        ) %>%  
        purrr::reduce(dplyr::bind_rows)
      
      return(gsea_go_results)
      
    },
    plot_gsea = function(gsea_result, test){
      
      if(length(test)==1){
        
        is_zero <- gsea_result %>% 
          dplyr::filter(identifier == test) %>%
          nrow()
        
        if(is_zero == 0){
          data <- data.frame(x="No enrichment terms",y="No enrichment terms")
          p <- echarts4r::e_chart(data, x) %>% 
            echarts4r::e_bar(y) %>% 
            echarts4r::e_title(text = test, left = "center") %>% 
            echarts4r::e_legend(show = FALSE)
        }else{
          p <- gsea_result %>% 
            dplyr::filter(identifier == test) %>% 
            dplyr::select(ID, Description, NES, pvalue, p.adjust) %>% 
            dplyr::mutate(p.adjust = -log10(p.adjust)) %>%
            dplyr::mutate(p.adjust = round(p.adjust, 2)) %>% 
            dplyr::mutate(NES = round(NES, 2)) %>% 
            dplyr::mutate(color = dplyr::if_else(NES > 0, "#cf4446", "#0d0887")) %>% 
            dplyr::arrange(p.adjust) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(Description = stringr::str_wrap(Description, width = 50)) %>%
            echarts4r::e_chart(Description, renderer = "svg") %>%
            echarts4r::e_bar(NES, bind = p.adjust) %>%
            echarts4r::e_add_nested("itemStyle", color) %>%
            echarts4r::e_flip_coords() %>%
            echarts4r::e_grid(containLabel = TRUE) %>%
            echarts4r::e_title(text = test, left = "center") %>% 
            echarts4r::e_legend(show = FALSE) %>%
            echarts4r::e_tooltip(
              formatter = htmlwidgets::JS(
                          "
                function(params){
                  return('<strong>' + params.value[1] +
                          '</strong><br />NES: ' + params.value[0] +
                          '<br />-log(p.adj): ' + params.name)
                          }
              "
              )
            ) %>%  
            echarts4r::e_x_axis(
              name = "NES",
              nameLocation = "center",
              nameTextStyle = list(
                fontWeight = "bold",
                fontSize = 16,
                lineHeight = 60
              )
            ) %>% 
            echarts4r::e_toolbox_feature(feature = c("saveAsImage", "dataView"))
        }
        
      }else{
        if(nrow(gsea_result) == 0){
          data <- data.frame(x="No enrichment terms",y="No enrichment terms")
          p <- echarts4r::e_chart(data, x, renderer = "svg") %>% 
            echarts4r::e_bar(y) %>% 
            echarts4r::e_title(text = "No enrichment", left = "center") %>% 
            echarts4r::e_legend(show = FALSE)
        }else{
          p <-
            ggplot(gsea_result, aes(x = identifier, y = stringr::str_wrap(Description, width = 50))) +
            geom_point(aes(color = NES, size = -log(p.adjust)), alpha = 0.9) +
            scale_colour_gradient2(
              low = "#0d0887",
              high = "#cf4446",
              mid = "white",
              midpoint = 0
            ) +
            scale_size(range = c(3, 9)) + # Adjust the range of points size
            theme_minimal() +
            labs(color = "NES", y = NULL, x = NULL) +
            theme(axis.text = element_text(size = 10))
        }
        
      }
      
      
      
      return(p)
      
    },
    go_ora_simplify = function(gene_vector, id, org_db, ontology, alpha, max_row){
      
      ora_result <- clusterProfiler::enrichGO(
        gene          = gene_vector,
        OrgDb         = org_db,
        keyType       = 'SYMBOL',
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.01,
        qvalueCutoff  = 0.05,
        readable      = TRUE
      )
      
      ora_result <- clusterProfiler::filter(ora_result, p.adjust < alpha)
      
      row <- nrow(ora_result@result) 
      cutoff <- 0.7
      simplifed <- FALSE
      
      while(row >= max_row){
        simplifed <- TRUE
        ora_result_simp <- clusterProfiler::simplify(ora_result, cutoff = cutoff)
        cutoff <- cutoff - 0.1
        row <- nrow(ora_result_simp@result)
      }
      
      if(simplifed){
        ora_result <- ora_result_simp
      }
      
      ora_result_tibble <- ora_result@result %>%
        tibble::as_tibble() %>%
        dplyr::mutate(identifier = id)
      
      return(ora_result_tibble)
    },
    ora_go_univariate = function(test, org_db, ontology, alpha, max_row){
      
      table <- self$stat_table %>%
        dplyr::select(gene_names, starts_with(test)) %>%
        rename_at(vars(matches(test)), ~ str_remove(., paste0(test, "_")))
      
      gene_vec_up <- table %>%
        dplyr::filter(significant & fold_change > 0) %>%
        dplyr::pull(gene_names)
      
      gene_vec_down <- table %>%
        dplyr::filter(significant & fold_change < 0) %>%
        dplyr::pull(gene_names)
      
      ora_result_up <-
        self$go_ora_simplify(
          gene_vector = gene_vec_up,
          id = test,
          org_db = org_db,
          ontology = ontology,
          alpha = alpha,
          max_row = max_row
        ) %>% dplyr::mutate(direction = "UP")
      
      ora_result_down <-
        self$go_ora_simplify(
          gene_vector = gene_vec_down,
          id = test,
          org_db = org_db,
          ontology = ontology,
          alpha = alpha,
          max_row = max_row
        ) %>% dplyr::mutate(direction = "DOWN")
      
      
      ora_result <- dplyr::bind_rows(ora_result_up, ora_result_down)
      
      return(ora_result)
      
    },
    ora_go_clusters = function(clusters, org_db, ontology, alpha, max_row){
      
      cluster_vector <- self$cluster_table %>%
        dplyr::filter(cluster == clusters) %>%
        dplyr::pull(gene_names)
      
      ora_result <-
        self$go_ora_simplify(
          gene_vector = cluster_vector,
          id = clusters,
          org_db = org_db,
          ontology = ontology,
          alpha = alpha,
          max_row = max_row
        )
      
      return(ora_result)
      
    },
    go_ora = function(org_db, ontology = "BP", alpha = 0.05, max_row = 20){
      
      if(self$univariate){
        ora_table <- 
          purrr::map(
            .x = self$tested_condition,
            .f = ~ self$ora_go_univariate(
              test = .x,
              org_db = org_db,
              ontology = ontology,
              alpha = alpha,
              max_row = max_row
            )
          ) %>% 
          purrr::reduce(dplyr::bind_rows)
        
      }else if(!self$univariate & !self$clusters_def){
        
        anova_vector <- self$anova_table %>%
          dplyr::filter(significant) %>%
          dplyr::pull(gene_names)
        
        ora_table <- self$go_ora_simplify(
          gene_vector = anova_vector,
          id = "ANOVA significant",
          org_db = org_db,
          ontology = ontology,
          alpha = alpha,
          max_row = max_row
        )
        
      }else if(!self$univariate & self$clusters_def){
        
        ora_table <- 
          purrr::map(
            .x = self$clusters_number,
            .f = ~ self$ora_go_clusters(
              clusters = .x,
              org_db = org_db,
              ontology = ontology,
              alpha = alpha,
              max_row = max_row
            )
          ) %>% 
          purrr::reduce(dplyr::bind_rows)
        
      }
      
      return(ora_table)
      
    },
    plot_ora = function(ora_result, test){
      
      if(length(test)==1){
        
        is_zero <- ora_result %>% 
          dplyr::filter(identifier == test) %>%
          nrow()
        
        if(is_zero == 0){
          data <- data.frame(x="No enrichment terms",y="No enrichment terms")
          p <- echarts4r::e_chart(data, x, renderer = "svg") %>% 
            echarts4r::e_bar(y) %>% 
            echarts4r::e_title(text = test, left = "center") %>% 
            echarts4r::e_legend(show = FALSE)
        }else{
          p <- ora_result %>% 
            dplyr::filter(identifier == test) %>% 
            dplyr::mutate(p.adjust = -log10(p.adjust)) %>%
            dplyr::mutate(p.adjust = round(p.adjust, 2)) %>% 
            {
              if (self$univariate)
                dplyr::mutate(., color = dplyr::if_else(direction == "UP", "#cf4446", "#0d0887"))
              else
                dplyr::mutate(., color = "#440154")
            } %>% 
            dplyr::arrange(p.adjust) %>%
            echarts4r::e_chart(ID, renderer = "svg") %>% 
            echarts4r::e_bar(p.adjust, bind = Description) %>% 
            echarts4r::e_add_nested("itemStyle", color) %>%
            echarts4r::e_flip_coords() %>%
            echarts4r::e_grid(containLabel = TRUE) %>%
            echarts4r::e_legend(show = FALSE) %>% 
            echarts4r::e_tooltip(
              formatter = htmlwidgets::JS(
                          "
                function(params){
                  return('<strong>' + params.name +
                          '</strong><br />-log(p.adj): ' + params.value[0])
                          }
              "
              )
            ) %>%  
            echarts4r::e_x_axis(
              name = "-log(p.adj)",
              nameLocation = "center",
              nameTextStyle = list(
                fontWeight = "bold",
                fontSize = 16,
                lineHeight = 60
              )
            ) %>% 
            echarts4r::e_toolbox_feature(feature = c("saveAsImage", "dataView"))
        }
        
      }else{
        if(nrow(ora_result) == 0){
          data <- data.frame(x="No enrichment terms",y="No enrichment terms")
          p <- echarts4r::e_chart(data, x, renderer = "svg") %>% 
            echarts4r::e_bar(y) %>% 
            echarts4r::e_title(text = "No enrichment", left = "center") %>% 
            echarts4r::e_legend(show = FALSE)
        }else{
          p <-
            ggplot(ora_result, aes(x = identifier, stringr::str_wrap(Description, width = 60))) +
            {
              if (self$univariate)
                geom_point(aes(color = direction, size = -log(p.adjust)))
              else
                geom_point(aes(color = identifier, size = -log(p.adjust)))
            } +
            {
              if (self$univariate)
                scale_color_manual(values = c(DOWN = "#0d0887", UP = "#cf4446")) 
              else
                scale_color_manual(values = head(viridis::viridis(
                  n = length(self$clusters_number) + 1, direction = 1
                ),-1))
            } +
            theme_minimal() +
            labs(size = "-log(p.adjust)", y = NULL, x = NULL) +
            theme(axis.text = element_text(size = 10),
                  legend.text = element_text(size = 10))
        }
      }
      
      return(p)
      
    },
    make_nodes = function(ora_table, test){
      
      nodes <- ora_table %>% 
        {
          if (self$univariate)
            dplyr::filter(., identifier == test)
          else .
        } %>% 
        tidyr::separate_rows(geneID, sep = "/") %>%
        {
          if (self$univariate)
            dplyr::select(.,Description, geneID, direction) %>%
            tidyr::pivot_longer(.,!direction, names_to = "category", values_to = "id")
          else
            dplyr::select(.,Description, geneID, identifier) %>%
            tidyr::pivot_longer(.,!identifier, names_to = "category", values_to = "id")
        } %>% 
        dplyr::distinct(id, .keep_all = TRUE) %>% 
        dplyr::mutate(shape = if_else(category == "geneID", "circle", "diamond")) %>% 
        {
          if (self$univariate)
            dplyr::mutate(.,color = case_when(category == "geneID" & direction == "UP" ~ "#cf444640",
                                              category == "Description" & direction == "UP" ~ "#cf4446",
                                              category == "geneID" & direction == "DOWN" ~ "#0d088740",
                                              TRUE ~ "#0d0887"))
          else .
        } %>%
        dplyr::mutate(label = id) %>% 
        dplyr::mutate(title = id) %>% 
        dplyr::mutate(order = if_else(category == "geneID", 2, 1)) %>% 
        dplyr::arrange(order) %>% 
        {
          if (self$univariate)
            dplyr::select(.,id, group = direction, shape, color, label, title)
          else 
            dplyr::select(.,id, group = identifier, shape, label, title)  
        } %>%
        as.data.frame()
      
      return(nodes)
      
    },
    make_edges = function(ora_table, test){
      
      edges <- ora_table %>% 
        {
          if (self$univariate)
            dplyr::filter(., identifier == test) 
          else .
        } %>% 
        tidyr::separate_rows(geneID, sep = "/") %>% 
        as.data.frame() %>% 
        dplyr::select(from = Description, to = geneID) %>% 
        dplyr::mutate(width = 0.5) 
      
      return(edges)
    } 
  )
)