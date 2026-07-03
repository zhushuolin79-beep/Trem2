run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}


# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  data <- as.matrix(GetAssayData(visium.slide, assay = assay)) %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}


# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.tmp$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias, cv.folds = 5)
}



seurat.vs <- IR72
gene.expression <- GetAssayData(seurat.vs, assay = "SCT")
coverage <- rowSums(gene.expression > 0) / ncol(gene.expression)
slide.markers <- names(which(coverage >= 0.05))


# estrogen.footprints <- getModel(top = 15) %>%
#   rownames_to_column("gene") %>%
#   filter(Estrogen != 0, gene %in% slide.markers) %>%
#   pull(gene)
# 
# hypoxia.footprints <- getModel(top = 15) %>%
#   rownames_to_column("gene") %>%
#   filter(Hypoxia != 0, gene %in% slide.markers) %>%
#   pull(gene)


view.assays <- list(
  "main" = "giotto",
  "para" = "giotto"
)

features <- rownames(IR72[['giotto']])
# Define features for each view
view.features <- list(
  "main" = features,
  "para" = features
)

# Define spatial context for each view
view.types <- list(
  "main" = "intra",
  "para" = "para"
)

# Define additional parameters (l in the case of paraview)
view.params <- list(
  "main" = NULL,
  "para" = 10
)

misty.out <- "misty_out"

misty.results <- run_misty_seurat(
  visium.slide = seurat.vs,
  view.assays = view.assays,
  view.features = view.features,
  view.types = view.types,
  view.params = view.params,
  spot.ids = NULL, # Using the whole slide
  out.alias = misty.out
) %>% collect_results()



for (region in unique(IR72@meta.data$ST_zone)){
  sub_slide <- subset(IR72, idents = region)
  print (sub_slide)
  
  view.assays <- list(
    "main" = "giotto",
    "para" = "giotto"
  )

  features <- rownames(IR72[['giotto']])
  # Define features for each view
  view.features <- list(
    "main" = features,
    "para" = features
  )

  # Define spatial context for each view
  view.types <- list(
    "main" = "intra",
    "para" = "para"
  )

  # Define additional parameters (l in the case of paraview)
  view.params <- list(
    "main" = NULL,
    "para" = 10
  )

  misty.out <- paste0(region, '__misty_out')

  misty.results <- run_misty_seurat(
    visium.slide = sub_slide,
    view.assays = view.assays,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    spot.ids = NULL, # Using the whole slide
    out.alias = misty.out
  ) 

}

for (region in unique(IR72@meta.data$ST_zone)){
  misty.results <- collect_results(paste0(region, '__misty_out'))
  saveRDS(misty.results, paste0(region, '_misty_res.rds'))
}

