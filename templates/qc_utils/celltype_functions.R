# Create data frame of cell types


#' Create `celltype_df` data frame for use in cell type QC reports
#'
#' @param processed_sce The processed sce object
#' @param has_singler Boolean for whether SingleR annotations are present
#' @param has_cellassign Boolean for whether CellAssign annotations are present
#'
#' @return `celltype_df` with column of cell types, as factors, for each annotation method
create_celltype_df <- function(processed_sce, has_singler, has_cellassign) {
  celltype_df <- colData(processed_sce) |>
    as.data.frame() |>
    # barcodes to a column
    tibble::rownames_to_column(var = "barcode") |>
    # keep only cell name, celltyping, and clusters
    dplyr::select(
      barcode,
      clusters,
      contains("singler"),
      contains("cellassign"),
      contains("submitter")
    )

  if (has_singler) {
    celltype_df <- celltype_df |>
      prepare_annotation_values(singler_celltype_annotation)
  }
  if (has_cellassign) {
    celltype_df <- celltype_df |>
      prepare_annotation_values(cellassign_celltype_annotation)
  }

  return (celltype_df)
}

# Functions for working with cell type annotations in QC reports

#' Prepare and reformat cell type annotation values for use in QC reports
#'  Unknown cell types are updated with the label "Unknown cell type", and
#'  cell types are ordered in order of descending frequency, but with
#'  "Unknown cell type" as the last level
#'
#' @param df The data frame containing cell type annotations, one row per cell
#' @param annotation_column The column (written plainly, not a string) containing annotations to reformat
#'
#' @return Updated data frame with the `annotation_column` reformatted
prepare_annotation_values <- function(df, annotation_column) {
  df |>
    dplyr::mutate(
      {{ annotation_column }} := dplyr::case_when(
        # singler condition
        is.na({{ annotation_column }}) ~ "Unknown cell type",
        # cellassign conditon
        {{ annotation_column }} == "other" ~ "Unknown cell type",
        # otherwise, keep it
        TRUE ~ {{ annotation_column }}
      ) |>
        forcats::fct_infreq() |>
        forcats::fct_relevel("Unknown cell type", after = Inf)
    )
}



#' Create tables of cell type annotation counts
#'
#' @param df Data frame with cell types
#' @param celltype_column Column with cell type annotations, not a string.
#'
#' @return kableExtra table with cell type counts.
create_celltype_n_table <- function(df, celltype_column) {
  df |>
    dplyr::count({{ celltype_column }}) |>
    # Add percentage column
    dplyr::mutate(
      `Percent of cells` = paste0(round(n / sum(n) * 100, digits = 2), "%")
    ) |>
    # set column order & rename
    dplyr::select(
      `Annotated cell type` = {{ celltype_column }},
      `Number of cells` = n,
      `Percent of cells`
    ) |>
    # kable formatting
    knitr::kable(align = "r") |>
    kableExtra::kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    ) |>
    kableExtra::column_spec(2, monospace = TRUE)
}


# Function to add a column of lumped cell types in a data frame
#'
#' @param df Data frame to manipulate
#' @param celltype_df Data frame that contains column of cell types,
#'   named `{celltype_method}_celltype_annotation`
#' @param celltype_method Cell type method
#' @param n_celltypes Number of groups to lump into, with rest put into "Other" group. Default is 7.
#'
#' @return Updated df with new column of lumped celltypes called `{celltype_method}`
lump_celltypes <- function(df,
                           celltype_df,
                           celltype_method,
                           n_celltypes = 7) {
  df |>
    dplyr::mutate(
      celltypes = celltype_df[[glue::glue("{celltype_method}_celltype_annotation")]] |>
        forcats::fct_lump_n(
          n_celltypes,
          other_level = "Other cell type"
        )
    ) |>
    # finally, rename temporary `celltypes` column to the provided method
    dplyr::rename({{ celltype_method }} := celltypes)
}




#' Make UMAP  colored by given variable
#'
#' @param umap_df Data frame with UMAP1 and UMAP2 columns
#' @param color_variable Column in data frame to color by, not a string.
#' @param legend_title Title for legend.
#' @param legend_nrow Number of rows in legend. Default is 2.
#'
#' @return UMAP plot as a ggplot2 object
plot_umap <- function(umap_df,
                      color_variable,
                      legend_title,
                      legend_nrow = 2) {
  ggplot(umap_df) +
    aes(
      x = UMAP1,
      y = UMAP2,
      color = {{ color_variable }}
    ) +
    geom_point(
      size = 0.3,
      alpha = 0.5
    ) +
    # remove axis numbers and background grid
    scale_x_continuous(labels = NULL, breaks = NULL) +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    coord_fixed() +
    guides(
      color = guide_legend(
        title = legend_title,
        nrow = legend_nrow,
        # more visible points in legend
        override.aes = list(
          alpha = 1,
          size = 1.5
        )
      )
    ) +
    theme(legend.position = "bottom")
}



#' Create a heatmap of cell type annotations with log1p transformation
#'
#' @param x_vector Vector of values for the x-axis (rows)
#' @param y_vector Vector of values for the y-axis (columns)
#' @param x_label x-axis label. Default is no label.
#' @param y_label y-axis label. Default is no label.
#' @param y_title_location location of the y-axis title. Default is on the bottom.
#' @param column_names_rotation degree to rotate column names. Default is 0.
#' @param row_font_size Size of row font. Default is 8.
#' @param column_font_size. Size of column font. Default is 10.
#'
#' @return
#' @export
#'
#' @examples
create_celltype_heatmap <- function(x_vector,
                                    y_vector,
                                    x_label = "",
                                    y_label = "",
                                    y_title_location = "bottom",
                                    column_names_rotation = 0,
                                    row_font_size = 8,
                                    column_font_size = 10) {
  # build a matrix for making a heatmap
  celltype_mtx <- table(
    x_vector,
    y_vector
  ) |>
    log1p() # log transform for visualization

  # Define CVD-friendly palette
  heatmap_palette <- viridisLite::inferno(7, alpha = 1, begin = 0, end = 1, direction = 1)

  # heatmap
  heat <- ComplexHeatmap::Heatmap(celltype_mtx,
                                  # Overall heatmap parameters
                                  col = heatmap_palette,
                                  # Column parameters
                                  column_title = y_label,
                                  column_title_side = y_title_location,
                                  column_dend_side = "top",
                                  column_names_rot = column_names_rotation,
                                  column_names_gp = grid::gpar(fontsize = column_font_size),
                                  # Row parameters
                                  row_dend_side = "left",
                                  row_title_side = "right",
                                  row_title = x_label,
                                  row_names_gp = grid::gpar(fontsize = row_font_size),
                                  # Legend parameters
                                  heatmap_legend_param = list(
                                    title = "Log(Number of cells)",
                                    title_position = "leftcenter-rot",
                                    legend_height = unit(4, "cm")
                                  )
  )
  # draw with legend on left for spacing
  ComplexHeatmap::draw(heat, heatmap_legend_side = "left")
}