#######################################################################
##  FULL SCRIPT – heat‑map + circle overlay + 5‑level size legend
##  * three tunable parameters at top
##  * automatically finds the two CSV files in the folder
##  * saves a PDF whose name includes thresholds and dataset name
#######################################################################

## ======================= USER PARAMETERS =========================== ##
data_dir       <- "Y:/long/10X_datasets/Xenium/Xenium_5K/t_by_c_result/" %>%
  file.path("t_by_c_Xenium_Prime_Mouse_Pup_FFPE_outs")

t_and_c_thresh <- 0.05      # keep rows with t_and_c < this value
pct_thresh     <- 10        # keep rows with pct      < this value

cell_mm        <- 3         # tile edge length in mm (≈ circle diameter at pct = 100)
base_size      <- 7         # global font size (pt)
## =================================================================== ##


## -------- derive dataset name, input file paths and PDF name -------- ##
folder_name  <- basename(normalizePath(data_dir, winslash = "/"))
dataset_name <- sub("^t_by_c_", "", folder_name)

t_and_c_file <- file.path(data_dir,
                          sprintf("t_and_c_result_%s.csv", dataset_name))
pct_file     <- file.path(data_dir,
                          sprintf("transcript_table_percentage_%s.csv",
                                  dataset_name))

pdf_file <- sprintf("%s_heatmap_t%.3g_pct%.3g.pdf",
                    dataset_name, t_and_c_thresh, pct_thresh)
## -------------------------------------------------------------------- ##


## 1) load required packages ------------------------------------------ ##
required_pkgs <- c("dplyr", "reshape2", "tidyr",
                   "ggplot2", "ggforce", "tibble")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}


## 2) read CSV files --------------------------------------------------- ##
t_and_c <- read.csv(t_and_c_file, row.names = 1, check.names = FALSE)
pct     <- read.csv(pct_file,     row.names = 1, check.names = FALSE)


## 3) melt to long format and apply filters --------------------------- ##
t_long <- melt(as.matrix(t_and_c),
               varnames   = c("Gene", "Cell"),
               value.name = "t_and_c") %>%
  mutate(across(c(Gene, Cell), as.character))

p_long <- melt(as.matrix(pct),
               varnames   = c("Gene", "Cell"),
               value.name = "pct") %>%
  mutate(across(c(Gene, Cell), as.character))

res <- inner_join(t_long, p_long, by = c("Gene", "Cell")) %>%
  filter(t_and_c < t_and_c_thresh, pct < pct_thresh)

genes <- unique(res$Gene)            # genes passing the filter
cells <- colnames(t_and_c)           # keep all cells


## 4) build sub‑matrices --------------------------------------------- ##
mat_sub <- t_and_c[genes, cells, drop = FALSE]
pct_sub <- pct     [genes, cells, drop = FALSE]

mat_sub[mat_sub < 0] <- 0
mat_sub[mat_sub > 1] <- 1
pct_sub[is.na(pct_sub)] <- 0          # replace NAs before plotting


## 5) hierarchical clustering ---------------------------------------- ##
gene_order <- hclust(dist(mat_sub))$labels[hclust(dist(mat_sub))$order]
cell_order <- hclust(dist(t(mat_sub)))$labels[hclust(dist(t(mat_sub)))$order]


## 6) long tables for ggplot + circle parameters ---------------------- ##
df_t <- as.data.frame(mat_sub) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell", values_to = "t_and_c")

df_p <- as.data.frame(pct_sub) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell", values_to = "pct")

plot_df <- inner_join(df_t, df_p, by = c("Gene", "Cell")) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Cell = factor(Cell, levels = cell_order)
  )

circle_df <- plot_df %>%
  mutate(
    x0 = as.numeric(Cell),
    y0 = as.numeric(Gene),
    r  = pct / 100 * 0.5     # pct = 100 → diameter = 1 tile (radius = 0.5)
  )


## -------- legend settings for circle size -------------------------- ##
legend_breaks <- c(5, 25, 45, 65, 85)        # pct values shown in legend
legend_range  <- c(0.1, cell_mm)             # diameter(mm) mapping
## ------------------------------------------------------------------- ##


## 7) build ggplot object -------------------------------------------- ##
p <- ggplot(plot_df, aes(x = Cell, y = Gene)) +
  
  # background tiles
  geom_tile(aes(fill = t_and_c), colour = "white", linewidth = 0.2) +
  
  # circles via radius (data units), no legend
  ggforce::geom_circle(
    data = circle_df,
    aes(x0 = x0, y0 = y0, r = r),
    fill = "white", colour = NA
  ) +
  
  # invisible points for the size legend (alpha=0)
  geom_point(
    aes(size = pct),
    alpha = 0
  ) +
  
  # colour scale (order=2: below size legend), half size
  scale_fill_gradient2(
    low    = "red", mid = "white", high = "blue",
    midpoint = 0.5, limits = c(0, 1), oob = scales::squish,
    name = "SSS",
    guide = guide_colorbar(
      order     = 2,
      barwidth  = unit(cell_mm, "mm"),      # half width
      barheight = unit(cell_mm * 5, "mm")   # half height
    )
  ) +
  
  # size scale for legend keys (5 stacked rows, order=1), half spacing
  scale_size_continuous(
    breaks = legend_breaks,
    range  = legend_range,
    limits = range(legend_breaks),
    name   = "Percentage",
    guide  = guide_legend(
      nrow        = 5,
      byrow       = TRUE,
      order       = 1,
      keyheight   = unit(cell_mm / 2, "mm"),   # half key height
      override.aes = list(
        shape   = 21,
        fill    = "white",
        colour  = "black",
        stroke  = 0.2,
        alpha   = 1
      )
    )
  ) +
  
  coord_fixed(ratio = 1, clip = "off") +
  
  labs(
    title = sprintf("%s (thresholds: SSS < %.3g, Percentage < %.3g)",
                    dataset_name, t_and_c_thresh, pct_thresh),
    x = "Cell", y = "Gene", fill = "t_and_c"
  ) +
  
  theme_minimal(base_size = base_size) +
  theme(
    plot.title    = element_text(size = base_size, hjust = 0.5),
    axis.title.x  = element_text(size = base_size),
    axis.title.y  = element_text(size = base_size),
    axis.text.x   = element_text(size = base_size, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y  = element_text(size = base_size, face = "italic"),
    legend.title  = element_text(size = base_size),
    legend.text   = element_text(size = base_size),
    panel.grid    = element_blank(),
    panel.border  = element_blank()
  )


## 8) compute PDF size and save -------------------------------------- ##
n_r <- nrow(mat_sub);  n_c <- ncol(mat_sub)
width_mm  <- n_c * cell_mm + 20        # small extra margins
height_mm <- n_r * cell_mm + 0

ggsave(pdf_file, plot = p,
       width  = width_mm  / 25.4,        # convert mm to inches
       height = max(2, height_mm / 25.4),# at least 3" tall
       units  = "in", device = cairo_pdf)

cat("Saved:", pdf_file, "\n")
