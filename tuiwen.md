# Day 7: Data Visualization with ggplot2

## 1. Introduction

In the previous sections, we introduced R data structures, data manipulation (tidyverse), control flow, and function writing. Now let's learn how to create elegant and professional data visualizations using `ggplot2`.


![](https://fastly.jsdelivr.net/gh/bucketio/img17@main/2026/03/21/1774102978772-dab98b3e-7165-448b-a89a-2fc902a42e41.png)


`ggplot2` is a core package of the tidyverse ecosystem, based on Leland Wilkinson's *Grammar of Graphics*. Its core idea is: **Deconstruct graphs into independent layers**, and build complex visualizations by combining them.

<br>

## 2. Grammar of Graphics in ggplot2

### 2.1 Core Concepts

`ggplot2` decomposes a plot into the following components:

1. **Data**: The data frame to visualize
2. **Aesthetics (aes)**: Mapping of data variables to visual properties (x, y, color, size, etc.)
3. **Geometries (Geoms)**: The geometric representation of data (points, lines, bars, etc.)
4. **Statistics (Stats)**: Statistical transformations of the data (counting, smoothing, etc.)
5. **Scales**: Controls the details of aesthetic mappings (color schemes, axis limits, etc.)
6. **Coordinate systems**: How data points are mapped to the 2D plane
7. **Facets**: Grouping data to create multiple subplots
8. **Themes**: Controls the overall appearance of the plot


![](https://fastly.jsdelivr.net/gh/bucketio/img11@main/2026/03/21/1774103000285-19909b68-04bb-4ba2-b293-243267547ffe.png)


<br>

### 2.2 Basic Syntax Structure

```r
ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +
  <GEOM_FUNCTION>() +
  <SCALE_FUNCTION>() +
  <FACET_FUNCTION>() +
  <THEME_FUNCTION>()
```

**Key Points**:
- Use the `+` sign to chain layers together (not the pipe operator `%>%` or `|>`)
- `ggplot()` initializes the plot object and specifies global data/aesthetics
- Each `geom_*()` function adds a new geometric layer

<br>

## 3. Your First ggplot2 Plot

### 3.1 Installation and Loading

```r
# Install (if not already installed)
# install.packages("ggplot2")

library(ggplot2)
library(dplyr)  
```

<br>

### 3.2 Using Built-in Datasets - Scatter Plot

We use the built-in `mpg` dataset (car fuel efficiency data) as an example.

```r
# View data
head(mpg)

# Basic scatter plot: Engine displacement vs Highway MPG
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point()
```

**Explanation**:
- `ggplot()`: Initializes the plot
- `aes(x = displ, y = hwy)`: Maps `displ` to the x-axis and `hwy` to the y-axis
- `geom_point()`: Represents data as points

<br>

### 3.3 Adding Color Mapping

```r
# Color by vehicle class
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point()
```

Colors are automatically grouped by the `class` variable, generating a legend.

<br>

### 3.4 Adding Size and Transparency

```r
# Map dot size to number of cylinders, add transparency
ggplot(mpg, aes(x = displ, y = hwy, color = class, size = cyl)) +
  geom_point(alpha = 0.6)
```

**Note**:
- `alpha` controls transparency (0-1). It can be mapped in `aes()` or passed as a fixed parameter.
- When `size` is mapped to a numeric variable, the point size will scale with the values.

<br>

## 4. Common Geometries (Geoms)

### 4.1 Scatter Plot (geom_point)

Ideal for displaying the relationship between two continuous variables.

```r
set.seed(123)
gene_expr <- data.frame(
  gene = paste0("Gene", 1:100),
  sample_A = rnorm(100, mean = 5, sd = 2),
  sample_B = rnorm(100, mean = 5, sd = 2),
  significant = sample(c("Sig", "Not Sig"), 100, replace = TRUE, prob = c(0.3, 0.7))
)

ggplot(gene_expr, aes(x = sample_A, y = sample_B, color = significant)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Gene Expression Correlation",
       x = "Sample A (log2 expression)",
       y = "Sample B (log2 expression)") +
  theme_bw()
```

<br>

### 4.2 Line Plot (geom_line)

Ideal for time series or ordered data.

```r
# Example: qPCR timeseries
time_series <- data.frame(
  time = rep(0:10, 3),
  expression = c(1 * 2^(0:10 * 0.3),  # Gene1
                 1 * 2^(0:10 * 0.5),  # Gene2
                 1 * 2^(0:10 * 0.2)), # Gene3
  gene = rep(c("Gene1", "Gene2", "Gene3"), each = 11)
)

ggplot(time_series, aes(x = time, y = expression, color = gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +  # Log scale
  labs(title = "Gene Expression Time Series",
       x = "Time (hours)",
       y = "Expression (log scale)") +
  theme_minimal()
```

<br>

### 4.3 Bar Plot (geom_bar / geom_col)

- `geom_bar()`: Automatically counts occurrences, used for categorical variables
- `geom_col()`: Uses raw y values directly, used for pre-summarized data

```r
# geom_bar: Automatic counting
ggplot(mpg, aes(x = class)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Vehicle Class Distribution", x = "Class", y = "Count") +
  theme_classic()
```

```r
# geom_col: Using pre-calculated values
pathway_counts <- data.frame(
  pathway = c("MAPK", "PI3K-AKT", "WNT", "Notch", "TGF-beta"),
  gene_count = c(45, 38, 25, 18, 22)
)

ggplot(pathway_counts, aes(x = reorder(pathway, gene_count), y = gene_count)) +
  geom_col(fill = "coral") +
  coord_flip() +  # Horizontal bar plot
  labs(title = "Signal Pathway Gene Counts", x = NULL, y = "Gene Count") +
  theme_minimal()
```

<br>

### 4.4 Boxplot (geom_boxplot)

Displays data distribution, median, quartiles, and outliers.

```r
# Example: Gene expression across tissues
tissue_expr <- data.frame(
  tissue = rep(c("Liver", "Brain", "Muscle", "Heart"), each = 50),
  expression = c(rnorm(50, 8, 2), rnorm(50, 6, 1.5), 
                 rnorm(50, 7, 2.5), rnorm(50, 9, 1.8))
)

ggplot(tissue_expr, aes(x = tissue, y = expression, fill = tissue)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +  # Add raw data points
  labs(title = "Gene Expression by Tissue",
       x = "Tissue Type", y = "Expression (log2)") +
  theme_bw() +
  theme(legend.position = "none")
```

<br>

### 4.5 Violin Plot (geom_violin)

Combines advantages of boxplots and density plots.

```r
ggplot(tissue_expr, aes(x = tissue, y = expression, fill = tissue)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(title = "Gene Expression Distribution (Violin Plot)",
       x = "Tissue Type", y = "Expression (log2)") +
  theme_minimal() +
  theme(legend.position = "none")
```

<br>

### 4.6 Histogram (geom_histogram)

Displays the distribution of a single continuous variable.

```r
ggplot(gene_expr, aes(x = sample_A)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(gene_expr$sample_A), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Gene Expression Distribution",
       x = "Expression (log2)", y = "Frequency") +
  theme_classic()
```

<br>

### 4.7 Density Plot (geom_density)

A smoothed version of the histogram.

```r
ggplot(tissue_expr, aes(x = expression, fill = tissue)) +
  geom_density(alpha = 0.5) +
  labs(title = "Gene Expression Density by Tissue",
       x = "Expression (log2)", y = "Density") +
  theme_minimal()
```

<br>

## 5. Aesthetics Mappings

### 5.1 Global vs Local Mapping

```r
# Global mapping: Shared across all layers
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)  # Inherits color mapping

# Local mapping: Applied to a specific layer
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point(aes(color = class)) +  # Only points have colors
  geom_smooth(method = "lm", se = FALSE, color = "black")  # Fixed black line
```

<br>

### 5.2 Fixed vs Mapped Properties

```r
# Incorrect example: Fixed value inside aes()
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point(aes(color = "blue"))  # Error! Will be treated as a variable mapping

# Correct example: Fixed value outside aes()
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point(color = "blue", size = 3)
```

<br>

### 5.3 Common Aesthetic Properties

| Property | Description | Applicable Geom |
|------|------|-----------|
| `x`, `y` | Coordinate axis position | All |
| `color` | Point/line color | point, line, text |
| `fill` | Fill color | bar, boxplot, violin, area |
| `size` | Point/line size | point, line, text |
| `alpha` | Transparency (0-1) | All |
| `shape` | Point shape | point |
| `linetype` | Line pattern | line, smooth |

<br>

## 6. Facets

Faceting creates multiple subplots based on a categorical variable.

### 6.1 facet_wrap()

Wraps subplots based on a single variable into multiple rows/cols.

```r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_wrap(~ class, nrow = 2) +
  labs(title = "Faceted by Vehicle Class") +
  theme_bw()
```

<br>

### 6.2 facet_grid()

Creates a grid based on two categorical variables.

```r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_grid(drv ~ cyl) +
  labs(title = "Faceted by Drive Type and Cylinders") +
  theme_minimal()
```

<br>

## 7. Scales

Scales define how data values are mapped to visual properties.

### 7.1 Axis Scales

```r
# Log scale
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(10, 50, by = 5)) +
  labs(title = "Log10 x-axis") +
  theme_bw()
```

<br>

### 7.2 Color Scales

```r
# Manually assigned colors
ggplot(tissue_expr, aes(x = tissue, y = expression, fill = tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Liver" = "#E64B35", 
                                "Brain" = "#4DBBD5",
                                "Muscle" = "#00A087",
                                "Heart" = "#F39B7F")) +
  theme_minimal()

# ColorBrewer palettes
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_classic()
```

<br>

### 7.3 Gradient Color Scales

```r
# Continuous variable color gradients
gene_heatmap_data <- expand.grid(
  gene = paste0("Gene", 1:20),
  sample = paste0("Sample", 1:10)
)
gene_heatmap_data$expression <- rnorm(200, mean = 5, sd = 2)

ggplot(gene_heatmap_data, aes(x = sample, y = gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 5) +
  labs(title = "Gene Expression Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<br>

## 8. Themes

### 8.1 Built-in Themes

```r
p <- ggplot(mpg, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Different Themes Example")

# theme_gray() (default)
p + theme_gray()

# theme_bw() (black and white, good for papers)
p + theme_bw()

# theme_minimal() (minimalist)
p + theme_minimal()

# theme_classic() (classic plot with no gridlines)
p + theme_classic()
```

<br>

### 8.2 Custom Theme Elements

```r
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 3) +
  labs(title = "Custom Theme Example",
       x = "Engine Displacement (L)",
       y = "Highway MPG") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  )
```

<br>

## 9. Bioinformatics Real-world Examples

### 9.1 Volcano Plot

```r
# Simulate differentially expressed gene data
set.seed(456)
n_genes <- 5000
volcano_data <- data.frame(
  gene = paste0("Gene", 1:n_genes),
  log2FC = rnorm(n_genes, mean = 0, sd = 1.5),
  pvalue = runif(n_genes, 0, 1)
)
volcano_data$padj <- p.adjust(volcano_data$pvalue, method = "BH")
volcano_data$neg_log10_padj <- -log10(volcano_data$padj)

# Categorization: Significant Up/Down/Not Sig
volcano_data <- volcano_data %>%
  mutate(
    diff_expressed = case_when(
      log2FC > 1 & padj < 0.05 ~ "Up-regulated",
      log2FC < -1 & padj < 0.05 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Draw Volcano Plot
ggplot(volcano_data, aes(x = log2FC, y = neg_log10_padj, color = diff_expressed)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red",
                                 "Down-regulated" = "blue",
                                 "Not significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  labs(title = "Volcano Plot - Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Status") +
  theme_bw() +
  theme(legend.position = "top")
```

<br>

### 9.2 MA Plot

```r
# Simulated Data
ma_data <- data.frame(
  gene = paste0("Gene", 1:n_genes),
  baseMean = 10^runif(n_genes, 0, 4),
  log2FC = rnorm(n_genes, 0, 1.2)
)
ma_data$significant <- abs(ma_data$log2FC) > 1

ggplot(ma_data, aes(x = log10(baseMean), y = log2FC, color = significant)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(title = "MA Plot",
       x = "Log10 Mean Expression",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none")
```

<br>

### 9.3 PCA Plot (Principal Component Analysis)

```r
# Simulated PCA Results
set.seed(789)
pca_data <- data.frame(
  sample = paste0("S", 1:30),
  PC1 = c(rnorm(15, -2, 1), rnorm(15, 2, 1)),
  PC2 = c(rnorm(15, 1, 1.5), rnorm(15, -1, 1.5)),
  group = rep(c("Control", "Treatment"), each = 15)
)

ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +  # Add confidence ellipses
  labs(title = "PCA Plot",
       x = "PC1 (45.3% variance)",
       y = "PC2 (23.7% variance)") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank())
```

<br>

## 10. Combining Plots

### 10.1 Using patchwork

```r
# install.packages("patchwork")
library(patchwork)

p1 <- ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  labs(title = "Plot 1") +
  theme_bw()

p2 <- ggplot(mpg, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Plot 2") +
  theme_bw() +
  theme(legend.position = "none")

p3 <- ggplot(mpg, aes(x = hwy)) +
  geom_histogram(bins = 20, fill = "steelblue") +
  labs(title = "Plot 3") +
  theme_bw()

# Combine layout
(p1 | p2) / p3 + 
  plot_annotation(title = "Combined Plots Layout Example", 
                  tag_levels = "A")
```

<br>

## 11. Saving Plots

```r
# Save current plot
ggsave("my_plot.png", width = 8, height = 6, dpi = 300)
ggsave("my_plot.pdf", width = 8, height = 6)

# Save specific plot object
p <- ggplot(mpg, aes(x = displ, y = hwy)) + geom_point()
ggsave("scatter.png", plot = p, width = 10, height = 8, dpi = 300)
```

<br>

## 12. Best Practices & Tips

### 12.1 Data Preparation

```r
# Prepare data using tidyverse
library(tidyr)

# Wide format vs Long format
wide_data <- data.frame(
  sample = c("S1", "S2", "S3"),
  GeneA = c(5.2, 6.1, 5.8),
  GeneB = c(7.3, 6.9, 7.5)
)

# Convert to long format for ggplot2
long_data <- wide_data %>%
  pivot_longer(cols = c(GeneA, GeneB), 
               names_to = "gene", 
               values_to = "expression")

ggplot(long_data, aes(x = sample, y = expression, fill = gene)) +
  geom_col(position = "dodge") +
  theme_minimal()
```

<br>

### 12.2 Code Organization

```r
# Formatting ggplot code across lines for readability 
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Engine Displacement vs Highway MPG",
    subtitle = "Colored by vehicle class",
    x = "Engine Displacement (L)",
    y = "Highway Miles per Gallon",
    color = "Vehicle Class"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )
```

<br>

### 12.3 Frequently Asked Questions

**Issue 1**: Rendering CJK Fonts
```r
# macOS / Linux
theme(text = element_text(family = "STHeiti"))

# Windows
theme(text = element_text(family = "SimHei"))
```

**Issue 2**: Overlapping axis labels
```r
ggplot(mpg, aes(x = manufacturer, y = hwy)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

**Issue 3**: Legend Position Adjustment
```r
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point() +
  theme(legend.position = "bottom")  # "top", "left", "right", "none"
```

<br>

## 13. Recommended Extension Packages

- **ggpubr**: Publication-ready plots (with statistical testing)
- **ggsci**: Color palettes matching scientific journals
- **ggrepel**: Automatically repels text labels away from data points
- **patchwork**: Combining multiple plots
- **plotly**: Interactive graphics (`ggplotly()`)
- **gganimate**: Animated plots

<br>

## 14. Summary

This section covered the core concepts and usage of ggplot2:

1. **Grammar of Graphics**: Data + Aesthetics + Geometries + ...
2. **Common Geoms**: Points, lines, bars, boxplots, violins, heatmaps
3. **Aesthetics Mappings**: Global vs Local, Fixed vs Mapped
4. **Facets**: `facet_wrap()` and `facet_grid()`
5. **Scales**: Axis, color palettes, gradients
6. **Themes**: Built-in themes and custom themes
7. **Bioinformatics Cases**: Volcano plots, MA plots, PCA plots

Next Steps:
- Read the [ggplot2 Official Documentation](https://ggplot2.tidyverse.org/)
- Browse the [R Graph Gallery](https://r-graph-gallery.com/) for inspiration
- Practice by visualizing real datasets
- Explore the ggplot2 extension ecosystem

---

**References**:
- [ggplot2: Elegant Graphics for Data Analysis](https://ggplot2-book.org/)
- [R for Data Science - Data Visualization](https://r4ds.hadley.nz/data-visualize)
- [ggplot2 Cheatsheet](https://github.com/rstudio/cheatsheets/blob/master/data-visualization.pdf)
