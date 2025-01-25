# Load required libraries
library(WGCNA)
library(networkD3)
library(igraph)
library(tidyverse)
library(edgebundleR)
library(vroom)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(webshot)
library(ggraph)

# Install PhantomJS for webshot if not already installed
webshot::install_phantomjs()

# Define function to process a single module
process_module <- function(module, evaluation_results, variable_importance_summary, mydata, moduleColors, TOM, deg, num) {
  
  # Create a folder for the current num if it doesn't exist
  output_folder <- paste0("./", num)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
    cat("Created directory for num:", output_folder, "\n")
  }
  
  cat("-------------------------------------------------\n")
  cat("Processing module:", module, "\n")
  
  # Step 1: Extract module-specific genes and TOM
  cat("Step 1: Extracting target genes and constructing TOM for module:", module, "\n")
  target_gene <- variable_importance_summary[variable_importance_summary$Module == module, ]$Gene
  gene <- colnames(mydata)
  inModule <- moduleColors == module
  modgene <- gene[inModule]
  
  if (length(modgene) == 0) {
    cat("No nodes available for module:", module, ". Skipping...\n")
    return()
  }
  
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modgene, modgene)
  filter_modTOM <- modTOM[target_gene, target_gene]
  
  # Step 2: Export Cytoscape files
  cat("Step 2: Exporting Cytoscape node and edge files for module:", module, "\n")
  edge_file <- file.path(output_folder, paste(num, "-step8_CytoscapeInput-edges-", module, ".txt", sep=""))
  node_file <- file.path(output_folder, paste(num, "-step8_CytoscapeInput-nodes-", module, ".txt", sep=""))
  
  exportNetworkToCytoscape(filter_modTOM,                                 
                           edgeFile = edge_file,                                 
                           nodeFile = node_file,                                 
                           weighted = TRUE,                                 
                           threshold = 0,  #weighted权重筛选阈值，可调整 nodeNames = modgene[top],                                  
                           nodeAttr = moduleColors[inModule][target_gene])
  
  # Step 3: Load and process network data
  cat("Step 3: Loading and processing network data for visualization...\n")
  df_node <- vroom::vroom(node_file)
  df_edge <- vroom::vroom(edge_file)
  
  # Check if nodes or edges are empty
  if (nrow(df_node) == 0 || nrow(df_edge) == 0) {
    cat("No valid edges or nodes available for module:", module, ". Skipping...\n")
    return()
  }
  
  df_node$altName <- df_node$nodeName
  df_node[,3] <- module
  
  # Remove self-loops and select top edges
  df_edge <- df_edge %>% arrange(-weight) %>% head(100)
  df_edge <- df_edge[df_edge$fromNode != df_edge$toNode, ]
  
  # Prepare data for visualization
  df_edge_net <- df_edge[, c(1, 2, 3)] %>% as.data.frame()
  df_node_net <- df_node[, c(1, 3)] %>% as.data.frame()
  colnames(df_edge_net) <- c("source", "target", "value")
  colnames(df_node_net) <- c("name", "group")
  
  # Map node IDs
  cat("Mapping node IDs and preparing node/edge data...\n")
  merged_elements <- union_all(df_edge_net$source, df_edge_net$target) %>% unique()
  element_numbers <- seq_along(merged_elements) - 1
  result_df <- data.frame(merged_elements, element_numbers)
  df_edge_net$source <- result_df$element_numbers[match(df_edge_net$source, result_df$merged_elements)]
  df_edge_net$target <- result_df$element_numbers[match(df_edge_net$target, result_df$merged_elements)]
  
  # Add additional node information
  df_node_net <- left_join(result_df, variable_importance_summary[variable_importance_summary$Module == module, ], by = c("merged_elements" = "Gene"))
  df_node_net <- left_join(df_node_net, deg, by = c("merged_elements" = "symbol"))
  df_node_net <- df_node_net[, c("merged_elements", "element_numbers", "Importance", "change")]
  colnames(df_node_net) <- c("name", "group", "size", "type")
  
  # Save processed node data
  cat("Saving processed node file for module:", module, "\n")
  node_output_file <- file.path(output_folder, paste(num, "-df_node_net-", module, ".csv", sep=""))
  write.csv(df_node_net, file = node_output_file, row.names = FALSE)
  
  # Step 4: Generate network plots
  cat("Step 4: Generating network plots for module:", module, "\n")
  generate_network_plots(df_edge_net, df_edge, df_node_net, module, num, output_folder)
  
  # Step 5: Perform GO/KEGG enrichment analysis and save plots
  cat("Step 5: Performing GO/KEGG enrichment analysis for module:", module, "\n")
  generate_GO_KEGG_plots(df_node_net$name, module, num, output_folder)
  
  cat("Module", module, "processed successfully!\n")
}

# Define function to generate network plots
generate_network_plots <- function(df_edge_net, df_edge, df_node_net, module, num, output_folder) {
  cat("  Generating force-directed network plot...\n")
  # Force-directed network
  p1 <- forceNetwork(Links = df_edge_net,
                     Nodes = df_node_net,
                     Source = "source",
                     Target = "target",
                     arrows = TRUE,
                     legend = TRUE,
                     Value = "value",
                     NodeID = "name",
                     Group = "type",
                     colourScale = JS('d3.scaleOrdinal().range(["#eac4d5","#b8e0d4", "#809bce"])'), # Custom color scale
                     bounded = TRUE,
                     opacityNoHover = 0.5,
                     linkDistance = 50,
                     charge = -200,
                     Nodesize = 'size',
                     opacity = 0.9,
                     zoom = TRUE,
                     fontFamily = "Aril",
                     fontSize = 12)
  
  html_output_file <- file.path(output_folder, paste0(num, "-", module, "-network.html"))
  saveNetwork(p1, file = html_output_file)
  
  # Create graph
  cat("  Generating edge-bundled network plot...\n")
  
  # Prepare the igraph object
  g <- graph_from_data_frame(d = df_edge, vertices = df_node_net, directed = TRUE)
  
  # Define a color palette for node groups
  group_colors <- c("up" = "#ffa0c5",   # Red for "up"
                    "down" = "#b8e0d4", # Blue for "down"
                    "stable" = "#88a4c9") # Green for "stable"
  
  # Add group information for nodes
  V(g)$group <- df_node_net$type
  V(g)$color <- group_colors[V(g)$group]
  
  # Add edge weights and gradient colors
  E(g)$weight <- df_edge$weight
  edge_colors <- apply(
    get.edges(g, seq_len(ecount(g))), 1, function(edge) {
      source_group <- V(g)$group[edge[1]]
      target_group <- V(g)$group[edge[2]]
      list(group_colors[source_group], group_colors[target_group])
    }
  )
  E(g)$source_color <- sapply(edge_colors, function(colors) colors[1])
  E(g)$target_color <- sapply(edge_colors, function(colors) colors[2])
  
  # Generate curved edges with gradient colors
  cat("  Generating circular plot with gradient edges...\n")
  df_igraph <- graph_from_data_frame(df_edge, directed = FALSE)
  
  # Create the node attributes (group and size)
  node <- tibble(
    gene = V(df_igraph)$name,
    size = df_node_net$size[match(V(df_igraph)$name, df_node_net$name)],  # Node size from "Importance"
    group = df_node_net$type[match(V(df_igraph)$name, df_node_net$name)]  # Node group based on up/down
  )
  
  # Map node attributes back to the igraph object
  V(df_igraph)$size <- node$size[match(V(df_igraph)$name, node$gene)]
  V(df_igraph)$color <- node$group[match(V(df_igraph)$name, node$gene)]
  
  # Add edge weights (numeric) for gradient colors
  E(df_igraph)$weight <- df_edge$weight  # Assuming 'value' represents edge weights
  
  # Generate the layout
  portraits <- create_layout(df_igraph, layout = 'linear', circular = TRUE)
  
  # Draw the circular graph
  plot <- ggraph(df_igraph, layout = "linear", circular = TRUE) +
    # Draw edges with gradient color based on weight
    geom_edge_arc(aes(color = weight), edge_width = 0.5) +
    scale_edge_color_gradientn(
      colors = c("#809bce","#95b8d1","#b8e0d4","#d6eadf","#f5e2ea","#eac4d5"),  # Gradient colors
      name = "Edge Gradient"
    ) +
    
    # Draw nodes with size and fill color based on group
    geom_node_point(aes(size = size, fill = color), shape = 21) +
    scale_fill_manual(values = group_colors) +
    
    # Add node labels with proper alignment and rotation
    geom_node_text(
      aes(
        label = name, x = x * 1.05, y = y * 1.05,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        vjust = 0.5, hjust = ifelse(x > 0, 0, 1)
      ),
      size = 3
    ) +
    
    # Customize legend and theme
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    theme_graph() + 
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      text = element_text(family = "sans")  # Use a default font to avoid invalid font type error
    ) + 
    coord_fixed(clip = "off")
  
  # Save the circular graph as a PDF
  circular_plot_file <- file.path(output_folder, paste0(num, "-", module, "-circular-graph.pdf"))
  ggsave(plot, filename = circular_plot_file, height = 8, width = 8, device = cairo_pdf)
  cat("  Circular graph saved to:", circular_plot_file, "\n")
  
  # Save edge-bundled network
  output <- edgebundle(g, tension = 0.8)
  html_bundle_file <- file.path(output_folder, paste0(num, "-", module, "-edgebundlenetwork.html"))
  htmlwidgets::saveWidget(widget = output, file = html_bundle_file)
  webshot(html_bundle_file, file.path(output_folder, paste0(num, "-", module, "-edgebundlenetwork.pdf")))
  cat("  Edge-bundled network saved to:", file.path(output_folder, paste0(num, "-", module, "-edgebundlenetwork.pdf")), "\n")
}

# Define function to perform GO/KEGG enrichment and generate plots
generate_GO_KEGG_plots <- function(all_nodes, module, num, output_folder) {
  cat("  Performing GO enrichment analysis...\n")
  # GO enrichment
  enrich_GO <- enrichGO(
    gene = all_nodes,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    keyType = "SYMBOL",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH"
  )
  
  # Save GO results
  go_output_file <- file.path(output_folder, paste0(module, "-", num, "-GO-results.csv"))
  write.csv(as.data.frame(enrich_GO), go_output_file, row.names = FALSE)
  cat("  GO enrichment results saved to:", go_output_file, "\n")
  
  cat("  Performing KEGG enrichment analysis...\n")
  # KEGG enrichment
  entrez_ids <- bitr(all_nodes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  enrich_KEGG <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH"
  )
  enrich_KEGG <- setReadable(enrich_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Save KEGG results
  kegg_output_file <- file.path(output_folder, paste0(module, "-", num, "-KEGG-results.csv"))
  write.csv(as.data.frame(enrich_KEGG), kegg_output_file, row.names = FALSE)
  cat("  KEGG enrichment results saved to:", kegg_output_file, "\n")
  
  cat("  Generating GO/KEGG enrichment plot...\n")
  # Generate GO/KEGG plot
  plot <- GO_KEGG_plot(enrich_GO, enrich_KEGG)
  plot_output_file <- file.path(output_folder, paste0(module, "-GOKEGGplot.pdf"))
  ggsave(plot = plot, filename = plot_output_file, height = 8, width = 8)
  cat("  GO/KEGG enrichment plot saved to:", plot_output_file, "\n")
}

# Main loop to process all modules
module_eva <- unique(evaluation_results$Module)[1:(length(unique(evaluation_results$Module)))]
module_eva

if (inherits(TOM, "dist")) {
  cat("Converting TOM from dist to matrix...\n")
  TOM <- as.matrix(TOM)
}

# Set row and column names if missing
rownames(TOM) <- colnames(TOM) <- colnames(mydata)  # Use gene names from the dataset

cat("Starting module processing...\n")
for (module in module_eva) {
  process_module(
    module = module,
    evaluation_results = evaluation_results,
    variable_importance_summary = variable_importance_summary,
    mydata = mydata,
    moduleColors = moduleColors,
    TOM = TOM,
    deg = deg,
    num = num
  )
}
cat(num, ": All modules processed successfully!\n")