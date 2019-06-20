# ------------------------------------------------------------------------------
# Plot one-mode adjacency matrix of politicians connected by shared followers
# ------------------------------------------------------------------------------

library(igraph)
library(ggraph)
library(readr)
library(dplyr)

# ------------------------------------------------------------------------------
# colors for party affiliations
# ------------------------------------------------------------------------------

p <- read_csv("data/parties-june2019.csv")
colors = p$color
names(colors) = p$party

# ------------------------------------------------------------------------------
# politicians sample
# ------------------------------------------------------------------------------

d <- cols(
  .default         = col_character(),
  dep17            = col_integer(),
  followers_sample = col_integer(),
  id_str           = col_integer(),
  followers_count  = col_integer(),
  statuses_count   = col_integer(),
  friends_count    = col_integer(),
  verified         = col_logical(),
  last_tweeted     = col_date(format = ""),
  age              = col_integer(), # number of days
  created          = col_integer()  # year
)
d <- read_csv("data/politicians-june2019.csv") %>% 
  rename(
    name2015 = name.x, # 2015, standardized to first_name FAMILY_NAME
    name2019 = name.y,
    description2015 = description.x,
    description2019 = description.y
  )

# ------------------------------------------------------------------------------
# create network object
# ------------------------------------------------------------------------------

f <- "data/network2019.rds"
if (!file.exists(f)) {
  
  y <- readRDS("data/matrix2019.rds")
  
  # one-mode (symmetric) adjacency matrix
  a <- t(y) %*% y
  
  # remove self-loops
  diag(a) <- NA_real_
  
  # sanity check: collapsed on politicians
  stopifnot(colnames(a) %in% d$twitter)
  stopifnot(rownames(a) %in% d$twitter)
  
  # total number of (informative) followers
  s <- colSums(y)
  
  p <- txtProgressBar(min = 1, max = nrow(a), style = 3)
  for (i in 1:nrow(a)) {
    
    setTxtProgressBar(p, i)
    
    for (j in 1:ncol(a)) {
      
      # skip diagonal celles
      if (i == j) next
      
      # make matrix asymmetric:
      #
      # A(i, j) = proportion of followers shared by i and j /
      #           (total followers of i)
      #
      a[ i, j ] <- a[ i, j ] / s[ rownames(a)[ i ] ]

    }
    
  }
  
  saveRDS(a, f)
  
}

a <- readRDS(f)

# ------------------------------------------------------------------------------
# convert to directed one-mode network
# ------------------------------------------------------------------------------

n <- "directed"
n <- graph_from_adjacency_matrix(a, weighted = TRUE, diag = FALSE, mode = n)

# add parties (for node coloring)
x <- d$party
names(x) <- d$twitter
V(n)$party <- x[ V(n)$name ]

# ------------------------------------------------------------------------------
# subset edges and vertices
# ------------------------------------------------------------------------------

# remove ties weighted below 0.5
n <- delete_edges(n, which(E(n)$weight < 0.5))

# keep only main component
n <- delete_vertices(n, which(components(n)$membership > 1))

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------

ggraph(n, layout = "fr") +
  geom_edge_link(color = "grey50", aes(alpha = weight)) +
  geom_node_point(size = 2, alpha = 0.85, aes(color = party)) +
  scale_color_manual("", values = colors, breaks = names(colors)) +
  scale_edge_alpha_continuous(range = c(0.1, 0.5)) +
  guides(edge_alpha = FALSE) +
  ggraph::theme_graph(base_family = "")

ggsave("plots/network2019.pdf", width = 10, height = 9)
ggsave("plots/network2019.png", width = 10, height = 9)

# have a nice day
