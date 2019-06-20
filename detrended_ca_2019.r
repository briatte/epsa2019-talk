#===============================================================================
#
# A3_detrended_ca.r -- detrended correspondence analysis of the model matrix
#
# The script runs detrended correspondence analysis (CA) on the same matrix used
# in the rest of the analysis. The detrending suppresses the arch effect created
# by FN politicians and followers, returns at least two sensible dimensions, and
# correlates to the estimated ideal points obtained through the Bayesian model.
#
# See Hill, M.O. and H.G. Gauch, Jr. 1980. Detrended Correspondence Analysis:
# an improved ordination technique. Vegetatio 42:47-58.
#
#===============================================================================

library(dplyr)
library(ggplot2)
library(readr)

# correspondence analysis
library(ca)
stopifnot(packageVersion("ca") >= 0.61)

# detrended correspondence analysis
library(vegan)

# ------------------------------------------------------------------------------
# colors for party affiliations
# ------------------------------------------------------------------------------

p <- read_csv("data/parties-june2019.csv")
colors = p$color
names(colors) = p$party

# ------------------------------------------------------------------------------
# plot theme
# ------------------------------------------------------------------------------

theme_paper <- theme_bw(14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = rel(1)),
    axis.title = element_text(size = rel(1)),
    strip.text = element_text(size = rel(1)),
    strip.background = element_rect(fill = "grey90"),
    legend.text = element_text(size = rel(1)),
    axis.ticks.x = element_blank()
  )

# ------------------------------------------------------------------------------
# load politicians sample
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

cat("[TOFIX] Removing", sum(is.na(d$party)), "politicians with no `party`\n")
d <- filter(d, !is.na(party))

# ------------------------------------------------------------------------------
# load follower-followee incidence matrix
# ------------------------------------------------------------------------------

y <- readRDS("data/matrix2019.rds")

cat("Selected matrix:", nrow(y), "rows,", ncol(y), "cols\n")

min(colSums(y)) # politicians: min. number of informative followers (200)
min(rowSums(y)) # users: min. number of politicians followed (4)

# ------------------------------------------------------------------------------
# correspondence analysis + detrended CA
# ------------------------------------------------------------------------------

# works only with ca >= 0.61
ca_2d <- ca(y, nd = 2)

r <- left_join(
  tibble(
    twitter = ca_2d$colnames,
    d1 = ca_2d$colcoord[, 1],
    d2 = ca_2d$colcoord[, 2]
  ),
  select(d, twitter, party),
  by = "twitter"
)

saveRDS(r, "models/ca2019.rds")

# problem -- FN + DVD (DLF, PCD) politicians are orthogonal to all others
ggplot(r, aes(d1, d2, color = party)) +
  geom_point() +
  scale_color_manual("", breaks = names(colors), values = colors) +
  labs(y = "Dimension 2\n", x = "\nDimension 1") +
  theme_paper +
  theme(legend.key = element_blank())

ggsave("plots/ca_2d_2019.png", width = 7, height = 5)
ggsave("plots/ca_2d_2019.pdf", width = 7, height = 5)

# solution -- switch to detrended correspondence analysis
ord <- decorana(
  t(y),
  iweigh = 0, # no downweighting of rare groups
  iresc = 4,  # number of rescaling cycles
  ira = 0,    # detrended
  mk = 26,    # number of segments
  short = 0   # rescale shortest gradient
)

res <- left_join(
  # first two dimensions of detrended CA
  tibble(
    twitter = rownames(scores(ord)),
    DCA1 = scores(ord)[, 1],
    DCA2 = scores(ord)[, 2],
    DCA3 = scores(ord)[, 3],
    DCA4 = scores(ord)[, 4]
  ),
  # politicians' party and Twitter account attributes
  select(d, twitter, party, followers_count, statuses_count, age),
  by = "twitter"
) %>%
  arrange(DCA1)

cat("[TOFIX] Removing", sum(is.na(res$party)), "DCA results with no `party`\n")
res <- filter(res, !is.na(party))

group_by(res, party) %>%
  tally %>% 
  print

saveRDS(res, "models/dca2019.rds")

# ------------------------------------------------------------------------------
# correlations
# ------------------------------------------------------------------------------

# correlations of covariates to each dimension

cat("Dim. 1 correlations:\n")

cat("- Popularity (followers): ")
with(res, cor(DCA1, log10(followers_count), use = "pairwise.complete.obs")) %>% 
  cat("\n") # ~ 0
cat("- Activity (statuses): ")
with(res, cor(DCA1, log10(statuses_count), use = "pairwise.complete.obs")) %>% 
  cat("\n") # ~ 0
cat("- Age (opening date): ")
with(res, cor(DCA1, age, use = "pairwise.complete.obs")) %>% 
  cat("\n") # -.24, weak negative

cat("Dim. 2 correlations:\n")

cat("- Popularity (followers): ")
with(res, cor(DCA2, log10(followers_count), use = "pairwise.complete.obs")) %>% 
  cat("\n") # .36, moderate
cat("- Activity (statuses): ")
with(res, cor(DCA2, log10(statuses_count), use = "pairwise.complete.obs")) %>% 
  cat("\n") # .28, weak
cat("- Age (opening date): ")
with(res, cor(DCA2, age, use = "pairwise.complete.obs")) %>% 
  cat("\n") # .11, v. weak

# ------------------------------------------------------------------------------
# plots of dimensions 1 and 2
# ------------------------------------------------------------------------------

ggplot(res, aes(DCA1, DCA2, color = party)) +
  geom_point() +
  scale_color_manual("", breaks = names(colors), values = colors) +
  theme_paper

ggsave("plots/dca_2d_2019.png", width = 7, height = 5)
ggsave("plots/dca_2d_2019.pdf", width = 7, height = 5)

ggplot(res, aes(DCA1, DCA2), color = "grey75") +
  geom_point(data = select(res, -party), color = "grey50", alpha = .25) +
  geom_point(data = mutate(res, party = "(all)"), color = "grey50", alpha = .25) +
  geom_point() +
  geom_vline(xintercept = mean(res$DCA1), lty = "dotted") +
  geom_hline(yintercept = mean(res$DCA2), lty = "dotted") +
  facet_wrap(~ party) +
  theme_paper

ggsave("plots/dca_2d_by_party_2019.png", width = 10, height = 9)
ggsave("plots/dca_2d_by_party_2019.pdf", width = 10, height = 9)

# ------------------------------------------------------------------------------
# boxplots of dimension 1 for selected parties
# ------------------------------------------------------------------------------

# recode smaller parties
g <- res %>% 
  # absorb small parties
  mutate(party = recode(party, "DVD-PCD" = "DVD", "DVG-PRG" = "DVG")) %>%
  # remove independents
  filter(party != "IND")

ggplot(g, aes(reorder(party, DCA1, median), DCA1, color = party)) +
  geom_boxplot() +
  scale_color_manual("", breaks = names(colors), values = colors) +
  guides(color = FALSE) +
  theme_paper +
  labs(x = NULL, y = "First dimension of detrended CA\n")

ggsave("plots/dca_parties_2019.png", width = 10, height = 5)
ggsave("plots/dca_parties_2019.pdf", width = 10, height = 5)

# have a nice day
