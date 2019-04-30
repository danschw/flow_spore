source("src/functions.R")

#### Figure 1 - Separation of Cells vs Spores possible, overview plot ####
figure_1 <- function() {
  # sample variables
  sample.var <- c("strain", "ident", "stain", "mode", "run")

  df1 <- flowCreateFlowSet(filepath = "data/f1/", sample_variables = sample.var) %>%
    Subset(., norm2Filter("asinh.FSC.H", "asinh.FSC.W", filterId = "norm_ssc.fsc", scale = 1)) %>%
    Subset(., norm2Filter("asinh.SSC.H", "asinh.SSC.W", filterId = "norm_ssc.fsc", scale = 1)) %>%
    flowFcsToDf(.)

  SSC.plot <- df1 %>%
    dplyr::filter(mode == "SSC") %>%
    select(asinh.SSC.A, ident, stain) %>% # ,run
    # dplyr::filter(stain=="unstained",ident %in% c("spores","cells"))%>%
    ggplot(aes(asinh.SSC.A)) +
    geom_density(aes(fill = ident, y = ..scaled..), alpha = 0.4) +
    ylab("norm. count") +
    guides(fill = FALSE) + theme(legend.position = c(.05, .9)) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8, name = "", direction = -1)

  FSC.plot <- df1 %>%
    dplyr::filter(mode == "FSC") %>%
    select(asinh.FSC.A, ident, stain, run) %>%
    ggplot(aes(asinh.FSC.A)) +
    geom_density(aes(fill = ident, y = ..scaled..), alpha = 0.4) +
    ylab("") +
    guides(fill = FALSE) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8, direction = -1)

  PI.plot <- df1 %>%
    dplyr::filter(mode == "PI") %>%
    select(asinh.FL3.A, ident, stain, run) %>%
    # dplyr::filter(stain=="PI",run=="2x")%>%
    ggplot(aes(asinh.FL3.A)) +
    geom_density(aes(fill = ident, y = ..scaled..), alpha = 0.4) +
    guides(fill = FALSE) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8, direction = -1) +
    xlim(c(5, 12)) +
    # annotate(geom = "segment", x = 9, xend = 8.6, y = 0.3, yend = 0.1,
    #          arrow=arrow(),colour = "red",alpha=0.7,size=1)+
    ylab("")

  PI.plot

  S1.plot <- df1 %>%
    dplyr::filter(mode == "SYBR1") %>%
    select(asinh.FL1.A, ident, stain, run) %>%
    ggplot(aes(asinh.FL1.A)) +
    geom_density(aes(fill = ident, y = ..scaled..), alpha = 0.4) +
    guides(fill = FALSE) +
    ylab("norm. count") + xlim(c(5, 14)) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8, direction = -1)

  S2.plot <- df1 %>%
    dplyr::filter(mode == "SYBR2") %>%
    select(asinh.FL1.A, ident, stain, run) %>%
    ggplot(aes(asinh.FL1.A)) +
    geom_density(aes(fill = ident, y = ..scaled..), alpha = 0.4) +
    guides(fill = FALSE) +
    ylab("") + xlim(c(5, 14)) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8, direction = -1)

  legend.plot <- df1 %>%
    select(asinh.FL1.A, ident, stain, run) %>%
    ggplot(aes(asinh.FL1.A)) + geom_density(aes(fill = ident), alpha = 0.4) +
    ylab("") + xlab("") +
    scale_fill_viridis(
      discrete = TRUE, alpha = 0.7, begin = 0, end = 0.8,
      labels = c("non-sporulating cells", "purified spores"),
      name = "", direction = -1
    ) +
    theme(legend.position = c(0.01, 0.95)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
    ylim(0, 100) +
    annotate(geom = "text", x = -10, y = 65, parse = TRUE, label = "bold('A')~Side~Scatter~(no~stain)", hjust = 0, size = 4) +
    annotate(geom = "text", x = -10, y = 55, parse = TRUE, label = "bold('B')~Forward~Scatter~(no~stain)", hjust = 0, size = 4) +
    annotate(geom = "text", x = -10, y = 45, parse = TRUE, label = "bold('C')~PI~stain", hjust = 0, size = 4) +
    annotate(geom = "text", x = -10, y = 35, parse = TRUE, label = "bold('D')~SYBR1~stain", hjust = 0, size = 4) +
    annotate(geom = "text", x = -10, y = 25, parse = TRUE, label = "bold('E')~SYBR2~stain", hjust = 0, size = 4) +
    geom_rect(aes(xmin = -10, xmax = 100, ymin = 0, ymax = 1), color = "white", size = 2, fill = "white")

  plot_grid(SSC.plot, FSC.plot, PI.plot, S1.plot, S2.plot, legend.plot, ncol = 3, labels = "AUTO")
  ggsave("fig/Figure_1.pdf", width = 9, height = 5, units = "in", dpi = 300)
}

#### Figure 2: - Measuring cells and spores together: ability to separate with GMM ####

figure_2 <- function() {
  set.seed(1)
  sample.var <- c("ident", "type", "stain", "t1", "t2", "t3", "time")
  fcsset2.1 <- flowCreateFlowSet(
    filepath = "data/f2/set5_cells+spores_30/",
    sample_variables = sample.var,
    additional_variable = 30, transformation = TRUE
  )
  fcsset2.2 <- flowCreateFlowSet(
    filepath = "data/f2/set5_cells+spores_90/",
    sample_variables = sample.var,
    additional_variable = 90, transformation = TRUE
  )

  ## ???
  # pData(fcsset2.1[[10]])
  # fcsset2.1[[10]]%>%
  #       Subset(., norm2Filter("FSC-H", "FSC-W", filterId = "norm_ssc.fsc", scale = 2)) %>%
  #       autoplot("asinh.FSC.A")


  # combinations of stain, channel, concentration and time to examine
  cases <- data.frame(
    type = rep(c(rep("PI", 3), rep("SYBR1", 3), rep("SYBR2", 3), rep("unstained", 2)), 2),
    stain = c(rep(c(1, 2, 4), 3), 0, 0),
    time = c(rep(30, 11), rep(90, 11)),
    channel = rep(c(
      rep("asinh.FL3.A", 3), rep("asinh.FL1.A", 6),
      "asinh.FSC.A", "asinh.SSC.A"
    ), 2)
  ) %>%
    unite(data = ., sep = ".", col = "merge") %>%
    with(merge) %>%
    "["(-c(21, 22)) # removing unstained at second time point, this was not measured

  f2.df <- rbind2(fcsset2.1, fcsset2.2) %>%
    Subset(rectangleGate(
      "asinh.FSC.A" = c(0, 15),
      "asinh.SSC.A" = c(0, 15),
      "asinh.FL1.A" = c(0, 15),
      "asinh.FL3.A" = c(0, 15)
    )) %>%
    Subset(., norm2Filter("FSC-H", "FSC-W", filterId = "norm_ssc.fsc", scale = 2)) %>%
    flowFcsToDf(.) %>%
    dplyr::filter(ident == "Cells+spores") %>%
    select(type, stain, time, asinh.FL1.A, asinh.FL3.A, asinh.FSC.A, asinh.SSC.A) %>%
    gather("channel", "value", 4:7)

  models.list <- f2.df %>%
    split(x = ., interaction(.$type, .$stain, .$time, .$channel)) %>%
    Filter(function(x) nrow(x) != 0, .) %>%
    lapply(., function(x) {
      x %>%
        dplyr::filter(value > 5) %>%
        with(value) %>%
        Mclust(., 2)
    })

  models.list.clean <- models.list[cases]

  # get distribution parameters
  distr.values <- lapply(models.list.clean, function(x) {
    c(x$parameters$mean, sqrt(firstsecondPair(x$parameters$variance$sigmasq)), x$parameters$pro)
  }) %>%
    do.call("rbind", .) %>%
    as_tibble() %>%
    magrittr::set_colnames(c("mu_1", "mu_2", "sd_1", "sd_2", "pro_1", "pro_2")) %>%
    mutate(
      diffM = mu_2 - mu_1, # difference between means
      pooledSD = sqrt((sd_1^2 + sd_2^2) / 2), # pooled standard deviation
      cutoff = findCutoffs(mu_1, sd_1, pro_1, mu_2, sd_2, pro_2) # cutoff values
    ) %>%
    cbind(., cases) %>%
    separate(cases, into = c("type", "stain", "time", "z", "channel"))

  plot_sup3 <- function(model.list, model.cases) {
    Map(function(x, y) {
      ggplot() +
        geom_density(aes(x = x$data), size = 1.2) +
        # geom_line(aes(x$data,x$classification-1),col="brown",alpha=0.3)+
        stat_function(aes(x$data),
          fun = sdnorm, n = 999,
          args = list(
            mean = x$parameters$mean[1],
            sd = sqrt(x$parameters$variance$sigmasq[1]),
            lambda = x$parameters$pro[1]
          ),
          col = "#F8766D", linetype = "dashed", size = 1.2
        ) +
        stat_function(aes(x$data),
          fun = sdnorm, n = 999,
          args = list(
            mean = x$parameters$mean[2],
            sd = sqrt(firstsecondElement(x$parameters$variance$sigmasq)),
            lambda = x$parameters$pro[2]
          ),
          col = "#00BFC4", linetype = "dashed", size = 1.2
        ) +
        scale_color_discrete(guide = FALSE) +
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)
        ) +
        xlab(y) +
        # xlab(paste("Stain:",cases$type[i],"Concentration:",cases$stain[i],
        #           "Time:",cases$time[i],"Channel:",cases$channel[i]))+
        ylab("")
    }, model.list, model.cases)
  }

  # splitting supplemental as requested by reviewer
  plot.list <- plot_sup3(models.list.clean, cases)
  plot.list1 <- plot_sup3(models.list.clean[1:11], cases[1:11])
  plot.list2 <- plot_sup3(models.list.clean[12:20], cases[12:20])

  ggsave(
    plot = do.call(plot_grid, c(plot.list1, ncol = 3)),
    filename = "suppl/Supplemental3A.pdf", width = 12, height = 15
  )

  ggsave(
    plot = do.call(plot_grid, c(plot.list2, ncol = 3)),
    filename = "suppl/Supplemental3B.pdf", width = 12, height = 11.25
  )

  ggsave(
    plot = do.call(plot_grid, c(plot.list, ncol = 4)),
    filename = "suppl/Supplemental3.pdf", width = 16, height = 24
  )

  PS2.A <- distr.values %>%
    ggplot(aes(stain, diffM, fill = interaction(channel, type))) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
    geom_errorbar(aes(ymin = diffM - pooledSD, ymax = diffM + pooledSD, group = interaction(channel, type)),
      position = position_dodge(), alpha = 0.8
    ) +
    ylab("Difference spore/cell distributions") + xlab("") +
    facet_grid(time ~ .) +
    scale_fill_discrete(name = "") +
    scale_x_discrete(label = c("unstained", "1x", "2x", "4x"))

  PS2.B <- read_csv("data/viability2.csv") %>%
    gather("tripl", "value", 5:7) %>%
    dplyr::filter(!is.na(value)) %>%
    ggplot(aes(fct_inorder(conc), as.numeric(value) / 2.10, col = stain)) +
    geom_point(aes(shape = type),
      alpha = 0.9,
      position = position_dodge(width = 0.3),
      size = 3
    ) +
    geom_line(
      stat = "summary", fun.y = "mean",
      aes(as.numeric(fct_inorder(conc)), as.numeric(value) / 2.10, linetype = type),
      alpha = 0.8, size = 1, position = position_dodge(width = 0.3)
    ) +
    facet_grid(time ~ .) +
    ylab("% CFU") + xlab("") +
    scale_color_discrete(name = "") +
    scale_shape_discrete(name = "", labels = c("non-sporulating cells", "purified spores")) +
    scale_linetype_discrete(name = "", labels = c("non-sporulating cells", "purified spores"))

  plot_grid(PS2.A, PS2.B, nrow = 2, labels = c("A", "B"))
  ggsave("suppl/Supplemental1.pdf", width = 8, height = 10, dpi = 300, units = "in")

  # Figure 2, predicted means with pooled standard deviations
  FIG2.1 <- distr.values %>%
    dplyr::filter(time == 30 & stain == 2 | time == 30 & stain == 0) %>%
    ggplot(aes(interaction(type, channel), diffM, fill = channel)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_errorbar(aes(ymin = diffM - pooledSD, ymax = diffM + pooledSD), width = 0.5) + xlab("") +
    ylab("Difference spore/cell distributions") + scale_color_discrete(name = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = c(0.6, 0.8), axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_x_discrete(label = c("SYBR1", "SYBR2", "PI", "FSC", "SSC"))

  # Figure 2, raw data panel
  f2.2 <- f2.df %>%
    group_by(type, stain, time) %>%
    mutate(id = row_number()) %>%
    tidyr::spread(channel, value, fill = NA) %>%
    ungroup() %>%
    dplyr::select(-id)

  # translate so that cutoff is at 0
  f2.2$asinh.FL1.A[f2.2$type == "SYBR1"] <- f2.2$asinh.FL1.A[f2.2$type == "SYBR1"] - distr.values$cutoff[5]
  f2.2$asinh.FL1.A[f2.2$type == "SYBR2"] <- f2.2$asinh.FL1.A[f2.2$type == "SYBR2"] - distr.values$cutoff[8]
  f2.2$asinh.FSC.A <- f2.2$asinh.FSC.A - distr.values$cutoff[10]
  f2.2$asinh.SSC.A <- f2.2$asinh.SSC.A - distr.values$cutoff[11]
  f2.2$asinh.FL3.A <- f2.2$asinh.FL3.A - distr.values$cutoff[2]

  FIG2.2 <- f2.2 %>%
    gather("channel", "value", 4:7) %>%
    dplyr::mutate(channel = paste(channel, type, sep = ".")) %>%
    dplyr::select(-type) %>%
    dplyr::filter(channel %in% c(
      "asinh.FL1.A.SYBR1",
      "asinh.FL1.A.SYBR2",
      "asinh.FL3.A.PI",
      "asinh.SSC.A.unstained",
      "asinh.FSC.A.unstained"
    )) %>%
    dplyr::filter(!is.na(value)) %>%
    ggplot(aes(x = value, y = channel)) +
    geom_density_ridges(alpha = 0.5) +
    geom_vline(aes(xintercept = 0), col = "red") +
    xlim(c(-5, 5)) +
    xlab("translated scatter/fluorescence signal") + ylab("")

  plot_grid(FIG2.2, FIG2.1, rel_widths = c(0.60, 0.4), labels = c("A", "B"))
  ggsave("fig/Figure_2.pdf", width = 10, height = 5)
}

#### Figure 3+4: Clustering ####

figure_34 <- function() {
  set.seed(1)
  sample.var <- c("strain", "time", "tripl")
  fcsset3 <- flowCreateFlowSet(filepath = "data/f3/set9/", sample_variables = sample.var, transformation = TRUE)

  df3 <- fcsset3 %>%
    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.W", scale.factor = 2)) %>%
    Subset(rectangleGate("asinh.FL1.A" = c(0, 15), "asinh.FL3.A" = c(0, 15))) %>%
    flowFcsToDf(.)

  # Load reference containing all subpopulations
  df3.ref <- ref.ds(strn = "Bs02003", time = "24h", tripl = 2, df = df3)
  df3.mix <- mclust::Mclust(data = df3.ref, G = 3)

  # getting centers for visualization and export
  centers.list.df <- t(df3.mix$parameters$mean)
  write.csv(centers.list.df, "suppl/centers_f3.csv")
  center.locs <- factor(df3.mix$classification, levels = c(1, 2, 3))

  # Prediction for other cases
  df3.list <- df3 %>%
    dplyr::select(asinh.FSC.A, asinh.SSC.A, asinh.FL1.A, strain, time, tripl) %>%
    split(., df3$strain)

  clusterB.predict <- c(
    predict.Mclust(df3.mix, newdata = df3.list[["Bs02003"]][, 1:3])$classification,
    predict.Mclust(df3.mix, newdata = df3.list[["Bs02025"]][, 1:3])$classification
  )

  clusterB.annot <- do.call(rbind, list(
    df3.list[["Bs02003"]],
    df3.list[["Bs02025"]]
  ))

  clusterB.pred <- cbind(clusterB.annot, cluster = clusterB.predict)

  clplot1 <- data.frame(df3.ref, cluster = center.locs) %>%
    ggplot(aes(asinh.SSC.A, asinh.FL1.A)) +
    geom_hex(aes(fill = factor(cluster,levels=c(2,1,3))), bins = 300) + # ,alpha=..ncount.. #order= ?
    geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
    xlim(c(10, 15)) + ylim(c(2.5, 15)) +
    scale_fill_viridis(
      discrete = TRUE, end = 0.8, label = c("Cells", "Forespores", "Spores"), name = "", direction = -1,
      guide = FALSE
    ) +
    scale_alpha_continuous(guide = FALSE) +
    theme_bw() +
    geom_point(aes(centers.list.df[1, 2], centers.list.df[1, 3]), col = "blue", size = 1) +
    geom_point(aes(centers.list.df[2, 2], centers.list.df[2, 3]), col = "blue", size = 1) +
    geom_point(aes(centers.list.df[3, 2], centers.list.df[3, 3]), col = "blue", size = 1)

  clplot2 <- data.frame(df3.ref, cluster = center.locs) %>%
    ggplot(aes(x = asinh.SSC.A, y = asinh.FSC.A)) +
    geom_hex(aes(fill = factor(cluster,levels=c(2,1,3))), bins = 300) +
    geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
    xlim(c(10, 15)) + ylim(c(9, 12)) +
    scale_fill_viridis(
      discrete = TRUE, end = 0.8, label = c("Cells", "Forespores", "Spores"), name = "", direction = -1,
      guide = FALSE
    ) +
    scale_alpha_continuous(guide = FALSE) +
    theme_bw() +
    geom_point(aes(centers.list.df[1, 2], centers.list.df[1, 1]), col = "blue", size = 1) +
    geom_point(aes(centers.list.df[2, 2], centers.list.df[2, 1]), col = "blue", size = 1) +
    geom_point(aes(centers.list.df[3, 2], centers.list.df[3, 1]), col = "blue", size = 1)

  plot_grid(clplot1, clplot2, align = "h")
  ggsave("fig/Figure_3.png", width = 8, height = 4, dpi = 900)

  # supplemental 4
  clusterB.pred %>%
    ggplot(aes(asinh.SSC.A, asinh.FL1.A)) +
    geom_hex(bins = 200) +
    # geom_hex(aes(col=as.factor(cluster)),bins=50)+
    facet_grid(tripl ~ strain) +
    xlim(c(10, 15)) + ylim(c(2.5, 14)) +
    scale_fill_viridis(guide = FALSE) +
    theme_bw()

  ggsave("suppl/Supplemental_4.pdf", width = 4, height = 6)

  c.count.B <- clusterB.pred %>%
    group_by(strain, time, cluster, tripl) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    group_by(strain, time, tripl) %>%
    mutate(perc.mean.count = 100 * count / sum(count))

  ccount2 <- c.count.B %>%
    ggplot(aes(strain, perc.mean.count,
      col = factor(cluster, levels = c(2, 1, 3)),
      shape = factor(cluster, levels = c(2, 1, 3))
    )) +
    geom_point(size = 3, position = position_dodge(1), stat = "identity") +
    ylab("Proportion / %") + xlab("time / h") +
    scale_color_viridis(
      labels = c("Cells", "Forespores", "Spores"), discrete = TRUE, end = c(0.8), direction = -1,
      alpha = 0.8
    ) +
    scale_shape_discrete(labels = c("Cells", "Forespores", "Spores")) +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"), legend.title = element_blank(),
      legend.position = c(0.20, 0.9)
    ) +
    xlab("")

  sample.var <- c("strain", "time", "stain", "tripl")
  fcsset3A <- flowCreateFlowSet(
    filepath = "data/f3/set8/",
    sample_variables = sample.var, transformation = TRUE
  )

  df3A <- fcsset3A %>%
    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.W", scale.factor = 2)) %>%
    Subset(rectangleGate("asinh.FL1.A" = c(0, 15), "asinh.FL3.A" = c(0, 15))) %>%
    flowFcsToDf() %>%
    dplyr::filter(stain != "unstained")

  # splitting and applying gmm
  centers.list <- list(
    df3.mix$classification,
    df3.mix$classification,
    df3.mix$classification
  )

  df3A.list <- df3A %>%
    select(asinh.FSC.A, asinh.SSC.A, asinh.FL1.A, strain, stain, time) %>%
    split(.$strain)

  cluster.predict.A <- c(
    predict(df3.mix, newdata = df3A.list[["Bs02003"]][, 1:3])$classification,
    predict(df3.mix, newdata = df3A.list[["Bs02018"]][, 1:3])$classification,
    predict(df3.mix, newdata = df3A.list[["Bs02020"]][, 1:3])$classification
  )

  cluster.annot.A <- do.call(rbind, list(
    df3A.list[["Bs02003"]],
    df3A.list[["Bs02018"]],
    df3A.list[["Bs02020"]]
  ))

  cluster.pred.A <- cbind(cluster.annot.A, cluster = cluster.predict.A)

  c.count.A <- cluster.pred.A %>%
    group_by(strain, time, cluster, stain) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    separate(time, c("time", "h"), 2) %>%
    mutate(time = as.numeric(time)) %>%
    group_by(strain, time, stain) %>%
    mutate(perc.count = 100 * count / sum(count))

  c.count.A %>% write.csv(file = "suppl/population_count.csv")

  ccount1 <- c.count.A %>%
    ggplot(aes(time, perc.count)) +
    geom_point(aes(col = as.factor(cluster), shape = as.factor(cluster), group = time),
      size = 2
    ) +
    stat_summary(aes(col = as.factor(cluster)),
      fun.y = "mean", geom = "line", size = 2
    ) +
    ylab("Proportion / %") + xlab("time / h") +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(strain ~ .) +
    scale_color_viridis(
      labels = c("Cells", "Forespores", "Spores"), discrete = TRUE, end = c(0.8), direction = -1,
      guide = FALSE, alpha = 0.8
    ) +
    scale_shape_discrete(guide = FALSE) +
    theme_minimal() +
    theme(panel.spacing = unit(1, "lines"), legend.title = element_blank()) +
    xlab("time / h")

  plot_grid(ccount1, ccount2, nrow = 1, rel_widths = c(0.6, 0.4), labels = c("A", "B"))

  ggsave("fig/Figure_4.pdf", width = 8, height = 5)
}

#### Supplemental: Combinations ####

supplemental_2 <- function() {
  remove(list = ls())
  source("src/functions.R")
  sample.var <- c("ident", "type", "stain", "conc", "trash")
  fcs.df1 <- flowCreateFlowSet(filepath = "data/reg/", sample_variables = sample.var) %>%
    Subset(., norm2Filter("FSC-H", "FSC-W", filterId = "norm_ssc.fsc", scale = 1)) %>%
    Subset(rectangleGate("asinh.FL1.A" = c(0, 15), "asinh.FL3.A" = c(0, 15))) %>%
    flowFcsToDf(.)

  fcs.df1 %>%
    select(type, stain, conc, asinh.FL1.A, asinh.FL3.A) %>%
    gather("channel", "value", 4:5) %>%
    ggplot(aes(value, interaction(stain, conc), fill = type)) +
    geom_density_ridges(alpha = 0.6) + ylab("") +
    scale_fill_discrete(name = "", labels = c("non-sporulating cells", "purified spores")) +
    theme(legend.position = c(0.45, 0.97)) +
    facet_grid(. ~ channel)

  ggsave("suppl/Supplemental_2.pdf", width = 10, height = 10)

  fcs.df2 <- fcs.df1 %>%
    group_by(type, stain, conc) %>%
    summarize(mean.fl1 = mean(asinh.FL1.A), mean.fl3 = mean(asinh.FL3.A)) %>%
    ungroup() %>%
    mutate(
      PI = c(0, 2, 1, 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0),
      SYBR1 = c(0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 0, 2, 1, 0),
      SYBR2 = c(0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 0, 1, 2)
    )

  fcs.df2 %>%
    dplyr::filter(PI == 0) %>%
    lm(as.numeric(mean.fl1) ~ -1 + as.numeric(SYBR2) * as.numeric(SYBR1) +
      as.numeric(SYBR2), data = .) %>%
    summary()

  fcs.df2 %>%
    dplyr::filter(SYBR2 == 0) %>%
    lm(as.numeric(mean.fl3) ~ -1 + as.numeric(SYBR1) * as.numeric(PI), data = .) %>%
    summary()

  fcs.df2 %>%
    dplyr::filter(SYBR1 == 0) %>%
    lm(as.numeric(mean.fl3) ~ -1 + as.numeric(SYBR2) * as.numeric(PI), data = .) %>%
    summary()

  write.csv2(fcs.df2, "suppl/mixdesign.csv")
}

#### Additional Figures ####

adfig <- function() {
  suppl6 <- lapply(1:3, function(x) {
    cluster.pred.A %>%
      # dplyr::filter(!stain=="unstained")%>%
      dplyr::filter(stain == x) %>%
      ggplot(aes(asinh.SSC.A, asinh.FL1.A)) +
      # geom_density_ridges()+
      geom_hex(aes(), bins = 200) +
      scale_fill_viridis() +
      facet_grid(time ~ strain) +
      xlim(c(10, 15)) + ylim(c(2.5, 14)) +
      # scale_fill_viridis()+
      theme_bw()
  })

  lapply(1, function(x) {
    cluster.pred.A %>%
      # dplyr::filter(!stain=="unstained")%>%
      dplyr::filter(stain == x) %>%
      ggplot(aes(asinh.SSC.A, asinh.FL1.A)) +
      # geom_density_ridges()+
      # geom_polygon()+
      geom_hex(aes(col = as.factor(cluster), fill = log(..count..)), bins = 30) +
      facet_grid(strain ~ time) +
      xlim(c(10, 15)) + ylim(c(2.5, 14)) +
      # scale_fill_viridis()+
      theme_bw() +
      scale_fill_viridis()
  })
}
