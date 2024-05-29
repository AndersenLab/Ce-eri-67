
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))

library(tidyverse)


`%notin%` <- Negate(`%in%`)

### ggplot theme ####


font_size = 7

theme_cust <- theme_bw() +
  theme(
    plot.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = font_size,  color = "black"),
    legend.title =  ggplot2::element_text(size = font_size,  color = "black"),
    axis.title =  ggplot2::element_text(size = font_size,  color = "black"),
    axis.text =  ggplot2::element_text(size = font_size,  color = "black"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black"
    ),
    strip.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    text = ggplot2::element_text(family = "Helvetica")
  )

part_labels_size = 9

### colors ####



fix.color = c(
  "WBTransposon00000738" = "red",
  "TIR" = "darkblue",
  "-" = "lightskyblue",
  "+" = 'magenta',
  # inv group
  "N2-like" = "orange",
  #FFA500
  "single_eri67" = "springgreen4",
  "CB4856-like" = "blue",
  #0000FF
  "INVsosi1" = "purple",
  #8000FF
  
  "DEL1" = "gray",
  #808080
  "DEL2" = "#800020",
  # Burgundy
  "DEL3" = "chartreuse",
  #7FFF00
  "DEL1_DUP" = "lightpink2" # EEA2AD
)

#### general info ####

date_f = "20240528"



N2cds_exon_eri67_1 <-
  data.table::fread("../processed_data/ref_cds_exon_eri67.tsv")


pxg_data <- data.table::fread("../processed_data/pxg_data.tsv") %>%
  dplyr::mutate(
    tx = dplyr::case_when(
      grepl("C41D11.1", transcript) ~ paste0("eri-6", "[", sub("(C41D11.1)(.*)", "\\2", transcript), "]"),
      transcript == "C41D11.6.1" ~ "sosi-1",
      transcript == "C41D11.7.1" ~ "eri-7",
      TRUE ~ transcript
    ),
    tx = ifelse(tx == "eri-6[abcd]", "ERI-6 exons", tx)
  ) %>%
  dplyr::mutate(top_marker_allele = fct_relevel(top_marker_allele, c("REF", "ALT")))


strain_SVassign <-
  data.table::fread("../processed_data/strain_SVassign_20231124.tsv") %>%
  dplyr::select(strain, isotype, SV_type, sv2)


##############################
# Main Figs #####
##############################

####  Figure 1  #####


######  Figure 1b ####
eri67_exon_ws283 <-
  data.table::fread("../processed_data/eri67_exon_ws283.tsv")
eri67_CDS_ws283 <-
  data.table::fread("../processed_data/eri67_CDS_ws283.tsv")

plt_eri6_ws283 <- ggplot(data = eri67_exon_ws283, aes(
  xstart = start / 1e6,
  xend = end / 1e6,
  y = tx
)) +
  ggtranscript::geom_range(fill = "gray",
                           height = 0.25,
                           size = 0.1) +
  ggtranscript::geom_range(data = eri67_CDS_ws283 , aes(fill = strand), size =
                             0.1) +
  ggtranscript::geom_intron(
    data = ggtranscript::to_intron(eri67_exon_ws283, "tx"),
    arrow.min.intron.length = 1,
    size = 0.1,
    aes(strand = strand)
  ) +
  theme_cust +
  scale_x_continuous(
    breaks = seq(4.450000, 4.470000, 0.010000) ,
    limits = c(4.450000, 4.470000)
  ) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(4, 1, 0, 0), "mm"),
    axis.title.y = element_blank(),
    axis.text.y =  ggplot2::element_text(
      size = font_size,
      color = "black",
      face = "italic"
    )
  ) +
  scale_fill_manual(values = c("lightskyblue", 'magenta')) +
  ylab("Gene / isoform") +
  xlab("Genomic position (Mb)\non Chromosome I")





######  Figure 1c #####


inbred_manha_files  <-
  data.table::fread("../processed_data/summarized_AGGREGATE_mapping_inbred_chrI.tsv.gz") %>%
  dplyr::filter(CHROM == "I")


inbred_manha_files_local <- inbred_manha_files %>%
  dplyr::filter(transcript %in% c("WBGene00016561", "eri-6[e.1]", "eri-6[f.1]")) %>%
  dplyr::mutate(tx = ifelse(transcript == "WBGene00016561", "ERI-6 exons", transcript))

plt_manha_local <-
  ggplot(inbred_manha_files_local, aes(x = POS / 1e6, y = log10p, color =
                                         sig)) +
  geom_point(size = 0.1) +
  theme_cust +
  facet_wrap(. ~ tx, nrow = 1) +
  scale_color_manual(values = c("BF" = "gold2", "NONSIG" = "gray69")) +
  ylab(expression(-log[10](italic(p)))) +
  theme_cust +
  xlab("Local eQTL position (Mb) on Chromosome I") +
  theme(
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )



######  Figure 1d #####


eri_local_exp <-
  subset(pxg_data,
         tx %in% c("ERI-6 exons", "eri-6[e.1]", "eri-6[f.1]", "sosi-1", "eri-7")) %>%
  dplyr::distinct(strain, value, tx, top_marker_allele)

plt_4464670_pxg_eri6 <- ggplot(data = eri_local_exp,
                               aes(x = top_marker_allele,
                                   y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    size = 0.5,
    color = "gray"
  ) +
  geom_jitter(position = position_jitter(0.2),
              aes(color = top_marker_allele),
              size = 0.2)  +
  scale_color_manual(values = c(
    "REF" = "orange",
    "ALT" = "blue",
    "No" = "gray69",
    "Yes" = "black"
  )) +
  facet_grid(. ~ tx) +
  theme_cust +
  theme(
    legend.position = "none",
    strip.text.x = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    panel.spacing = unit(0.1, "line")
  ) +
  labs(x = "Genotype at I:4464670", y = "Expression") +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(-2, 6))



stats_plt_4464670_pxg_eri6 <-
  ggpubr::compare_means(
    value ~ top_marker_allele,
    data = eri_local_exp,
    group.by = c("tx"),
    p.adjust.method = "bonferroni",
    label = "p.signif",
    method = "wilcox.test"
  )

######  Figure 1e #####

inbred_manha_files_distant <- inbred_manha_files %>%
  dplyr::filter(
    transcript %in% c(
      "F56D6.16",
      "F56D6.17",
      "K02E2.6.1",
      "W04A4.2.1",
      "Y105C5A.14",
      "Y105C5A.25",
      "Y59H11AR.6",
      "Y82E9BL.18.1",
      "Y82E9BL.18.2",
      "ZK795.6"
    )
  ) %>%
  dplyr::group_by(transcript) %>%
  dplyr::mutate(nn = row_number()) %>%
  dplyr::mutate(biotype = ifelse(
    nn == 1 &
      transcript %in% c(
        "F56D6.16",
        "F56D6.17",
        "Y105C5A.14",
        "Y105C5A.25",
        "Y59H11AR.6",
        "ZK795.6"
      ),
    "pseudo",
    NA
  ))

plt_manha_distant <-
  ggplot(inbred_manha_files_distant,
         aes(x = POS / 1e6, y = log10p, color = sig)) +
  geom_point(size = 0.1) +
  theme_cust +
  facet_wrap(. ~ transcript, nrow = 2) +
  scale_color_manual(values = c("BF" = "plum4", "NONSIG" = "gray69")) +
  ylab(expression(-log[10](italic(p)))) +
  theme_cust +
  xlab("Distant eQTL position (Mb) on Chromosome I") +
  theme(
    plot.margin = unit(c(0, 0, 0, 2), "mm"),
    
    legend.position = "none",
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    panel.spacing = unit(0.1, "line")
  ) +
  geom_text(aes(label = biotype, x = 12.5, y = 10),
            size = font_size * 5 / 14,
            color = "black")


######  Figure 1f #####


exp_6abcd <- pxg_data %>%
  dplyr::filter(transcript == "C41D11_1abcd") %>%
  dplyr::distinct(strain, value) %>%
  dplyr::rename(value_abcd = value)

exp_6e_exp <- pxg_data %>%
  dplyr::filter(transcript == "C41D11.1e.1") %>%
  dplyr::distinct(strain, value) %>%
  dplyr::rename(value_e = value)

pxg_data_eri6_target2 <- pxg_data %>%
  dplyr::distinct(strain,
                  value,
                  GeneName,
                  transcript,
                  top_marker,
                  top_marker_allele) %>%
  dplyr::filter(transcript %in% c("F56D6.16", "Y105C5A.25")) %>%
  dplyr::left_join(exp_6abcd) %>%
  dplyr::left_join(exp_6e_exp)


fig_1_cor_abcd <- ggplot(data = pxg_data_eri6_target2,
                         aes(x = value_abcd,
                             y = value,
                             color = top_marker_allele)) +
  geom_point(size = 0.3) +
  facet_grid(transcript ~ .  , scales = "free_y") +
  theme_cust +
  scale_x_continuous(breaks = c(0, 2, 3, 4, 5)  , limits = c(2, 5)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  scale_color_manual(values = c("REF" = "orange", "ALT" = "blue")) +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    plot.title =   ggplot2::element_text(size = font_size,  color = "black"),
    plot.margin = unit(c(0, 1, 0, 2), "mm")
  ) +
  ggtitle("ERI-6 exons") +
  labs(y = "Expression", x = "Expression")

fig_1_cor_e <- ggplot(data = pxg_data_eri6_target2,
                      aes(x = value_e,
                          y = value,
                          color = top_marker_allele)) +
  geom_point(size = 0.3) +
  facet_grid(transcript ~ .  , scales = "free_y") +
  theme_cust  +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  scale_color_manual(values = c("REF" = "orange", "ALT" = "blue")) +
  theme(
    legend.position = "none",
    
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    plot.title =   ggplot2::element_text(size = font_size,  color = "black"),
    plot.margin = unit(c(0, 1, 0, 2), "mm")
  ) +
  ggtitle(expression(paste(italic('eri-6[e.1]'))))

#####

fig_1ab <- cowplot::plot_grid(
  NULL,
  plt_eri6_ws283,
  labels = c('', 'B'),
  rel_widths =  c(2.2, 1),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "tb",
  nrow = 1
)



fig_1cd <-
  cowplot::plot_grid(
    plt_manha_local,
    plt_4464670_pxg_eri6 ,
    labels = c('', 'D'),
    rel_widths =  c(1.2, 2),
    label_size = part_labels_size,
    label_fontfamily = "Helvetica",
    axis = "tb",
    align = "h",
    nrow = 1
  )


fig_1f <- cowplot::plot_grid(
  fig_1_cor_abcd ,
  fig_1_cor_e,
  labels = c('', ''),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "tb",
  align = "h",
  nrow = 1
)


fig_1ef <- cowplot::plot_grid(
  plt_manha_distant,
  fig_1f,
  labels = c('', 'F'),
  rel_widths =  c(2.2, 1),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  align = "h",
  nrow = 1
)


fig_1 <- cowplot::plot_grid(
  fig_1ab,
  fig_1cd,
  fig_1ef,
  labels = c('A', 'C', 'E'),
  rel_heights = c(5.5, 3, 4.5),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "lr",
  nrow = 3
)

ggsave(
  fig_1,
  filename = paste("../figures/Fig_1_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 140,
  width = 180
)
####  Figure 2  #####

cds_exon_eri67_3 <-
  data.table::fread("../processed_data/cds_exon_eri67_141023.tsv") %>%
  dplyr::mutate(name = gsub("Ss", "ss", name))

eri7_orien <- cds_exon_eri67_3 %>%
  dplyr::filter(gene_name == "WBGene00016566" & exon_id == 6) %>%
  dplyr::mutate(start = end + 500,
                end = start + 50)

eri6_orien <- cds_exon_eri67_3 %>%
  dplyr::filter(gene == "Ce-eri-6") %>%
  dplyr::group_by(sp_st_gene) %>%
  dplyr::filter(start == max(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(end = start - 1300,
                start = end - 50)


cds_exon_eri67_3poly3 <-
  data.table::fread("../processed_data/shade_eri67.tsv")

cds_exon_eri67_4 <-
  data.table::fread("../processed_data/TE_eri67.tsv")


gene_sturc_plt <- ggplot() +
  geom_polygon(
    data = cds_exon_eri67_3poly3,
    aes(x = Vx, y = V8, group = V9),
    fill = "grey",
    alpha = .5
  ) +
  geom_line(data = subset(subset(
    cds_exon_eri67_3,
    gene %in% c("Cbn-ref" , "Cbr-ref", "Ce-eri-7", "Ce-sosi-1" , "INV-sosi-1")
  )), aes(
    x = start ,
    y = V8 + 0.1,
    group = interaction(gene, sp_st_gene)
  )) +
  geom_line(data = subset(subset(cds_exon_eri67_3, gene %in% c("Ce-eri-6"))), aes(
    x = start ,
    y = V8 - 0.1,
    group = interaction(gene, sp_st_gene)
  )) +
  ggplot2::geom_rect(
    data = subset(subset(
      cds_exon_eri67_3,
      gene %in% c("Cbn-ref", "Cbr-ref", "Ce-eri-7", "Ce-sosi-1", "INV-sosi-1")
    )),
    ggplot2::aes(
      xmin = start ,
      xmax = end ,
      ymin = V8 + .2,
      ymax = V8  ,
      fill = strand
    ) ,
    size = 0.1,
    color = "black"
  ) +
  ggplot2::geom_rect(
    data = subset(subset(cds_exon_eri67_3, gene %in% c("Ce-eri-6"))),
    ggplot2::aes(
      xmin = start ,
      xmax = end ,
      ymin = V8 ,
      ymax = V8 - .2  ,
      fill = strand
    ) ,
    size = 0.1,
    color = "black"
  ) +
  #inv sosi1, dl238, eca396
  geom_segment(
    aes(
      x = 14476,
      xend = 20332,
      y = 5.2,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  geom_segment(
    aes(
      x = 16911,
      xend = 17897,
      y = 5.2,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  geom_segment(
    aes(
      x = 19067,
      xend = 20332,
      y = 5.2,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  geom_segment(
    aes(
      x = 16989,
      xend = 18252,
      y = 5.2,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  #inv eri6ad, dl238, eca396
  geom_segment(
    aes(
      x = 6650,
      xend = 12964,
      y = 5,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  geom_segment(
    aes(
      x = 7908,
      xend = 11687,
      y = 5,
      yend = 6
    ),
    linetype = 1,
    color = "gray",
    size = 0.3
  ) +
  #^^^TE
  geom_segment(
    data = subset(cds_exon_eri67_4, TE_STRAND == "+" & te_le != 40) ,
    aes(
      x = start,
      xend = end,
      y = V8 + 0.32,
      yend = V8 + 0.32,
      color = TE_family
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    data = subset(cds_exon_eri67_4, TE_STRAND == "-" & te_le != 40),
    aes(
      xend = start,
      x = end,
      y = V8 + 0.3,
      yend = V8 + 0.3,
      color = TE_family
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    data = subset(cds_exon_eri67_4, TE_STRAND == "+" & te_le == 40),
    aes(
      x = start,
      xend = end,
      y = V8 + 0.18,
      yend = V8 + 0.18,
      color = TE_family
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    data = subset(cds_exon_eri67_4, TE_STRAND == "-" & te_le == 40) ,
    aes(
      xend = start,
      x = end,
      y = V8 + 0.18,
      yend = V8 + 0.18,
      color = TE_family
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  #^^^TE
  #^^^orientation
  geom_segment(
    data = eri7_orien ,
    aes(
      xend = start,
      x = end,
      y = V8 + .1,
      yend = V8 + .1
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    data = eri6_orien ,
    aes(
      xend = end,
      x = start,
      y = V8 - .1,
      yend = V8 - .1
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    aes(
      xend = 1900,
      x = 1950,
      y = 11 + .1,
      yend = 11 + .1
    ),
    size = 0.5 ,
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  #^^^orientation
  theme_cust +
  theme(
    legend.position = c(0.8, 0.815),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.text = ggplot2::element_text(
      size = font_size,
      color = "black",
      face = "italic"
    )
  )  +
  ylab("Species - Strain") +
  xlab("bp")   +
  scale_y_continuous(
    breaks = c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
    labels = c(
      "Cbn-ref",
      "Cbr-ref",
      "Ce-XZ1516",
      "Ce-ECA36",
      "Ce-NIC526",
      "Ce-CB4856",
      "Ce-DL238",
      "Ce-ECA396",
      "Ce-JU2526",
      "Ce-ref-N2",
      "Ce-9 strains",
      "Ce-JU1400"
    )
  ) +
  scale_fill_manual(values = c("-" = "lightskyblue", "+" = 'magenta')) +
  guides(fill = "none",
         color = guide_legend(ncol = 2)) +
  scale_color_manual(
    values = c(
      "Ce000179" = "gold2",
      "CELE45" = "lightskyblue2",
      "CELETC2" = "orange",
      "CEREP1A" = "plum4",
      "Polinton-1_CB" = "red",
      "Polinton_CE_TIR" = "darkblue",
      "Tc4v" = "springgreen4",
      "Tc4" = "lightpink2"
    ),
    name = "TE family"
  ) +
  #JU1400
  geom_curve(
    aes(
      x = 6443,
      y = .7,
      xend = 9276,
      yend = .7
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature = 0.3,
    angle = 90
  ) +
  geom_curve(
    aes(
      x = 9276,
      y = .7,
      xend = 6443,
      yend = .7
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature =  -0.3,
    angle = 90
  ) +
  #cb4856
  geom_curve(
    aes(
      x = 14000,
      y = 7.2,
      xend = 15500,
      yend = 7.2
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature = 0.3,
    angle = 90
  ) +
  geom_curve(
    aes(
      x = 15500,
      y = 7.2,
      xend = 14000,
      yend = 7.2
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature =  -0.3,
    angle = 90
  ) +
  geom_segment(aes(
    xend = 18042,
    x = 19478,
    y = 7.2,
    yend = 7.2
  ),
  size = 0.5 ,
  color = "gray")  +
  geom_curve(
    aes(
      x = 18700,
      y = 7.15,
      xend = 15800,
      yend = 7.15
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature =  -0.3,
    angle = 90
  ) +
  geom_curve(
    aes(
      x = 15800,
      y = 7.15,
      xend = 18700,
      yend = 7.15
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature =  0.3,
    angle = 90
  ) +
  #   N2
  geom_curve(
    aes(
      x = 5200,
      y = 3.25,
      xend = 8000,
      yend = 3.25
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature = 0.3,
    angle = 90
  ) +
  geom_curve(
    aes(
      x = 8000,
      y = 3.25,
      xend = 5200 ,
      yend = 3.25
    ),
    arrow = arrow(length = unit(0.01, "npc"),),
    colour = "gray",
    curvature =   -0.3,
    angle = 90
  ) +
  geom_text(aes(
    y = 2.8,
    x = 2000,
    label = "eri-7",
    fontface = 4
  ), size = font_size * 5 / 14) +
  geom_text(aes(
    y = 3.12,
    x = 9930,
    label = "sosi-1",
    fontface = 4
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 6.12,
    x = 18930,
    label = "sosi-1",
    fontface = 4
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 2.6,
    x = 12851,
    label = "eri-6[e]",
    fontface = 4
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 2.6,
    x = 15031,
    label = "eri-6[f]",
    fontface = 4
  ) , size = font_size * 5 / 14) +
  geom_text(
    data = subset(cds_exon_eri67_3, name != "ss6"),
    aes(
      y = V8 - 0.3,
      x = (start + end) / 2,
      label = name
    ),
    size = font_size * 5 / 14
  ) +
  geom_text(aes(y = 4.35, x = 17250, label = "ss6"), size = font_size *
              5 / 14) +
  geom_text(
    data = subset(cds_exon_eri67_4, name != "tir2"),
    aes(
      label = name,
      x = (end + start) / 2,
      y = V8 + 0.5
    ),
    size = font_size * 5 / 14
  ) +
  geom_text(
    data = subset(cds_exon_eri67_4,  name == "tir2"),
    aes(
      label = name,
      x = -100 + (end + start) / 2,
      y = V8 - 0.1
    ),
    size = font_size * 5 / 14
  )   +
  geom_text(
    aes(
      y = 12.4,
      x = 3900,
      label = "ERI-6 exons",
      fontface = 2
    ),
    size = font_size * 5 / 14,
    color = "black"
  ) +
  geom_text(
    aes(
      y = 12.7,
      x = 2850,
      label = "Crick",
      fontface = 2
    ),
    size = font_size * 5 / 14,
    color = "lightskyblue"
  ) +
  geom_text(aes(y = 12.7, x = 3700, label = "/"),
            size = font_size * 5 / 14,
            color = "black") +
  geom_text(
    aes(
      y = 12.7,
      x = 4750,
      label = "Watson",
      fontface = 2
    ),
    size = font_size * 5 / 14,
    color = "magenta"
  )

ggsave(
  gene_sturc_plt,
  filename = paste("../figures/Fig_2_", date_f, ".png", sep = ""),
  units = "mm",
  height = 120,
  width = 180,
  dpi = 1200
)


####  Figure 3  ####
## Network ##
####  Figure 4  #####
#####  Figure 4a ####


inv_pxg_data3 <- pxg_data %>%
  dplyr::select(tx, strain, value) %>%
  dplyr::distinct() %>%
  dplyr::left_join(strain_SVassign, by = c("strain" = "isotype")) %>%
  dplyr::mutate(SV_type2 = fct_relevel(
    SV_type,
    c("single_eri67", "CB4856-like", "N2-like", "INVsosi1")
  ))



inv_pxg_data3_eri <-
  subset(inv_pxg_data3,
         tx %in% c("ERI-6 exons",  "eri-6[e.1]", "eri-6[f.1]", "sosi-1", "eri-7")) %>%
  dplyr::mutate(tx = fct_relevel(
    tx,
    c("eri-7", "ERI-6 exons", "sosi-1", "eri-6[e.1]", "eri-6[f.1]")
  ))



inv_pxg_plt2 <- ggplot() +
  geom_boxplot(
    data = inv_pxg_data3_eri,
    outlier.shape = NA,
    alpha = 0.5 ,
    aes(x = SV_type2, y = value, color = SV_type)
  ) +
  geom_jitter(
    data = subset(inv_pxg_data3_eri, SV_type == "single_eri67"),
    aes(
      x = SV_type2,
      y = value,
      color = sv2,
      fill = sv2
    ) ,
    position = position_jitter(0.2),
    shape = 21,
    size = 1
  ) +
  geom_jitter(
    data = subset(inv_pxg_data3_eri, SV_type != "single_eri67"),
    aes(
      x = SV_type2,
      y = value,
      color = sv2,
      fill = sv2
    ) ,
    position = position_jitter(0.2),
    shape = 21,
    size = .5
  )  +
  theme_cust +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    panel.spacing = unit(0.01, "line"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  facet_wrap(. ~ tx , scales = "free_y", nrow = 1) +
  xlab("Structural variation") +
  ylab("Expression\n(normalized TPM)") +
  scale_color_manual(values = fix.color) +
  scale_fill_manual(values = fix.color) +
  scale_y_continuous(breaks =  seq(-1, 5, 1))


#####  Figure 4b ####
rnaseq_trans_cis_split <-
  data.table::fread("../processed_data/rnaseq_trans_cis_split_20231101.tsv") %>%
  dplyr::left_join(strain_SVassign)

spliced_perc_plt <-
  ggplot(rnaseq_trans_cis_split,
         aes(
           x = fct_reorder(strain, mean_perc1),
           y = mean_perc1 ,
           color = SV_type
         )) +
  geom_point(size = .5) +
  theme_cust +
  xlab("207 strains") +
  ylab(paste("Percent of RNA-seq reads\nspliced to eri-7"))   +
  scale_color_manual(values = fix.color)   +
  scale_fill_manual(values = fix.color) +
  scale_y_continuous(breaks =  seq(0, 100, 10)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.margin = unit(0.3, "mm"),
    panel.grid.major.y = ggplot2::element_line(linewidth = .051, color =
                                                 "gray70"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.position = "none"
  )



#####  Figure 4c ####


inv_pxg_data3_targets <-
  subset(
    inv_pxg_data3,
    tx %in% c(
      "F56D6.16",
      "F56D6.17",
      "K02E2.6.1",
      "W04A4.2.1",
      "Y105C5A.14",
      "Y105C5A.25",
      "Y59H11AR.6",
      "Y82E9BL.18.1",
      "ZK795.6"
    )
  )

inv_pxg_plt4 <- ggplot() +
  geom_boxplot(
    data = inv_pxg_data3_targets,
    outlier.shape = NA,
    alpha = 0.5 ,
    aes(x = SV_type2, y = value, color = SV_type)
  ) +
  geom_jitter(
    data = subset(inv_pxg_data3_targets, SV_type == "single_eri67"),
    aes(
      x = SV_type2,
      y = value,
      color = sv2,
      fill = sv2
    ) ,
    position = position_jitter(0.2),
    shape = 21,
    size = 1
  ) +
  geom_jitter(
    data = subset(inv_pxg_data3_targets, SV_type != "single_eri67"),
    aes(
      x = SV_type2,
      y = value,
      color = sv2,
      fill = sv2
    ) ,
    position = position_jitter(0.2),
    shape = 21,
    size = .5
  )  +
  theme_cust +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 2, 0, 2), "mm"),
    panel.spacing = unit(0.01, "line"),
    legend.position = "none",
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  facet_grid(. ~ tx) +
  xlab("Structural variation")  +
  ylab("Expression\n(normalized TPM)") +
  scale_color_manual(values =  fix.color)  +
  scale_fill_manual(values = fix.color)




#####

fig_4a <- cowplot::plot_grid(
  inv_pxg_plt2,
  NULL,
  rel_widths = c(4, 1),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "b",
  align = "h",
  nrow = 1
)


fig_4b <- cowplot::plot_grid(
  NULL,
  spliced_perc_plt,
  labels = c('B', ''),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "b",
  align = "h",
  nrow = 1
)

fig_4 <- cowplot::plot_grid(
  fig_4a,
  fig_4b,
  inv_pxg_plt4,
  labels = c('A', "", "C"),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "lr",
  nrow = 3
)

ggsave(
  fig_4,
  filename = paste("../figures/Fig_4_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 150,
  width = 180,
  dpi = 1200
)





##############################
# Supp Figs #####
##############################
####  supp 1  #####
nemascan_eQTL <-
  data.table::fread("../processed_data/nemascan_eQTL.tsv")


nemascan_eQTL2 <- nemascan_eQTL %>%
  dplyr::mutate(
    eQTL_classification = ifelse(
      eQTL_Chr == "I" &
        transcript %in% c(
          "C36A4.11.1",
          "F52D2.5.1",
          "F56D6.16",
          "F56D6.17",
          "K02E2.6.1",
          "W04A4.2.1",
          "W04B5.1",
          "Y105C5A.14",
          "Y105C5A.25",
          "Y59H11AR.6",
          "Y82E9BL.18.1",
          "Y82E9BL.18.2",
          "ZK795.6"
        ),
      "dis21",
      eQTL_classification
    )
  )

nemascan_eQTL2$transcript_Chr[nemascan_eQTL2$transcript_Chr == "MtDNA"] <-
  "M"
nemascan_eQTL2$eQTL_Chr[nemascan_eQTL2$eQTL_Chr == "MtDNA"] <- "M"

nemascan_eQTL2$Chr_pos <-
  factor(nemascan_eQTL2$transcript_Chr,
         levels = c("M", "X", "V", "IV", "III", "II", "I"))
nemascan_eQTL2$eChr_pos <-
  factor(nemascan_eQTL2$eQTL_Chr,
         levels = c("I", "II", "III", "IV", "V", "X", "M"))

#####  Figure S1a1  #####
fig_S1a <- ggplot()  +
  geom_point(
    data = subset(
      nemascan_eQTL2,
      mapping == "LOCO" &
        eQTL_classification != "dis21"
    ),
    aes(x = eQTL_peak / 1E6, y = transcript_start, color = eQTL_classification),
    size = 0.25,
    alpha = 0.5
  )  +
  geom_point(
    data = subset(
      nemascan_eQTL2,
      mapping == "LOCO" &
        eQTL_classification == "dis21"
    ),
    aes(x = eQTL_peak / 1E6, y = transcript_start, color = eQTL_classification),
    size = 0.5
  )  +
  scale_color_manual(values = c(
    "Distant eQTL" = "plum4",
    "Local eQTL" = "gold2",
    "dis21" = "red"
  )) +
  facet_grid(
    cols = vars(eChr_pos),
    rows = vars(Chr_pos),
    scales = "free",
    switch = "both"
  ) +
  theme_cust +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 1
    ),
    plot.margin = unit(c(1, 1, 0, 1), "mm"),
    panel.spacing = unit(0.01, "line"),
    axis.text = element_blank(),
    legend.position = "none"  ,
    axis.ticks = element_blank()
  )  +
  ylab("Transcript position (Mb)") +
  xlab("LOCO eQTL position (Mb)")

gt_S1a = ggplot_gtable(ggplot_build(fig_S1a))

gt_S1a$widths[18] = 0.3 * gt_S1a$widths[10]

gt_S1a$heights[7] = 0.3 * gt_S1a$heights[9]


fig_S1a_gt <- ggplotify::as.ggplot(gt_S1a)



#####  Figure S1a2  #####

fig_S1b <- ggplot()   +
  geom_point(
    data = subset(
      nemascan_eQTL2,
      mapping == "INBRED" &
        eQTL_classification != "dis21"
    ),
    aes(x = eQTL_peak / 1E6, y = transcript_start, color = eQTL_classification),
    size = 0.25,
    alpha = 0.5
  )  +
  geom_point(
    data = subset(
      nemascan_eQTL2,
      mapping == "INBRED" &
        eQTL_classification == "dis21"
    ),
    aes(x = eQTL_peak / 1E6, y = transcript_start, color = eQTL_classification),
    size = 0.5
  )  +
  scale_color_manual(values = c(
    "Distant eQTL" = "plum4",
    "Local eQTL" = "gold2",
    "dis21" = "red"
  )) +
  facet_grid(
    cols = vars(eChr_pos),
    rows = vars(Chr_pos),
    scales = "free",
    switch = "both"
  ) +
  theme_cust +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 1
    ),
    plot.margin = unit(c(1, 1, 0, 1), "mm"),
    panel.spacing = unit(0.01, "line"),
    axis.text = element_blank(),
    legend.position = "none"  ,
    axis.ticks = element_blank()
  )  +
  ylab("Transcript position (Mb)") +
  xlab("INBRED eQTL position (Mb)")


gt_S1b = ggplot_gtable(ggplot_build(fig_S1b))

gt_S1b$widths[18] = 0.3 * gt_S1b$widths[10]

gt_S1b$heights[7] = 0.3 * gt_S1b$heights[9]


fig_S1b_gt <- ggplotify::as.ggplot(gt_S1b)


#####



fig_S1 <- cowplot::plot_grid(
  fig_S1a_gt,
  fig_S1b_gt ,
  label_fontfamily = "Helvetica",
  label_size = part_labels_size,
  axis = "tblr",
  nrow = 1
)


#####  Figure S1b  #####

pxg_data_candidate <- pxg_data %>%
  dplyr::distinct(GeneName,
                  transcript,
                  tx,
                  strain,
                  value,
                  top_marker,
                  top_marker_allele,
                  fine_4464670)

pxg_data_eri6_candidate <- pxg_data_candidate %>%
  dplyr::filter(GeneName %in% c("eri-6", "sosi-1", "Y119C1B.10", "eri-7"))

pxg_data_eri6_candidate$tx2 <- factor(
  pxg_data_eri6_candidate$tx,
  levels = c(
    "eri-6[a.1]",
    "eri-6[b.1]",
    "eri-6[b.2]",
    "eri-6[c.1]",
    "eri-6[d.1]",
    "ERI-6 exons",
    "eri-6[e.1]",
    "eri-6[f.1]",
    "sosi-1",
    "eri-7"
  )
)

#plt_mk_pxg_all
plt_4464670_box_eri6 <- ggplot(data = pxg_data_eri6_candidate,
                               aes(x = tx2,
                                   y = value)) +
  geom_jitter(position = position_jitter(0.2),
              size = 0.2,
              alpha = 0.5)  +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c(
    "REF" = "orange",
    "ALT" = "blue",
    "No" = "gray69",
    "Yes" = "black"
  )) +
  theme_cust +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x =  ggplot2::element_text(
      size = font_size,
      color = "black",
      face = "italic"
    ),
    strip.text.x = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(1, 1, 0, 1), "mm")
  ) +
  labs(x = "Genotype at I:4464670", y = "Expression") +
  scale_y_continuous(breaks =  seq(-2, 6, 1))





#####
plt_supp_1 <- cowplot::plot_grid(
  fig_S1,
  plt_4464670_box_eri6,
  labels = c('A', 'B', 'C'),
  rel_heights = c(8, 4),
  ncol = 1,
  label_fontfamily = "Helvetica",
  label_size = part_labels_size,
  axis = "lr"
)




ggsave(
  plt_supp_1,
  filename = paste("../figures/plt_supp_1_", date_f, ".svg", sep = ""),
  units = "mm",
  height = 110,
  width = 150
)



#### supp 2  ####


finemap_data <-
  data.table::fread("../supp_files/TableS3_fine_mappings.tsv")

finemap_data$impact <-
  factor(finemap_data$VARIANT_IMPACT,
         levels = c("HIGH", "LOW", "Intergenic"))


### eri fine
finemap_data_eri <- finemap_data %>%
  dplyr::filter(transcript %in% c("C41D11_1abcd", "C41D11.1e.1", "C41D11.1f.1"))   %>%
  dplyr::mutate(
    tx = dplyr::case_when(
      grepl("C41D11.1", transcript) ~ paste0("eri-6", "[", sub("(C41D11.1)(.*)", "\\2", transcript), "]"),
      TRUE ~ "ERI-6 exons"
    ),
    tx = ifelse(tx == "eri-6[abcd]", "ERI-6 exons", tx)
  )




fine_eri_LOCO <- subset(finemap_data_eri, mapping == "LOCO")


plt_fine_eri_LOCO <- ggplot() +
  geom_point(data = fine_eri_LOCO,
             aes(x = POS / 1e6, y = VARIANT_LOG10p, color = impact),
             size = 0.1) +
  geom_point(
    data = dplyr::distinct(fine_eri_LOCO, transcript, eQTL_peak),
    aes(x = eQTL_peak / 1e6, y = 0),
    size = 1 ,
    alpha = 1,
    shape = 25,
    stroke = 0.5,
    color = "black",
    fill = "gold2"
  ) +
  geom_segment(
    data = subset(fine_eri_LOCO, POS == 4464670),
    arrow = arrow(length = unit(5, "points")),
    aes(
      x = POS / 1e6,
      xend = POS / 1e6,
      y = VARIANT_LOG10p,
      yend = VARIANT_LOG10p - 0.5
    ),
    color = "red"
  )  +
  facet_wrap(. ~ tx, scales = "free", nrow = 3) +
  theme_cust +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c(
    "HIGH" = "orange",
    "LOW" = "gray50",
    "Intergenic" = "gray80"
  )) +
  labs(x = "Genomic position (Mb)",
       y = expression(paste(-log[10](italic(
         p
       )), " in LOCO")))  +
  scale_x_continuous(breaks = seq(3, 11, 1))





fine_eri_INBRED <- subset(finemap_data_eri, mapping == "INBRED")


plt_fine_eri_INBRED <- ggplot() +
  geom_point(data = fine_eri_INBRED,
             aes(x = POS / 1e6, y = VARIANT_LOG10p, color = impact),
             size = 0.1) +
  geom_point(
    data = dplyr::distinct(fine_eri_INBRED, transcript, eQTL_peak),
    aes(x = eQTL_peak / 1e6, y = 0),
    size = 1 ,
    alpha = 1,
    shape = 25,
    stroke = 0.5,
    color = "black",
    fill = "gold2"
  ) +
  geom_segment(
    data = subset(fine_eri_INBRED, POS == 4464670),
    arrow = arrow(length = unit(5, "points")),
    aes(
      x = POS / 1e6,
      xend = POS / 1e6,
      y = VARIANT_LOG10p,
      yend = VARIANT_LOG10p - 0.5
    ),
    color = "red"
  )  +
  facet_wrap(. ~ tx, scales = "free" , nrow = 3) +
  theme_cust +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c(
    "HIGH" = "orange",
    "LOW" = "gray50",
    "Intergenic" = "gray80"
  )) +
  labs(x = "Genomic position (Mb)",
       y = expression(paste(-log[10](italic(
         p
       )), " in INBRED")))  +
  guides(color = guide_legend(ncol = 1, title = "Predicted impact of variant")) +
  scale_x_continuous(breaks = seq(3, 11, 1))






plt_fine_eri <- cowplot::plot_grid(
  plt_fine_eri_LOCO,
  plt_fine_eri_INBRED,
  rel_widths =    c(1, 1.6),
  labels = c('A', 'B'),
  ncol = 2,
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "tb"
)



ggsave(
  plt_fine_eri,
  filename = paste("../figures/plt_supp_2_", date_f, ".svg", sep = ""),
  units = "mm",
  height = 150,
  width = 170
)





#### supp 3  ####

finemap_data_targets <- finemap_data %>%
  dplyr::filter(!GeneName %in% c("eri-6", "sosi-1", "eri-7"))



fine_target_LOCO <- subset(finemap_data_targets, mapping == "LOCO")


plt_fine_target_LOCO <- ggplot() +
  geom_point(data = fine_target_LOCO,
             aes(x = POS / 1e6, y = VARIANT_LOG10p, color = impact),
             size = 0.1) +
  geom_point(
    data = dplyr::distinct(fine_target_LOCO, transcript, eQTL_peak),
    aes(x = eQTL_peak / 1e6, y = 0),
    size = 1 ,
    alpha = 1,
    shape = 25,
    stroke = 0.5,
    color = "black",
    fill = "gold2"
  ) +
  geom_segment(
    data = subset(fine_target_LOCO, POS == 4464670),
    arrow = arrow(length = unit(5, "points")),
    aes(
      x = POS / 1e6,
      xend = POS / 1e6,
      y = VARIANT_LOG10p,
      yend = VARIANT_LOG10p - 0.5
    ),
    color = "red"
  )  +
  facet_wrap(. ~ transcript, scales = "free", nrow = 7) +
  theme_cust +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c(
    "HIGH" = "orange",
    "LOW" = "gray50",
    "Intergenic" = "gray80"
  )) +
  labs(x = "Genomic position (Mb)",
       y = expression(paste(-log[10](italic(
         p
       )), " in LOCO")))  +
  scale_x_continuous(breaks = seq(3, 11, 1))






fine_target_INBRED <-
  subset(finemap_data_targets, mapping == "INBRED")


plt_fine_target_INBRED <- ggplot() +
  geom_point(data = fine_target_INBRED,
             aes(x = POS / 1e6, y = VARIANT_LOG10p, color = impact),
             size = 0.1) +
  geom_point(
    data = dplyr::distinct(fine_target_INBRED, transcript, eQTL_peak),
    aes(x = eQTL_peak / 1e6, y = 0),
    size = 1 ,
    alpha = 1,
    shape = 25,
    stroke = 0.5,
    color = "black",
    fill = "gold2"
  ) +
  geom_segment(
    data = subset(fine_target_INBRED, POS == 4464670),
    arrow = arrow(length = unit(5, "points")),
    aes(
      x = POS / 1e6,
      xend = POS / 1e6,
      y = VARIANT_LOG10p,
      yend = VARIANT_LOG10p - 0.5
    ),
    color = "red"
  )  +
  facet_wrap(. ~ transcript, scales = "free" , ncol = 2) +
  theme_cust +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(0, 2, 0, 2), "mm"),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c(
    "HIGH" = "orange",
    "LOW" = "gray50",
    "Intergenic" = "gray80"
  )) +
  labs(x = "Genomic position (Mb)",
       y = expression(paste(-log[10](italic(
         p
       )), " in INBRED")))  +
  guides(color = guide_legend(ncol = 1, title = "Predicted impact of variant")) +
  scale_x_continuous(breaks = seq(3, 11, 1))


plt_fine_target <- cowplot::plot_grid(
  plt_fine_target_LOCO,
  plt_fine_target_INBRED,
  labels = c('A', 'B'),
  ncol = 2,
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "tb"
)



ggsave(
  plt_fine_target,
  filename = paste("../figures/plt_supp_3_", date_f, ".png", sep = ""),
  units = "mm",
  height = 200,
  width = 170
)



#### supp 4 ####

pxg_data_target_candidate <- pxg_data_candidate %>%
  dplyr::filter(!GeneName %in% c("eri-6", "sosi-1", "Y119C1B.10", "eri-7", "Y73B6BL.288"))


unique(pxg_data_target_candidate$transcript)


plt_4464670_pxg_all <- ggplot(data = pxg_data_target_candidate,
                              aes(x = top_marker_allele,
                                  y = value)) +
  geom_jitter(position = position_jitter(0.2),
              size = 0.2,
              alpha = 0.5)  +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.5,
               aes(color = top_marker_allele)) +
  scale_color_manual(values = c(
    "REF" = "orange",
    "ALT" = "blue",
    "No" = "gray69",
    "Yes" = "black"
  )) +
  facet_wrap(transcript ~ ., scales = "free", nrow = 2) +
  theme_cust +
  theme(
    legend.position = c(0.93, 0.2),
    axis.text.x = element_blank(),
    strip.text.x = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    panel.spacing = unit(0.1, "line"),
    plot.margin = unit(c(0, 1, 0, 1), "mm")
  ) +
  labs(x = "Genotype at I:4464670", y = "Expression", color = "Genotype") +
  scale_y_continuous(breaks =  seq(-2, 6, 1))


ggsave(
  plt_4464670_pxg_all,
  filename = paste("../figures/plt_supp_4_", date_f, ".png", sep = ""),
  units = "mm",
  height = 60,
  width = 150
)



#### supp 5 ####

eri6_mpt2022_23 <-
  data.table::fread("../processed_data/eri6_cirspr_2223.tsv")

eri6_mpt2022_23_sum <- eri6_mpt2022_23 %>%
  dplyr::filter(!transcript == "C41D11.1f.1") %>%
  dplyr::mutate(tx = ifelse(transcript == "C41D11.1e.1", "eri-6[e]", "eri-6[abcd]"))  %>%
  dplyr::group_by(tx, sample) %>%
  dplyr::mutate(exp_tpm = sum(tpm)) %>%
  dplyr::select(-tpm,-transcript) %>%
  dplyr::distinct() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(background = fct_relevel(
    background,
    c("JU2141", "JU3144", "JU2106", "JU642", "ECA3096", "ECA3095")
  )) %>%
  dplyr::mutate(
    genotype_top = fct_relevel(genotype_top, c("REF", "ALT")),
    genotype_2nd = fct_relevel(genotype_2nd, c("REF", "ALT"))
  )  %>%
  dplyr::mutate(strain = fct_relevel(
    strain,
    c(
      "JU2141",
      "JU3144",
      "JU2106",
      "JU642",
      "ECA3096",
      "ECA3097",
      "ECA3099",
      "ECA3100",
      "ECA3095",
      "ECA3098",
      "ECA3101",
      "ECA3102",
      "ECA3617",
      "ECA3618",
      "ECA3619",
      "ECA3620",
      "ECA3621",
      "ECA3622",
      "ECA3623",
      "ECA3624"
    )
  ))  %>%
  dplyr::mutate(tx = fct_relevel(tx, c("eri-6[e]", "eri-6[abcd]")))




eri6_mpt2022_sum <- eri6_mpt2022_23_sum  %>%
  dplyr::filter(seq_pool == "pool2022")


eri_edit_plt1 <- ggplot(eri6_mpt2022_sum, aes(x = strain, y = exp_tpm)) +
  geom_boxplot(aes(color = genotype_top),
               outlier.alpha = 0,
               alpha = 0.5) +
  geom_jitter(size = 0.5, color = "black") +
  facet_grid(tx ~ background , scales = "free") +
  theme_cust +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.03, "line"),
    legend.box = "vertical",
    strip.text.y = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c("REF" = "orange", "ALT" = "blue"),
                     name = "Genotype at I:4464670") +
  xlab("Strain") +
  ylab("Expression")



eri6_mpt2023_sum <- eri6_mpt2022_23_sum  %>%
  dplyr::filter(seq_pool == "pool2023")  %>%
  dplyr::filter(!strain %in% c(
    "ECA3095",
    "ECA3096",
    "ECA3621",
    "ECA3622",
    "ECA3623",
    "ECA3624"
  ))


eri_edit_plt2 <- ggplot(eri6_mpt2023_sum, aes(x = strain, y = exp_tpm)) +
  geom_boxplot(aes(color = genotype_2nd),
               outlier.alpha = 0,
               alpha = 0.5) +
  geom_jitter(size = 0.5, color = "black") +
  facet_grid(tx ~ background , scales = "free") +
  theme_cust +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing = unit(0.03, "line"),
    strip.text.y = element_blank()
  ) +
  scale_color_manual(values = c("REF" = "orange", "ALT" = "blue"),
                     name = "Genotype at I:4464857") +
  xlab("Strain") +
  ylab("Expression")



eri6_mpt2023_sum2 <- eri6_mpt2022_23_sum  %>%
  dplyr::filter(seq_pool == "pool2023")  %>%
  dplyr::filter(strain %in% c(
    "ECA3095",
    "ECA3096",
    "ECA3621",
    "ECA3622",
    "ECA3623",
    "ECA3624"
  ))


eri_edit_plt3 <- ggplot(eri6_mpt2023_sum2, aes(x = strain, y = exp_tpm)) +
  geom_boxplot(aes(color = genotype_2nd),
               outlier.alpha = 0,
               alpha = 0.5) +
  geom_jitter(size = 0.5, color = "black") +
  facet_grid(tx ~ background , scales = "free") +
  theme_cust +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.03, "line"),
    axis.title.y = element_blank(),
    legend.box = "vertical",
    strip.text.y = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    )
  ) +
  scale_color_manual(values = c("REF" = "orange", "ALT" = "blue"),
                     name = "Genotype at I:4464857") +
  xlab("Strain") +
  ylab("Expression")



fig_edit_bc <- cowplot::plot_grid(
  eri_edit_plt2,
  eri_edit_plt3,
  labels = c('', 'C'),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "tb",
  nrow = 1
)



fig_edit <- cowplot::plot_grid(
  eri_edit_plt1,
  fig_edit_bc,
  labels = c('A', 'B'),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "lr",
  nrow = 2
)


ggsave(
  fig_edit,
  filename = paste("../figures/plt_supp_5_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 190,
  width = 170
)





#### supp 6  #####


eri67_pacbio_SV <-
  data.table::fread("../processed_data/eri67_pacbio_SV.csv") %>%
  dplyr::mutate(mutation = ifelse(mutation == "insertion", "insertion (bp)", mutation)) %>%
  dplyr::filter(!mutation == "duplication")

plt_sv2N2 <-  ggplot() +
  geom_vline(
    xintercept = 4456377 / 1e6,
    color = "grey69" ,
    linetype = 2
  ) +
  geom_vline(xintercept = 4457096 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(xintercept = 4459211 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(
    xintercept = 4460135 / 1e6,
    color = "grey69",
    linetype = 2
  ) +
  ggplot2::geom_rect(
    data = subset(eri67_pacbio_SV, mutation %in% c("deletion", "inversion")),
    ggplot2::aes(
      xmin = start / 1e6,
      xmax = end / 1e6,
      ymin = V8 - 0.2,
      ymax = V8 - 0 ,
      fill = mutation ,
      color = mutation
    ),
    alpha = 1
  ) +
  ggplot2::geom_rect(
    data = subset(eri67_pacbio_SV, mutation %in% c("duplication")),
    ggplot2::aes(
      xmin = start / 1e6,
      xmax = end / 1e6,
      ymin = V8,
      ymax = V8 + 0.2 ,
      fill = mutation,
      color = mutation
    ),
    alpha = 1
  ) +
  geom_point(
    data = subset(eri67_pacbio_SV, mutation == "insertion (bp)"),
    aes(
      x = start / 1e6,
      y = V8 + 0.2,
      fill = mutation,
      color =
        mutation
    ),
    shape = 25
  ) +
  theme_cust +
  scale_x_continuous(
    breaks = seq(4.450000, 4.470000, 0.005000) ,
    limits = c(4.450000, 4.470000)
  ) +
  ylab("Strain") +
  xlab("Genomic position (Mb)")  +
  theme(
    legend.position = c(0.91, 0.6),
    # legend.margin=unit(10,"mm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.spacing = unit(0.5, 'mm'),
    legend.box.background = element_rect(colour = "gray"),
    legend.title = element_blank()
  ) +
  scale_y_continuous(
    breaks = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
    labels = c(
      "Ce-ref-N2",
      "9strains",
      "JU1400",
      "JU2526",
      "ECA396",
      "CB4856",
      "DL238",
      "NIC526",
      "ECA36",
      "XZ1516"
    )
  ) +
  scale_fill_manual(
    values = c(
      #"deletion"="black","duplication"="plum4","inversion"="springgreen4","insertion"="red",
      "insertion (bp)" = "#1e466e",
      "inversion" = "#ffd06f",
      "deletion" = "#e76254",
      "#528fad"
    )
  ) +
  scale_color_manual(
    values = c(
      #"deletion"="black","duplication"="plum4","inversion"="springgreen4","insertion"="red",
      "insertion (bp)" = "#1e466e",
      "inversion" = "#ffd06f",
      "deletion" = "#e76254",
      "#528fad"
    )
  ) +
  geom_text(
    data = subset(
      eri67_pacbio_SV,
      mutation == "insertion (bp)" &
        length %notin% c(1078, 2820, 5218, 5220, 42, 45, 4089, 4124)
    ),
    aes(
      label = length ,
      x =   (start + 500) / 1e6 ,
      y = V8 + 0.4
    ),
    color = "#1e466e" ,
    size = font_size * 5 / 14
  ) +
  geom_text(
    data = subset(eri67_pacbio_SV, mutation == "insertion (bp)" &
                    length == 1078),
    aes(
      label = length ,
      x =   (start) / 1e6 ,
      y = V8 + 0.8
    ),
    color = "#1e466e" ,
    size = font_size * 5 / 14
  ) +
  geom_text(
    data = subset(eri67_pacbio_SV, mutation == "insertion (bp)" &
                    length == 2820),
    aes(
      label = length ,
      x =   (start) / 1e6 ,
      y = V8 - 0.6
    ),
    color = "#1e466e" ,
    size = font_size * 5 / 14
  ) +
  geom_text(
    data = subset(
      eri67_pacbio_SV,
      mutation == "insertion (bp)" &
        length %in% c(5218, 5220, 4089, 4124)
    ),
    aes(
      label = length ,
      x =   (start - 500) / 1e6 ,
      y = V8 + 0.4
    ),
    color = "#1e466e" ,
    size = font_size * 5 / 14
  ) +
  geom_text(
    data = subset(
      eri67_pacbio_SV,
      mutation == "insertion (bp)" & length %in% c(42, 45)
    ),
    aes(
      label = length ,
      x =   (start + 300) / 1e6 ,
      y = V8 + 0.4
    ),
    color = "#1e466e" ,
    size = font_size * 5 / 14
  ) +
  ggplot2::geom_rect(
    data = subset(N2cds_exon_eri67_1, strand == "+"),
    ggplot2::aes(
      xmin = start / 1e6 ,
      xmax = end / 1e6,
      ymin = 10,
      ymax = 10
    ),
    color = "magenta" ,
    size = 3
  ) +
  ggplot2::geom_rect(
    data = subset(N2cds_exon_eri67_1, strand == "-"),
    ggplot2::aes(
      xmin = start / 1e6 ,
      xmax = end / 1e6,
      ymin = 10,
      ymax = 10
    ),
    color = "lightskyblue" ,
    size = 3
  )  +
  geom_text(aes(
    y = 10.4  ,
    x = (10305500 - 5850692) / 1e6,
    label = "eri-7",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 10.4  ,
    x = (10308800 - 5850692) / 1e6,
    label = "ERI-6 exons"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 10.4  ,
    x = (10312000 - 5850692) / 1e6,
    label = "sosi-1",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 10.4  ,
    x = (10316000 - 5850692) / 1e6,
    label = "eri-6[ef]",
    fontface = 3
  ) ,
  size = font_size * 5 / 14)

#####  supp 3bc  #####

binding_sites <-
  data.table::fread("../processed_data/binding_sites.csv")

plt_dr_tir <- ggplot() +
  geom_vline(
    xintercept = 4456377 / 1e6,
    color = "grey69" ,
    linetype = 2
  ) +
  geom_vline(xintercept = 4457096 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(xintercept = 4459211 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(
    xintercept = 4460135 / 1e6,
    color = "grey69",
    linetype = 2
  ) +
  ggplot2::geom_rect(
    data = subset(binding_sites, mutation != "TF"),
    ggplot2::aes(
      xmin = start / 1e6,
      xmax = end / 1e6,
      ymin = V8,
      ymax = V8 + 0.2 ,
      fill = strain ,
      color = strain
    ),
    
    alpha = 1
  ) +
  theme_cust +
  scale_x_continuous(
    breaks = seq(4.450000, 4.470000, 0.005000) ,
    limits = c(4.450000, 4.470000)
  ) +
  scale_y_continuous(limits = c(0.5, 2.5)) +
  ylab("Direct\nrepeats") +
  xlab("Genomic position (Mb)") +
  scale_fill_manual(values = c(
    "DR-930bp" = "darkblue",
    "DR-744bp" = "darkblue"
  )) +
  scale_color_manual(values = c(
    "DR-930bp" = "darkblue",
    "DR-744bp" = "darkblue"
  )) +
  theme(
    axis.text.y =   element_blank(),
    axis.ticks.y = element_blank(),
    axis.text  = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  geom_text(aes(y = 2.1, x = 4452000 / 1e6, label = "~ 930 bp direct repeats"),
            size = font_size * 5 / 14) +
  geom_text(aes(y = 1.2, x = 4453000 / 1e6, label = "~ 715 bp repeats with transposon origins"),
            size = font_size * 5 / 14)


plt_tf_bs <- ggplot() +
  geom_vline(
    xintercept = 4456377 / 1e6,
    color = "grey69" ,
    linetype = 2
  ) +
  geom_vline(xintercept = 4457096 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(xintercept = 4459211 / 1e6,
             color = "black",
             linetype = 2) +
  geom_vline(
    xintercept = 4460135 / 1e6,
    color = "grey69",
    linetype = 2
  ) +
  ggplot2::geom_rect(
    data = subset(binding_sites, mutation == "TF"),
    ggplot2::aes(
      xmin = start / 1e6,
      xmax = end / 1e6,
      ymin = V8,
      ymax = V8 + 0.2
    ),
    fill = "springgreen4",
    color = "springgreen4",
    alpha = 1
  ) +
  theme_cust +
  scale_x_continuous(
    breaks = seq(4.450000, 4.470000, 0.005000) ,
    limits = c(4.450000, 4.470000)
  ) +
  ylab("TF binding sites") +
  xlab("Genomic position (Mb)") +
  theme(axis.text.y  = element_blank(),  axis.ticks.y = element_blank())



####supp_Fig_3d ####
break_points4 <-
  data.table::fread("../processed_data/break_points4.tsv") %>%
  dplyr::left_join(strain_SVassign)



break_points4$SV_type2 <- factor(
  break_points4$SV_type,
  levels = c("single_eri67", "CB4856-like", "N2-like", "INVsosi1")
)


plt_break_point <-
  ggplot(data = subset(break_points4), aes(x = SV_type2, y = perc)) +
  geom_jitter(position = position_jitter(0.3),
              #shape=21,
              size = 1,
              aes(color = SV_type))  +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5 ,
    aes(color = SV_type),
    size = 0.5
  ) +
  theme_cust +
  facet_grid(. ~ pos) +
  scale_color_manual(values = fix.color, name = "Structural variation") +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )+
  ylab("Normalized counts\nof spanning reads")

#####

supp_Fig_6 <-
  cowplot::plot_grid(
    plt_sv2N2,
    plt_dr_tir,
    plt_tf_bs,
    plt_break_point,
    labels = c('A', 'B', 'C' , 'D'),
    rel_heights = c(3 , 0.7, 1.3, 1.8),
    label_size = part_labels_size,
    label_fontfamily = "Helvetica",
    axis = "lr",
    align = "v",
    nrow = 4
  )

ggsave(
  supp_Fig_6,
  filename = paste("../figures/plt_supp_6_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 160,
  width = 180,
  dpi=1200
)

#### supp 7 ####

##### supp 7a ####

sosi1_all <-
  data.table::fread("../processed_data/sosi1_all.tsv") %>%
  dplyr::mutate(nam = paste(Strain, Name, sep = "-"),)

sosi1_all_plt <- ggplot() +
  ggplot2::geom_rect(
    data = sosi1_all,
    ggplot2::aes(
      xmin = start ,
      xmax = end ,
      ymin = V8 + .2,
      ymax = V8  ,
      fill = strand
    ) ,
    size = 0.1,
    color = "black"
  ) +
  theme_cust +
  xlab("bp")   +
  ylab("") +
  scale_y_continuous(
    breaks = c(7, 6, 5, 4, 3, 2, 1),
    labels = c(
      "N2-ref"  ,
      "ECA396-ss1" ,
      "ECA396-ss2" ,
      "ECA396-ss3" ,
      "JU2526-ss4",
      "JU2526-ss5",
      "JU2526-ss6"
    )
  )  +
  geom_text(aes(
    y = 7.35,
    x = 1200,
    label = "sosi-1" ,
    fontface = 3
  ), size = font_size * 5 / 14) +
  scale_fill_manual(values = c("-" = "lightskyblue", "+" = 'magenta'))  +
  theme(legend.position = "none")



##### supp 7b ####

polinton1_all <-
  data.table::fread("../processed_data/polinton1_all.tsv") %>%
  dplyr::mutate(
    start = start - 2034 + 1,
    end = end - 2034 + 1 ,
    Direction = ifelse(
      Direction == "R" & Type == "TIR",
      "Y",
      ifelse(Direction == "Y" &
               Type == "TIR", "R", Direction)
    ) ,
    Name = gsub("Plt", "tir", Name),
    Name = dplyr::case_when(
      Type == "TE" & Name == "Tir1" ~ "plt1",
      Type == "TE"  & Name == "Tir5" ~ "plt2",
      Type == "TE"  & Name == "Tir6" ~ "plt3",
      Type == "TE"  & Name == "Tir7" ~ "plt4",
      TRUE ~ Name
    )
  )

polinton1_all_plt <- ggplot() +
  geom_segment(
    data = subset(polinton1_all, Direction == "Y"),
    aes(
      x = start,
      xend = end,
      y = V8,
      yend = V8 ,
      color = Type
    ),
    size = 1 ,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  geom_segment(
    data = subset(polinton1_all, Direction == "R"),
    aes(
      xend = start,
      x = end,
      y = V8,
      yend = V8 ,
      color = Type
    ),
    size = 1 ,
    arrow = arrow(length = unit(0.02, "npc"))
  )  +
  geom_segment(
    data = subset(polinton1_all, Direction == "N"),
    aes(
      x = start,
      xend = end,
      y = V8,
      yend = V8 ,
      color = Type
    ),
    size = 1
  ) +
  geom_segment(
    data = subset(polinton1_all, Type == "INTRON"),
    aes(
      x = start,
      xend = end,
      y = V8,
      yend = V8 ,
      color = Type
    ),
    size = 2
  ) +
  theme_cust +
  scale_color_manual(values = c("white", "red", "darkblue")) +
  scale_y_continuous(
    breaks = c(8, 7, 6, 5, 4, 3, 2, 1),
    labels = c(
      "N2-ref",
      "CB4856-plt1_tir1",
      "CB4856-tir2",
      "CB4856-tir3",
      "ECA396-tir4",
      "ECA396-plt2_tir5",
      "ECA396-plt3",
      "ECA396-plt4"
    )
  ) +
  theme(legend.position = "none") +
  labs(x = "bp", y = "") +
  geom_text(aes(y = 7.75, x = 9700, label = "(WBTransposon00000738)") ,
            size = font_size * 5 / 14) +
  geom_text(aes(
    y = 7.75,
    x = 6500,
    label = "Polinton-1_CB",
    fontface = 3
  ) ,
  size = font_size * 5 / 14) +
  geom_text(aes(y = 7.75, x = 370, label = "TIR") , size = font_size * 5 /
              14) +
  geom_text(aes(y = 7.75, x = 16872, label = "TIR") , size = font_size *
              5 / 14) +
  geom_text(aes(y = 8.15, x = 1370, label = "INT") , size = font_size *
              5 / 14) +
  geom_text(aes(y = 8.15, x = 13000, label = "pPolB1") , size = font_size *
              5 / 14) +
  geom_segment(
    aes(
      xend = 2018,
      x = 905,
      y = 8 + 0.28,
      yend = 8 + 0.28
    ),
    size = 1 ,
    color = "orange",
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  geom_segment(
    aes(
      xend = 15493,
      x = 11229,
      y = 8 + 0.28,
      yend = 8 + 0.28
    ),
    size = 1 ,
    color = "lightpink2",
    arrow = arrow(length = unit(0.01, "npc"))
  )


align_plt <- cowplot::plot_grid(
  sosi1_all_plt,
  polinton1_all_plt,
  labels = c('A', 'B'),
  #  rel_heights = c(2,1,1),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  axis = "lr",
  align = "v",
  nrow = 2
)

ggsave(
  align_plt,
  filename = paste("../figures/plt_supp_7_", date_f, ".png", sep = ""),
  units = "mm",
  height = 180,
  width = 170
)

#### supp 8 ####

load("../processed_data/treeNJ_550_eri67plus_250st.RData")

raw_tree <- treeNJ_550

tbl_tree_admix <- as.tibble(raw_tree) %>%
  dplyr::distinct(label) %>%
  na.omit() %>%
  dplyr::left_join(strain_SVassign, by = c("label" = "isotype"))  %>%
  dplyr::mutate(label2 = ifelse(
    label %in% c(
      "CB4856",
      "ECA396",
      "ECA36",
      "NIC526",
      "XZ1516",
      "JU2526",
      "N2",
      "JU1400",
      "JU1896",
      "ECA703",
      "ECA812"#, aa$strain
    ),
    label,
    NA
  ))

#haplot
tree24 <- left_join(raw_tree, tbl_tree_admix, by = 'label')

plt_tree_wgs <-
  ggtree::ggtree(tree24, layout = "rectangular", size = 0.1) +
  scale_color_manual(values = fix.color) +
  scale_fill_manual(values = fix.color) +
  ggtree::geom_tippoint(aes(color = sv2, fill = sv2), size = 0.5, shape =
                          21) +
  ggtree::geom_tiplab(
    size = 2,
    nudge_x = 0.004,
    aes(
      label = label2,
      color = sv2 ,
      subset = label %notin% c("XZ1516")
    )
  ) +
  ggtree::geom_tiplab(
    size = 2,
    nudge_y = 10,
    aes(
      label = label2 ,
      color = sv2 ,
      subset = label %in% c("XZ1516")
    )
  ) +
  theme(legend.position = "none",
        plot.margin = unit(c(2, 0, 2, 0), "mm"))


ggsave(
  plt_tree_wgs,
  filename = paste("../figures/plt_supp_8_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 140,
  width =  80
)

### supp 9 ####


polin_blast_pPolB_INT_all <-
  data.table::fread("../processed_data/polin_blast_pPolB_INT5.tsv")


polin_blast_pPolB_INT5 <- polin_blast_pPolB_INT_all %>%
  dplyr::filter(!strain %in% c("AF16", "VX34", "QX1410"))

polin_blast_pPolB_INT5$st2 <-
  factor(
    polin_blast_pPolB_INT5$strain,
    levels = c(
      "N2"  ,
      "DL226",
      "EG4725" ,
      "JU1395",
      "JU2600",
      "JU310" ,
      "MY2147",
      "MY2693",
      "NIC2" ,
      "QX1794" ,
      "JU1400" ,
      "JU2526",
      "ECA396",
      "DL238" ,
      "CB4856" ,
      "NIC526",
      "ECA36" ,
      "XZ1516"
    )
  )


plt_polin_blast_pPolB_INT <- ggplot() +
  geom_hline(data = polin_blast_pPolB_INT5,
             aes(yintercept = count, color = eri67),
             size = 0.3) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, TIR_strand == "+"),
    aes(
      xend = TIR_mend1,
      x = TIR_mstart1,
      y = count,
      yend = count
    ),
    color = "darkblue",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, TIR_strand == "-"),
    aes(
      xend = TIR_mstart1,
      x = TIR_mend1,
      y = count,
      yend = count
    ),
    color = "darkblue",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, pPolB_strand == "+"),
    aes(
      xend = pPolB_mend1,
      x = pPolB_mstart1,
      y = count,
      yend = count
    ),
    color = "lightpink2",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, pPolB_strand == "-"),
    aes(
      xend = pPolB_mstart1,
      x = pPolB_mend1,
      y = count,
      yend = count
    ),
    color = "lightpink2",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, pPolB_strand == "+"),
    aes(
      xend = pPolB_mend2,
      x = pPolB_mstart2,
      y = count,
      yend = count
    ),
    color = "red",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, pPolB_strand == "-"),
    aes(
      xend = pPolB_mstart2,
      x = pPolB_mend2,
      y = count,
      yend = count
    ),
    color = "red",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, INT_strand == "+"),
    aes(
      xend = INT_mend1,
      x = INT_mstart1,
      y = count,
      yend = count
    ),
    color = "orange",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INT5, INT_strand == "-"),
    aes(
      xend = INT_mstart1,
      x = INT_mend1,
      y = count,
      yend = count
    ),
    color = "orange",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  theme_cust +
  facet_wrap(. ~ st2, scales = "free_y" , ncol = 6) +
  theme(axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("kb") +
  ylab("copies") +
  scale_color_manual(values = c("yes" = "lightskyblue2", "no" = "gray69"))  +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 20000, 10000)  ,
                     labels = c("0", "10", "20"))


polin_blast_pPolB_INTcb <- polin_blast_pPolB_INT_all %>%
  dplyr::filter(strain %in% c("AF16", "VX34", "QX1410"))


plt_polin_blast_pPolB_cb <- ggplot() +
  geom_hline(
    data = polin_blast_pPolB_INTcb,
    aes(yintercept = count),
    color = "gray69",
    size = 0.3
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INTcb, pPolB_strand == "+"),
    aes(
      xend = pPolB_mend1 + 2000,
      x = pPolB_mstart1 + 2000,
      y = count,
      yend = count
    ),
    color = "lightpink2",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INTcb, pPolB_strand == "-"),
    aes(
      xend = pPolB_mstart1 + 2000,
      x = pPolB_mend1 + 2000,
      y = count,
      yend = count
    ),
    color = "lightpink2",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INTcb, INT_strand == "+"),
    aes(
      xend = INT_mend1 + 2000,
      x = INT_mstart1 + 2000,
      y = count,
      yend = count
    ),
    color = "orange",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  geom_segment(
    data = subset(polin_blast_pPolB_INTcb, INT_strand == "-"),
    aes(
      xend = INT_mstart1 + 2000,
      x = INT_mend1 + 2000,
      y = count,
      yend = count
    ),
    color = "orange",
    size = 0.7  ,
    arrow = arrow(length = unit(0.04, "npc"))
  ) +
  theme_cust +
  facet_wrap(. ~ strain, scales = "free_y" , ncol = 1) +
  theme(
    axis.text.y = element_blank(),
    axis.title  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  #  xlab("Kb")+
  scale_x_continuous(
    breaks = seq(0, 20000, 10000),
    limits = c(0, 25000) ,
    labels = c("0", "10", "20")
  )



polin_scan_plt <-
  cowplot::plot_grid(
    plt_polin_blast_pPolB_INT,
    plt_polin_blast_pPolB_cb,
    labels = c('A', 'B'),
    rel_widths = c(5, 1),
    label_size = part_labels_size,
    label_fontfamily = "Helvetica",
    axis = "tb",
    align = "h",
    nrow = 1
  )

ggsave(
  polin_scan_plt,
  filename = paste("../figures/plt_supp_9_", date_f, ".png", sep = ""),
  units = "mm",
  height = 120,
  width = 180
)

### supp 10 ####

wgs_split_reads_bin <-
  data.table::fread("../processed_data/wgs_SA_bin200bp.tsv")

wgs_split_reads_bin_I_inv <- wgs_split_reads_bin  %>%
  dplyr::mutate(
    split1_bin = ifelse(bin_start < SA_bin_start, bin_start, SA_bin_start),
    split2_bin = ifelse(bin_start < SA_bin_start, SA_bin_start, bin_start)
  ) %>%
  dplyr::filter(
    bin_chr == "I" &
      SA_bin_chr == "I" &
      split1_bin >= 4450000 &
      split2_bin <= 4470000 &
      strandness == "opposite"
  ) %>%
  dplyr::group_by(strain, SA_bin_start, bin_start, split1_bin, split2_bin) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 3) %>%
  dplyr::mutate(distance = abs(split2_bin - split1_bin)) %>%
  dplyr::filter(distance > 200) %>%
  dplyr::distinct(strain, SA_bin_start, bin_start, split1_bin, split2_bin) %>%
  dplyr::mutate(bin_range = paste(SA_bin_start, bin_start, sep = "_"))

wgs_split_reads_bin_I_inv_count  <- wgs_split_reads_bin_I_inv  %>%
  dplyr::group_by(SA_bin_start, bin_start) %>%
  dplyr::add_count(name = "strain_occur") %>%
  dplyr::ungroup() %>%
  dplyr::filter(strain_occur > 10) %>%
  dplyr::select(-strain) %>%
  dplyr::distinct()

wgs_split_reads_bin_I_inv_count2 <-
  wgs_split_reads_bin_I_inv_count %>%
  dplyr::select(-split1_bin, -split2_bin) %>%
  tidyr::gather(loc, pos, -strain_occur, -bin_range)  %>%
  dplyr::mutate(V8 = ifelse(loc == "bin_start", 0.97, 0.03)) %>%
  dplyr::mutate(SV_type = ifelse(strain_occur %in% c(13, 21), "INVsosi1", "single_eri67"))


# to polinton

SA_XA_bref_DL238_yALT3 <- wgs_split_reads_bin  %>%
  dplyr::filter(SA_bin_chr == "I",
                SA_read_pos > 10302516,
                SA_read_pos < 10319657) %>%
  dplyr::filter(!(read_pos >= 4456378 & read_pos <= 4457096)) %>%
  dplyr::filter(!(read_pos >= 4459211 & read_pos <= 4459920)) %>%
  dplyr::distinct(strain, read_id, SA_bin_chr, SA_bin_start, bin_start) %>%
  dplyr::group_by(strain, SA_bin_chr, SA_bin_start, bin_start) %>%
  dplyr::add_count(name = "alt_count") %>%
  dplyr::ungroup()  %>%
  dplyr::filter(alt_count > 1) %>%
  dplyr::left_join(wgs_split_reads_bin) %>%
  dplyr::mutate(read_pos = read_pos + 5850692) %>%
  dplyr::select(strain, read_pos, SA_read_pos) %>%
  dplyr::distinct() %>%
  dplyr::mutate(linegroup = row_number()) %>%
  tidyr::gather(read, pos, -strain,-linegroup) %>%
  dplyr::mutate(V8 = dplyr::case_when(read == "read_pos" ~ -0.03,
                                      TRUE ~ -1.01)) %>%
  dplyr::left_join(strain_SVassign)



WBTe738_fit <-
  data.table::fread("../processed_data/WBTransposon00000738.tsv") %>%
  dplyr::mutate(read_type = gsub(" ", "\n", read_type)) %>%
  dplyr::filter(read_type == "Chimeric\nalignment")

# plot


plt_wgs_inv <- ggplot() + theme_cust +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(2, .5, 0, 2), "mm"),
    axis.text.y =  ggplot2::element_text(
      size = font_size,
      color = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  geom_hline(
    yintercept = c(1, 0, -1),
    color = "grey",
    linewidth = .2
  ) +
  scale_color_manual(values = fix.color) +
  geom_line(
    data = wgs_split_reads_bin_I_inv_count2,
    aes(
      x = (pos + 5850692) / 1e6,
      y = V8,
      group = bin_range,
      color = SV_type
    ),
    linewidth = .5
  ) +
  geom_line(
    data = subset(SA_XA_bref_DL238_yALT3 , SV_type == "INVsosi1"),
    aes(
      x = pos / 1e6,
      y = V8,
      group = linegroup,
      color = SV_type
    ) ,
    size = .1
  ) +
  geom_line(
    data = subset(SA_XA_bref_DL238_yALT3, SV_type == "N2-like"),
    aes(
      x = pos / 1e6,
      y = V8,
      group = linegroup,
      color = SV_type
    ) ,
    size = .1
  ) +
  geom_line(
    data = subset(SA_XA_bref_DL238_yALT3, SV_type == "CB4856-like"),
    aes(
      x = pos / 1e6,
      y = V8,
      group = linegroup,
      color = SV_type
    ) ,
    size = .1
  ) +
  geom_segment(
    aes(
      xend = (4457096 + 5850692) / 1e6,
      x = (4456378 + 5850692) / 1e6,
      y = 0,
      yend = 0
    ),
    size = .5 ,
    color = "darkblue",
    arrow = arrow(length = unit(0.01, "npc"))
  )  +
  geom_segment(
    aes(
      xend = (4459920 + 5850692) / 1e6,
      x = (4459211 + 5850692) / 1e6,
      y = 0,
      yend = 0
    ),
    size = .5 ,
    color = "darkblue",
    arrow = arrow(length = unit(0.01, "npc"))
  ) +
  ggplot2::geom_rect(
    data = N2cds_exon_eri67_1,
    ggplot2::aes(
      xmin = start2 / 1e6 ,
      xmax = end2 / 1e6,
      ymin = 1,
      ymax = 1 ,
      color = strand
    ),
    size = 2
  ) +
  ggplot2::geom_rect(
    data = N2cds_exon_eri67_1,
    ggplot2::aes(
      xmin = start2 / 1e6 ,
      xmax = end2 / 1e6,
      ymin = 0,
      ymax = 0 ,
      color = strand
    ),
    size = 2
  ) +
  geom_segment(
    aes(
      xend = 10318008 / 1e6,
      x = 10313744 / 1e6,
      y = -1.05,
      yend = -1.05
    ),
    size = .5 ,
    color = "lightpink2",
    arrow = arrow(length = unit(0.01, "npc"))
  )  +
  geom_segment(
    data = WBTe738_fit,
    aes(
      x = start / 1e6,
      xend = end / 1e6,
      y = -1,
      yend = -1 ,
      color = name
    ),
    size = 1 ,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  ylab("Split WGS reads mapped to ref-N2") +
  xlab("Genomic regions on chromosome I") +
  scale_x_continuous(
    breaks = seq(4.450000 + 5.850692, 4.470000 + 5.850692, 0.005000)  ,
    expand = c(0.0050, 0)
  ) +
  scale_y_continuous(
    breaks = c(-1, 0, 1),
    limits = c(-1.1, 1.1),
    labels = c(
      "Chimeric\nalignment",
      "Primary\nalignment" ,
      "Chimeric\nalignment"
    )
  ) +
  geom_text(aes(
    y = 1.1,
    x = (2000 + 4451508 - 1 + 5850692) / 1e6,
    label = "eri-7",
    fontface = 3
  ), size = font_size * 5 / 14) +
  geom_text(aes(
    y = 1.1,
    x = (6543 + 4451508 - 1 + 5850692) / 1e6,
    label = "ERI-6 exons"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 1.1,
    x = (9930 + 4451508 - 1 + 5850692) / 1e6,
    label = "sosi-1",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 1.1,
    x = (13031 + 4451508 - 1 + 5850692) / 1e6,
    label = "eri-6[e]",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 1.1,
    x = (14831 + 4451508 - 1 + 5850692) / 1e6,
    label = "eri-6[f]",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -1.1  ,
    x = 10310800 / 1e6,
    label = "Polinton-1_CB",
    fontface = 3
  ) ,
  size = font_size * 5 / 14) +
  geom_text(aes(
    y = -1.1  ,
    x = 10315744 / 1e6,
    label = "pPolB1"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -1.1  ,
    x = 10302816 / 1e6,
    label = "TIR"
  ) , size = font_size * 5 / 14)  +
  geom_text(aes(
    y = -1.1  ,
    x = 10319257 / 1e6,
    label = "TIR"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 0.5  ,
    x = (4455800 + 5850692) / 1e6,
    label = "Single"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = 0.5  ,
    x = (4456800 + 5850692) / 1e6,
    label = "eri-6-7" ,
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(
    aes(
      y = 0.38  ,
      x = (4456100 + 5850692) / 1e6,
      label = "93 strains"
    ) ,
    size = font_size * 5 / 14,
    color = "springgreen4"
  ) +
  geom_text(aes(
    y = 0.5  ,
    x = (4461546 + 5850692) / 1e6,
    label = "INVsosi-1"  ,
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(
    aes(
      y = 0.38  ,
      x = (4462046 + 5850692) / 1e6,
      label = "13 + 21 strains"
    ) ,
    size = font_size * 5 / 14,
    color = "purple"
  )


#####  supp 10b  #####

load("../processed_data/I_4451194_4469460_perBase_coverage_SD_step200.RData")

zone_forDepth <-
  data.table::fread("../processed_data/zone_forDepth.tsv") %>%
  dplyr::mutate(
    V8 = ifelse(zone == "eri6f_3UTR", -24, -22),
    pos = (start + end) / 2e6,
    pos = ifelse(
      zone == "eri6f_3UTR",
      pos + 0.0008,
      ifelse(zone == "eri6f", pos - 0.0002, pos)
    )
  ) %>%
  dplyr::mutate(zone = gsub("TIR*", "DR", zone))

I_4451194_4469460_perBase_coverage3_mean <-
  I_4451194_4469460_perBase_coverage3 %>%
  dplyr::select(Pos_slide, mean_d2a, sd_d2a) %>%
  dplyr::distinct()

coverage_color <- I_4451194_4469460_perBase_coverage3 %>%
  dplyr::left_join(strain_SVassign) %>%
  dplyr::filter(Pos_slide < 4465193 | Pos_slide > 4465693)


coverage_color$SV_type2 <- factor(
  coverage_color$SV_type,
  levels = c("single_eri67", "CB4856-like", "N2-like", "INVsosi1"),
  labels = c(
    "single_eri67\n(92 strains)",
    "CB4856-like\n(20 strains)",
    "N2-like\n(394 strains)",
    "INVsosi-1\n(44 strains)"
  )
)



plt_coverage <- ggplot() +
  geom_line(
    data = subset(coverage_color, SV_type %in% c("N2-like")),
    aes(
      x = Pos_slide / 1e6,
      y = depth_window,
      group = isotype,
      color = SV_type
    ) ,
    size = .1,
    alpha = 0.5
  ) +
  geom_line(
    data = subset(coverage_color, SV_type %in% c("single_eri67")),
    aes(
      x = Pos_slide / 1e6,
      y = depth_window,
      group = isotype,
      color = SV_type
    ) ,
    size = .1,
    alpha = 0.5
  ) +
  geom_line(
    data = subset(coverage_color, SV_type %in% c("INVsosi1")),
    aes(
      x = Pos_slide / 1e6,
      y = depth_window,
      group = isotype,
      color = SV_type
    ) ,
    size = .1,
    alpha = 0.5
  ) +
  geom_line(
    data = subset(coverage_color, SV_type %in% c("CB4856-like")),
    aes(
      x = Pos_slide / 1e6,
      y = depth_window,
      group = isotype,
      color = SV_type
    ) ,
    size = .1,
    alpha = 0.5
  ) +
  # facet_grid( SV_type2~.)+
  ggplot2::geom_rect(
    data = subset(N2cds_exon_eri67_1, strand == "+"),
    ggplot2::aes(
      xmin = start / 1e6 ,
      xmax = end / 1e6,
      ymin = -20,
      ymax = -20
    ),
    color = "magenta" ,
    size = 3
  ) +
  ggplot2::geom_rect(
    data = subset(N2cds_exon_eri67_1, strand == "-"),
    ggplot2::aes(
      xmin = start / 1e6 ,
      xmax = end / 1e6,
      ymin = -20,
      ymax = -20
    ),
    color = "lightskyblue" ,
    size = 3
  )  +
  scale_color_manual(values = fix.color,
                     name = "Haplotype") +
  theme_cust +
  theme(legend.position = "none", plot.margin = unit(c(2, .5, 0, 2), "mm")) +
  scale_x_continuous(breaks = seq(4.450000, 4.470000, 0.005000)  ,
                     expand = c(0, 0)) +
  xlab("Genomic position (Mb) on chromosome I") +
  ylab("Relative coverage") +
  geom_segment(
    data = zone_forDepth,
    aes(
      xend = start / 1e6,
      x = end / 1e6,
      y = V8 - 14,
      yend = V8 - 14
    ),
    size = 0.5
  ) +
  geom_text(aes(
    y = -52,
    x = (2000 + 4451508 - 1) / 1e6,
    label = "eri-7",
    fontface = 3
  ), size = font_size * 5 / 14) +
  geom_text(aes(
    y = -52,
    x = (6543 + 4451508 - 1) / 1e6,
    label = "ERI-6 exons"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -52,
    x = (9930 + 4451508 - 1) / 1e6,
    label = "sosi-1",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -52,
    x = (13031 + 4451508 - 1) / 1e6,
    label = "eri-6[e]",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -52,
    x = (14631 + 4451508 - 1) / 1e6,
    label = "eri-6[f]",
    fontface = 3
  ) , size = font_size * 5 / 14) +
  geom_text(aes(
    y = -54,
    x = (15531 + 4451508 - 1) / 1e6,
    label = "3'UTR"
  ) , size = font_size * 5 / 14) +
  geom_text(aes(y = -54, x = 4.456737, label = "DR-L") , size = font_size *
              5 / 14) +
  geom_text(aes(y = -54, x = 4.459566, label = "DR-R") , size = font_size *
              5 / 14)

#####  supp 10c  #####
I_4451194_4469460_perBase_coverage_zone <-
  data.table::fread("../processed_data/I_4451194_4469460_zone_meanCoverage_step.tsv") %>%
  dplyr::mutate(zone = ifelse(zone == "eri6f3utr", "eri6f_3UTR", zone))  %>%
  dplyr::group_by(zone) %>%
  dplyr::arrange(zone_mean) %>%
  dplyr::mutate(st_ord = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    zone = gsub("TIR*", "DR", zone),
    zone = ifelse(zone == "sosi1", "sosi-1", zone),
    zone = ifelse(zone == "eri-6[e]", "sosi-1", zone),
    zone = ifelse(zone == "sosi1", "sosi-1", zone)
  ) %>%
  dplyr::left_join(strain_SVassign)

dr12_depth <- I_4451194_4469460_perBase_coverage_zone %>%
  dplyr::filter(grepl("DR", zone)) %>%
  dplyr::select(strain, zone, zone_mean, SV_type) %>%
  tidyr::spread(zone, zone_mean)

plt_dr12_depth <- ggplot(dr12_depth, aes(x = DR1, y = DR2, color = SV_type)) +
  geom_point(size = .5) +
  theme_cust +
  scale_color_manual(values = fix.color) +
  theme(legend.position = "none",
        plot.margin = unit(c(2, .5, 0, 2), "mm")) +
  ylab("Mean coverage in DR-R") +
  xlab("Mean coverage in DR-L")


I_4451194_4469460_perBase_coverage_zone$zone2 <-
  factor(
    I_4451194_4469460_perBase_coverage_zone$zone,
    levels =   c("DR1", "DR2", "sosi-1", "eri6e", "eri6f", "eri6f_3UTR"),
    labels = c(
      "DR-L",
      "DR-R",
      "sosi-1",
      "eri-6[e]",
      "eri-6[f]",
      "eri-6[f] 3'UTR"
    )
  )

plt_ext10c2 <-
  ggplot(
    data = subset(
      I_4451194_4469460_perBase_coverage_zone,
      zone %in% c("sosi-1", "eri6e", "eri6f", "eri6f_3UTR")
    ),
    aes(x = st_ord, y = zone_mean)
  ) +
  geom_jitter(position = position_jitter(0.3),
              #shape=21,
              size = 1,
              aes(color = SV_type))  +
  theme_cust +
  facet_wrap(. ~ zone2 , scales = "free_y", nrow = 1)  +
  scale_color_manual(values = fix.color) +
  scale_fill_manual(values = fix.color) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = ggplot2::element_text(
      size = font_size,
      vjust = 1,
      color = "black",
      face = "italic"
    ),
    legend.position = "none",
    plot.margin = unit(c(2, .5, 0, 2), "mm")
  ) +
  ylab("Mean coverage in each region") +
  xlab("550 strains")

plt_ext10c <- cowplot::plot_grid(
  plt_dr12_depth,
  plt_ext10c2,
  labels = c('', 'D'),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  rel_widths = c(1.4, 4),
  axis = "tb",
  align = "h",
  nrow = 1
)


#####

plt_supp10 <-
  cowplot::plot_grid(
    plt_wgs_inv,
    plt_coverage,
    plt_ext10c,
    labels = c('A', 'B' , 'C'),
    label_size = part_labels_size,
    label_fontfamily = "Helvetica",
    rel_heights = c(1, 1, 0.7),
    axis = "lr",
    # align = "v",
    nrow = 3
  )



ggsave(
  plt_supp10,
  filename = paste("../figures/plt_supp_10_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 170,
  width = 180
)


### supp 11 ####

inv_st_distribution <-
  data.table::fread("../processed_data/inv_st_distribution.tsv") %>%
  dplyr::left_join(strain_SVassign)

world <- map_data("world") %>%
  dplyr::filter(!region == "Antarctica")

strain_world <- inv_st_distribution %>%
  dplyr::mutate(strain2 = "other") %>%
  dplyr::filter(longitude > -155)


Veri6e7_wolrd <- ggplot() +
  geom_map(
    data = world,
    map = world,
    aes(x = long, y = lat, map_id = region),
    color = "gray51",
    fill = "white",
    size = 0.5
  ) +
  geom_point(
    data = strain_world,
    aes(x = as.numeric(longitude),
        y = as.numeric(latitude)),
    color = "black",
    shape = 21,
    size = 0.1,
    alpha = 0.1
  ) +
  ggrepel::geom_label_repel(
    data = strain_world,
    aes(
      x = as.numeric(longitude),
      y = as.numeric(latitude),
      label = strain2,
      fill = sv2,
      color = sv2
    ),
    box.padding = 0.01,
    label.padding = 0.15,
    max.overlaps = Inf,
    point.padding = 0.01,
    segment.alpha = 0.5,
    segment.size = 0.1,
    size = 0.01
  ) +
  cowplot::theme_map() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white' , color = NA),
    plot.background = element_rect(fill = 'white' , color = NA),
    legend.text = element_text(size = 12,  color = "black")
  ) +
  scale_color_manual(values =  fix.color)  +
  scale_fill_manual(values = fix.color)


## in Hawaii

hawaii <- map_data("world") %>%
  dplyr::filter(subregion == "Hawaii")

strain_hawaii <- inv_st_distribution %>%
  dplyr::mutate(strain2 = " ") %>%
  dplyr::filter(longitude < -155)


Veri6e7_hawaii <- ggplot() +
  geom_map(
    data = hawaii,
    map = hawaii,
    aes(x = long, y = lat, map_id = region),
    color = "gray51",
    fill = "white",
    size = 0.5
  ) +
  geom_point(
    data = strain_hawaii,
    aes(x = as.numeric(longitude),
        y = as.numeric(latitude)),
    color = "black",
    shape = 21,
    size = 0.1,
    alpha = 0.1
  ) +
  ggrepel::geom_label_repel(
    data = strain_hawaii,
    aes(
      x = as.numeric(longitude),
      y = as.numeric(latitude),
      label = strain2,
      fill = sv2,
      color = sv2
    ),
    box.padding = 0.01,
    label.padding = 0.15,
    max.overlaps = Inf,
    point.padding = 0.01,
    segment.alpha = 0.5,
    segment.size = 0.1,
    size = 0.01
  ) +
  scale_fill_manual(values = fix.color) +
  scale_color_manual(values = fix.color) +
  cowplot::theme_map() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(5, 0, 5, 0), "mm"),
    panel.background = element_rect(fill = 'white' , color = NA),
    plot.background = element_rect(fill = 'white' , color = NA),
    legend.text = element_text(size = 12,  color = "black"),
    legend.title =  element_blank()
  )   +
  geom_text(aes(label = "Hawaiian Islands", y = 20, x = -159), size = font_size *
              5 / 14)

geo_dist_plt <- cowplot::plot_grid(
  Veri6e7_wolrd,
  Veri6e7_hawaii,
  labels = c('A', 'B'),
  label_size = part_labels_size,
  label_fontfamily = "Helvetica",
  nrow = 2
)

ggsave(
  geo_dist_plt,
  filename = paste("../figures/plt_supp_11_", date_f, ".pdf", sep = ""),
  units = "mm",
  height = 170,
  width = 160
)

### supp 12 ####

st550_83targets_cover <-
  data.table::fread("../processed_data/st550_83targets_possess.tsv") %>%
  dplyr::filter(perc >= 50)


st550_83targets_cover_count <- st550_83targets_cover %>%
  dplyr::distinct(strain, ens_gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::count()   %>%
  dplyr::arrange(n)

target_con_plt <-
  ggplot(st550_83targets_cover_count, aes(n)) +
  geom_histogram(color = "black",
                 fill = "white" ,
                 size = .3) +
  theme_cust +
  scale_x_continuous(breaks =  c(59, 70, 83)) +
  scale_y_continuous(breaks =  c(0, 20, 40, 60, 80)) +
  ylab("Number of strains") +
  xlab("Number of ERI-6/7 siRNA target genes\nin the reference strain N2")

ggsave(
  target_con_plt,
  filename = paste("../figures/plt_supp_12_", date_f, ".png", sep = ""),
  units = "mm",
  height = 60,
  width = 70
)


#####