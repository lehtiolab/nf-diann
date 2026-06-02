#!/usr/bin/env Rscript
library(plotly)
library(glue)

args = commandArgs(trailingOnly=TRUE)
precursortable = args[1]
inputfnpath = args[2]
enzyme = args[3]

# ModSeq, charge, file combination is "a single precursor"
seqcol = 'Modified.Sequence'
chargecol = 'Precursor.Charge'
filenamecol = 'Run'
miscleavcol = 'missed_cleavages'


# For boxplots:
rtcol = 'RT'
precerrcol = 'Ms1.Apex.Mz.Delta'
precquantcol = 'Ms1.Area'
peakwidthcol = 'FWHM'
ms2quantcol = 'Precursor.Quantity'

# Sample col is joined from conditions
samplecol = 'sample'


feats = read.table(precursortable, header=T, sep="\t", comment.char = "", quote = "")
inputfn = read.table(inputfnpath, header=T, sep="\t", comment.char = "", quote = "")
inputfn$Run = tools::file_path_sans_ext(basename(inputfn$file_path))
feats_labeled = merge(feats, inputfn[,c('Run', samplecol)], by='Run', all=TRUE)

boxplot_stats = function(data, col) {
  summary_stats = data %>%
  summarise(lower = quantile({{col}}, 0.25, na.rm=T),
            middle = median({{col}}, na.rm=T),
            upper = quantile({{col}}, 0.75, na.rm=T),
            minval = min({{col}}, na.rm=T),
            maxval = max({{col}}, na.rm=T),
    )
  summary_stats$iqr = summary_stats$upper - summary_stats$lower
  summary_stats$whisk_min = pmax(summary_stats$minval, summary_stats$middle - summary_stats$iqr * 1.5)
  summary_stats$whisk_max = pmin(summary_stats$maxval, summary_stats$middle + summary_stats$iqr * 1.5)
  return(summary_stats)
}


amount_ms2 = read.table("filescans", sep="\t", header=F)
colnames(amount_ms2) = c('file', 'nr_scans')
amount_ms2$file = gsub('[^A-Za-z0-9_]', '.', amount_ms2$file)
nr_verts = length(unique(amount_ms2[[1]]))
vert_height = 200 * nr_verts + 200

colmap = list(
  sample=c(samplecol),
  file=c(filenamecol)
  )
## List of plot names and their axis label:
ptypes = list(
  retentiontime=c(rtcol, 'time(min)'),
  precerror=c(precerrcol, 'Precursor error (ppm)'),
  precquant=c(precquantcol, 'MS1 quant'),
  peakwidth=c(peakwidthcol, 'FWHM'),
  ms2quant=c(ms2quantcol, 'MS2 quant')
)


for (grouper in names(colmap)) {
  xcol =  colmap[[grouper]][1]
  feats_labeled[,xcol] = gsub('[^A-Za-z0-9_]', '.', feats_labeled[,xcol])
  precursors = aggregate(feats_labeled[c(seqcol)], by=feats_labeled[xcol], length)
  names(precursors) = c(grouper, 'precursorcount')

  miscleav = aggregate(feats_labeled[c(seqcol)], by=feats_labeled[c(miscleavcol, xcol)], length)
  names(miscleav) = c('missed_cleavage', grouper, 'nrprec')
  miscleav$text = glue('{miscleav$nrprec} precursors')
  # Subset so we can use this when plotting
  miscleav = subset(miscleav, missed_cleavage %in% c(0,1,2))
  if (grouper == 'file') {
    fake_3_mc = data.frame(missed_cleavage=3, file=unique(miscleav$file), nrprec=0, text=unique(miscleav$file))
    miscleav = rbind(miscleav, fake_3_mc)
    precursors = merge(precursors, amount_ms2, by=grouper)
    precursors$perc_id = precursors$precursorcount / precursors$nr_scans * 100
  }
  miscleav = merge(precursors, miscleav, by=grouper)
  miscleav$percent = miscleav$nrprec / miscleav$precursorcount * 100
  if (grouper == 'file') {
    mc_text_y = 50
    mc_text_size = 3
    write.table(precursors, glue('{grouper}__counttable_qc.txt'), row.names=F, quote=F, sep='\t')
    write.table(miscleav, glue('{grouper}__miscleav_qc.txt'), row.names=F, quote=F, sep='\t')
  } else {
    miscleav$percent = miscleav$nrprec / miscleav$precursorcount * 100
    mc_text_y = max(miscleav$percent) * 2/6
    mc_text_size = 4
  }
  

  ggp = ggplot(precursors) +
    geom_bar(aes(x=.data[[grouper]], y=precursorcount), stat='identity', position='dodge') + 
    coord_flip() + 
    ylab('# precursors') + 
    theme_bw() + 
    theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
          axis.text.x=element_text(size=10, angle=90),
          legend.text=element_text(size=10), legend.title=element_blank())
    if (grouper == 'file') {
      # plotly doesnt support hjust, so calculate position to be in middle
      ggp = ggp + theme(axis.text.y=element_blank()) + 
        geom_text(aes(x=.data[[grouper]], y=precursorcount / 2, label=.data[[grouper]]), size=3, colour="white")
    } else {
      ggp = ggp + theme(axis.text.y=element_text(size=10, angle=90))
    }

  p = ggplotly(ggp, width=600, height=vert_height) %>%
          layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  # Work around since plotly does not honor above legend.title=element_blank call
  p$x$layout$legend$title$text = ''
  htmlwidgets::saveWidget(p, glue('{grouper}__amount_precursors.html'), selfcontained=F)

  ## Missed cleavages plot
  miscleav$missed_cleavage = as.factor(miscleav$missed_cleavage)
  mcplot = ggplot(miscleav) +
      geom_bar(aes(x=.data[[grouper]], y=percent, fill=missed_cleavage, group=missed_cleavage), position='dodge', stat='identity') +
      # 0.9 is the default dodge (90% of 1, 1 used bc all same value) but when not spec -> no dodge at all?
      geom_text(position=position_dodge(width=0.9), aes(x=.data[[grouper]], y=mc_text_y, group=missed_cleavage, label=text), colour="black", size=mc_text_size, inherit.aes=T) +
      ylim(c(0, 100)) + ylab('% of precursors') +
      theme_bw() +
      theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
            legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
      coord_flip() 
    if (grouper == 'file') {
      # plotly doesnt support hjust, so calculate position to be in middle
      mcplot = mcplot + theme(axis.text.y=element_blank())
    } else {
      mcplot = mcplot + theme(axis.text.y=element_text(size=10, angle=90))
    }


  p = ggplotly(mcplot, width=600, height=vert_height) %>%
          layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  p$x$layout$legend$title$text = ''
  htmlwidgets::saveWidget(p, glue('{grouper}__missed_cleavages.html'), selfcontained=F)

  for (ptype in names(ptypes)) {
    if (ptypes[[ptype]][1] %in% colnames(feats_labeled)) {
      summary_stats <- feats_labeled[c(xcol, ptypes[[ptype]][1])] %>%
        group_by(.data[[xcol]]) %>%
        boxplot_stats(.data[[ptypes[[ptype]][1]]])

      # Plot using geom_crossbar, geom_errorbar
      p = ggplot(summary_stats, aes(x=.data[[ xcol ]], y=middle)) +
        geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), position=position_dodge(width=1) ) +
        geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15, position=position_dodge(width=1)) + xlab(xcol) + coord_flip()
      if (ptype == 'precerror') { p = p + geom_hline(yintercept=0, size=2) }
      ggp = p + ylab(ptypes[[ptype]][2]) + theme_bw() + 
        theme(axis.title=element_text(size=15), axis.text.x=element_text(size=10))
      if (grouper == 'file') {
        # plotly doesnt support hjust, so calculate position to be in middle
        ggp = ggp + theme(axis.text.y=element_blank()) + 
          geom_text(aes(x=.data[[xcol]], y=(get('whisk_min') + get('whisk_max')) / 2, label=.data[[xcol]]), position=position_nudge(x=0.1), size=3, colour="black")
      } else {
        ggp = ggp + theme(axis.text.y=element_text(size=10, angle=90))
      }
     p = ggplotly(ggp, width=600, height=vert_height)
     htmlwidgets::saveWidget(p, glue('{grouper}__{ptype}.html'), selfcontained=F)
    }
  }
}
