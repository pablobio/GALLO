#' Plot relationship between data and grouping factors
#'
#' Takes the output from find_genes_qtls_around_markers function and creates a chord plot with the relationship between groups
#' @param qtl_file The output from find_genes_qtls_around_markers function
#' @param x The first grouping factor, to be plotted in the left hand side of the chord plot
#' @param y The second grouping factor, to be plotted in the left hand side of the chord plot
#' @param grid.col A character with the grid color for the chord plot or a vector with different colors to be used in the grid colors. Note that when a color vector is provided, the length of this vector must be equal the number of sectors in the chord plot
#' @param gap A numeric value corresponding to the gap between the chord sectors
#' @param degree A numeric value corresponding to the starting degree from which the circle begins to draw. Note this degree is always reverse-clockwise
#' @param canvas.xlim The coordinate for the canvas in the x-axis. By default is c(-1,1)
#' @param canvas.ylim The coordinate for the canvas in the y-axis. By default is c(-1,1)
#' @param cex The size of the labels to be printed in the plot
#' @importFrom circlize circos.track
#' @importFrom circlize chordDiagram
#' @importFrom circlize circos.track
#' @importFrom circlize circos.clear
#' @importFrom circlize get.cell.meta.data
#' @importFrom circlize circos.text
#' @importFrom graphics par
#' @importFrom unbalhaar uh
#' @return A chords relating x and y
#' @examples
#' \donttest{data(QTLmarkers)
#' data(gffQTLs)
#' out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs,
#' marker_file=QTLmarkers, method = "qtl",
#' marker = "snp", interval = 500000, nThreads = NULL)
#' out.enrich<-qtl_enrich(qtl_db=gffQTLs, qtl_file=out.qtls,
#' qtl_type = "Name", enrich_type = "chromosome",
#' chr.subset = NULL, padj = "fdr",nThreads = NULL)
#' out.enrich$ID<-paste(out.enrich$QTL," - ","CHR",out.enrich$CHR,sep="")
#' out.enrich.filtered<-out.enrich[which(out.enrich$adj.pval<0.05),]
#' out.qtls$ID<-paste(out.qtls$Name," - ","CHR",out.qtls$CHR,sep="")
#' out.enrich.filtered<-out.enrich.filtered[order(out.enrich.filtered$adj.pval),]
#' out.qtls.filtered<-out.qtls[which(out.qtls$ID%in%out.enrich.filtered$ID[1:10]),]
#' out.qtls.filtered[which(out.qtls.filtered$Reference==
#' "Feugang et al. (2010)"), "color_ref"]<-"purple"
#' out.qtls.filtered[which(out.qtls.filtered$Reference==
#' "Buzanskas et al. (2017)"),"color_ref"]<-"pink"
#' color.grid<-c(rep("black",length(unique(out.qtls.filtered$Abbrev))),
#' unique(out.qtls.filtered$color_ref))
#' names(color.grid)<-c(unique(out.qtls.filtered$Abbrev),
#' unique(out.qtls.filtered$Reference))
#' relationship_plot(qtl_file=out.qtls.filtered, x="Abbrev",
#' y="Reference",cex=1,gap=5,degree = 90,
#' canvas.xlim = c(-5, 5), canvas.ylim = c(-3, 3),
#' grid.col = color.grid)}
#' @export

relationship_plot<-function (qtl_file, x, y, grid.col = "gray60", degree = 90,canvas.xlim = c(-2, 2), canvas.ylim = c(-2, 2),cex, gap){
requireNamespace("circlize")
chord.matrix <- matrix(data = 0, nrow = length(unique(qtl_file[,x])), ncol = length(unique(qtl_file[, y])), dimnames = list(unique(qtl_file[,x]), unique(qtl_file[, y])))

    for (i in seq_along(1:nrow(chord.matrix))) {
    pos.col <- which(colnames(chord.matrix) %in% qtl_file[which(qtl_file[,x] == rownames(chord.matrix)[i]), y])
    chord.matrix[i, pos.col] <- 1
    }
par(mar = c(0, 0, 0, 0))
circlize::circos.par(gap.after = c(rep(gap, nrow(chord.matrix) - 1), length(unique(qtl_file[, x])), rep(gap, ncol(chord.matrix) - 1), length(unique(qtl_file[, x]))), start.degree = degree, clock.wise = FALSE, canvas.xlim = canvas.xlim, canvas.ylim = canvas.ylim, track.height = 0.5)
circlize::chordDiagram(t(chord.matrix), order = c(rownames(chord.matrix), colnames(chord.matrix)), grid.col = grid.col, transparency = 0, annotationTrack = "grid", h.ratio = 0.8, diffHeight = -uh(1, "mm"), annotationTrackHeight = c(0.01, 0.05))
circlize::circos.track(track.index = 1, panel.fun = function(x,y) {sector.name = get.cell.meta.data("sector.index")
circlize::circos.text(circlize::CELL_META$xcenter,  circlize::CELL_META$ylim[1],  circlize::CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5), cex = cex, col = "black", font = 1)}, bg.border = NA)
circlize::circos.clear()
}
