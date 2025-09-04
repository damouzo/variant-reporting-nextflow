#!/usr/bin/env Rscript
# Script: plotQCsmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("script.R <clean_table> <gene_name> <prot_file> <exon_file> <part_metadata_file>")
}

clean_table     <- args[1] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants.rds"
gene_name       <- args[2] # "DHX34"
prot_file       <- args[3] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Protein\\DHX34.gff"
exon_file       <- args[4] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Exon\\DHX34.tsv"
p_metadata_file <- args[5] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants_participantMetadata.tsv"

# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(clean_table)) missing_files <- c(missing_files, paste("Clean table:", clean_table))
if (!file.exists(prot_file)) missing_files <- c(missing_files, paste("Protein file:", prot_file))
if (!file.exists(exon_file)) missing_files <- c(missing_files, paste("Exon file:", exon_file))
if (!file.exists(p_metadata_file)) missing_files <- c(missing_files, paste("Participant metadata file:", p_metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")



# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(paletteer)
library(grid)

# Settings ----------------------------------------------------------------------
#box::use(./utils_plotLollipop[...])
lolliplot <- function(SNP.gr, features=NULL, ranges=NULL,
                      type="circle",
                      newpage=TRUE, ylab=TRUE, ylab.gp=gpar(col="black"),
                      yaxis=TRUE, yaxis.gp=gpar(col="black"), 
                      xaxis=TRUE, xaxis.gp=gpar(col="black"), 
                      legend=NULL, legendPosition='top', cex=1, 
                      dashline.col="gray80", 
                      jitter=c("node", "label"), 
                      rescale=FALSE, 
                      label_on_feature=FALSE,
                      lollipop_style_switch_limit=10,
                      ...){
    ## Legends Function ##########
    checkLegendPosition <- function(legendPosition){
        positions <- c('top', 'left', 'right')
        if(!is.list(legendPosition)){
            legendPosition <- list(
            position = match.arg(legendPosition, choices = positions)
            )
        }else{
            legendPosition$position <- 
            match.arg(legendPosition$position,
                        choices = positions)
        }
        return(legendPosition)
    }

    getMaxYlimNchar <- function(SNP.gr, minVal, types){
        if(!is.list(SNP.gr)){
            SNP.gr <- list(SNP.gr)
        }
        m <- lapply(SNP.gr[types!='pie'], function(.ele){
            .e <- .ele$score
            if(length(.e)>0){
            return(max(.e[!is.infinite(.e)], na.rm=TRUE))
            }else{
            return(0)
            }
        })
        m <- max(c(1, unlist(m, recursive = TRUE)))
        m <- nchar(as.character(round(m))) + 3.5
        if(m<minVal){
            return(minVal)
        }else{
            return(min(10, m))
        }
    }
    
    handleLegend <- function(legend, len, dat){
        if(length(legend)>0){
            if(is.character(legend)){
            if(!missing(dat)){
                if(is.list(dat)){
                if(length(legend)<length(dat)){
                    legend <- rep(legend, length(dat))[seq_along(dat)]
                }
                }else{
                dat <- list(dat)
                legend <- legend[1]
                }
                para <- c("shape", "color", "border", "alpha")
                preset <- list("circle", "white", "black", 1)
                shapeMap <- c("circle"=21, "square"=22, "diamond"=23, 
                            "triangle_point_up"=24, "triangle_point_down"=25)
                names(preset) <- para
                legend <- mapply(function(.legend, .dat){
                coln <- colnames(mcols(.dat))
                if(.legend %in% coln){
                    labels <- mcols(.dat)[, .legend]
                    gp <- lapply(para, function(.ele){
                    if(.ele %in% coln){
                        mcols(.dat)[, .ele]
                    }else{
                        rep(preset[[.ele]], length(.dat))
                    }
                    })
                    names(gp) <- para
                    names(gp)[names(gp)=="color"] <- "fill"
                    if(is.list(gp[["shape"]])){
                    warning('multiple shape for one legend is not accept. Using default')
                    gp[["shape"]] <- shapeMap[vapply(gp[["shape"]], FUN=function(.ele) .ele[1], FUN.VALUE = character(1L))]
                    }else{
                    gp[["shape"]] <- shapeMap[gp[["shape"]]]
                    }
                    gp <- as.data.frame(gp, stringsAsFactors=FALSE)
                    gp <- cbind(labels=labels, gp)
                    names(gp)[names(gp)=="shape"] <- "pch"
                    gp <- gp[!duplicated(gp[, "labels"]), ]
                    gp <- gp[order(gp[, "labels"]), ]
                    gp <- as.list(gp)
                }
                }, legend, dat, SIMPLIFY = FALSE)
            }
            }
            if(!is.list(legend)){
            tmp <- legend
            legend <- vector(mode = "list", length = len)
            legend[[len]] <- tmp
            rm(tmp)
            }else{
            if(length(legend)==1){
                tmp <- legend[[1]]
                legend <- vector(mode = "list", length = len)
                legend[[len]] <- tmp
                rm(tmp)
            }else{
                if("labels" %in% names(legend)){
                tmp <- legend
                legend <- vector(mode = "list", length = len)
                legend[[len]] <- tmp
                rm(tmp)
                }else{
                if(length(legend)<len){
                    length(legend) <- len
                }
                }
            }
            }
        }
        return(legend)
    }

    handleRanges <- function(ranges, SNP.gr, features, len){
        if(length(ranges)>0){
            stopifnot(inherits(ranges, c("GRanges", "GRangesList", "list")))
            if(is(ranges, "GRanges")){
            if(length(ranges)==1){
                ranges <- split(rep(ranges, len)[seq.int(len)],
                                seq.int(len))
            }else{
                ranges <- split(rep(ranges, len),
                                rep(seq.int(len), each=len))[seq.int(len)]
            }
            }else{## GRangesList
            if(length(ranges)!=len){
                ranges <- rep(ranges, seq.int(len))[seq.int(len)]
            }
            }
            stopifnot(length(ranges)==len)
        }else{
            if(is(features, "GRanges")){
            ranges <- split(range(unname(features), ignore.strand=TRUE)[rep(1, len)],
                            seq.int(len))
            }else{
            if(length(features)!=len){
                stop("if both SNP.gr and features is GRangesList,",
                    " the lengthes of them should be identical.")
            }
            ranges <- GRangesList(lapply(features, function(.ele){
                range(unname(.ele), ignore.strand=TRUE)}))
            }
        }
        return(ranges)
    }

    cutSNP <- function(SNP.gr, ranges, len){
        if(is(ranges, "GRanges")){
            for(i in seq.int(len)){
            range <- ranges[i]
            stopifnot(all(width(SNP.gr[[i]])==1))
            SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
            }
        }else{
            if(inherits(ranges, c("GRangesList", "list"))){
            for(i in seq.int(len)){
                range <- ranges[[i]]
                stopifnot(all(width(SNP.gr[[i]])==1))
                SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
            }
            }
        }
        return(SNP.gr)
    }

    ################
    # Other Plot Functions
    jitterLables <- function(coor, xscale, lineW, weight=1.2){
        if(weight==1.2) {
        stopifnot("Please sort your inputs by start position"= 
                    order(coor)==1:length(coor))
        }
        if(weight<0.5) return(coor)
        stopifnot(length(xscale)==2)
        pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
        pos.diff <- diff(c(0, pos, 1))
        idx <- which(pos.diff < weight*lineW)
        if(length(idx)<1){
            return(coor)
        }
        if(all(idx %in% c(1, length(pos)+1))){
            return(coor)
        }
        idx.diff <- diff(c(-1, idx))
        idx.grp <- rle(idx.diff)
        idx.grp$values[idx.grp$values==1] <- length(pos) + 1:sum(idx.grp$values==1)
        idx.grp <- inverse.rle(idx.grp)
        idx.grp.w <- which(idx.grp>length(pos))-1
        idx.grp.w <- idx.grp.w[idx.grp.w>0]
        idx.grp[idx.grp.w] <- idx.grp[idx.grp.w+1]
        idx.grp <- split(idx, idx.grp)
        flag <- as.numeric(names(idx.grp))>length(pos)
        idx.grp.mul <- lapply(idx.grp[flag], function(.ele){
            c(.ele[1]-1, .ele)
        })
        idx.grp.sin <- lapply(idx.grp[!flag], function(.ele){
            lapply(as.list(.ele), function(.ele){c(.ele-1, .ele)})
        })
        idx.grp.sin <- unlist(idx.grp.sin, recursive = FALSE)
        idx.grp <- c(idx.grp.mul, idx.grp.sin)
        
        adj.pos <- lapply(idx.grp, function(.ele){
            .ele <- .ele[.ele>0 & .ele<=length(pos)]
            this.pos <- pos[.ele]
            names(this.pos) <- .ele
            if(length(this.pos)%%2==1){
                center <- ceiling(length(this.pos)/2)
            }else{
                center <- length(this.pos)/2 + .5
            }
            if(length(this.pos)>5){ ## too much, how to jitter?
                this.pos <- this.pos + 
                    ((1:length(this.pos))-center) * (weight-.1) * 
                    lineW/ceiling(log(length(this.pos), 5))
            }else{
                this.pos <- this.pos + 
                    ((1:length(this.pos))-center) * (weight-.1) * lineW
            }
            this.pos
        })
        names(adj.pos) <- NULL
        adj.pos <- unlist(adj.pos)
        coor[as.numeric(names(adj.pos))] <- adj.pos*diff(xscale)+xscale[1]
        
        Recall(coor, xscale=xscale, lineW=lineW, weight=weight-0.2)
    }

    getDrawLabelParam <- function(SNPs){
        label.parameter.draw <- rep(TRUE, length(SNPs))
        if(length(SNPs$label.parameter.draw)==length(SNPs)){
            label.parameter.draw <- vapply(SNPs$label.parameter.draw, `[`, i=1,
                                        FUN.VALUE = logical(1L))
        }
        label.parameter.draw
    }

    handleLabelParams <- function(SNPs, prefix="label.parameter.", cex=1,
                                    ...){
        labels <- list(...)
        ## change the parameter by use definitions.
        for(label.parameter in names(labels)){
            label.para <- paste0(prefix, label.parameter)
            if(label.para %in% colnames(mcols(SNPs))){
            labels[[label.parameter]] <- mcols(SNPs)[, label.para]
            }
        }
        for(label.parameter in c("col", "fill", "alpha", "lty", "lwd", "lex",
                                "lineend", "linejoin", "linemitre",
                                "fontsize", "cex", "fontfamily", "fontface",
                                "lineheight", "font")){
            label.para <- paste0(prefix, label.parameter)
            if(label.para %in% colnames(mcols(SNPs))){
            labels$gp[[label.parameter]] <- mcols(SNPs)[, label.para]
            }
        }
        if(!"cex" %in% names(labels$gp)){
            labels$gp$cex <- cex
        }
        mergeList <- function(.ele){
            .ele <- do.call(list, .ele)
            .n <- unique(unlist(lapply(.ele, names)))
            .out <- list()
            if(length(.n)>0){
            for(.name in .n){
                .out[[.name]] <- sapply(.ele, function(.e){
                if(.name %in% names(.e)){
                    .e[[.name]][1]
                }else{
                    NA
                }
                })
            }
            }else{
            .n <- unique(names(.ele))
            for(.name in .n){
                .out[[.name]] <- unlist(.ele[names(.ele) %in% .name])
            }
            }
            .out
        }
        labels$gp <- mergeList(labels$gp)
        labels$gp[duplicated(names(labels$gp))] <- NULL
        labels$gp <- do.call(gpar, labels$gp)
        return(labels)
        }
        filterThisLabel <- function(this.label){
        na_label <- is.na(this.label$label)
        for(key in c("x", "just", "hjust", "vjust", "rot",
                    "check.overlap", "default.units")){
            if(length(this.label[[key]])>1 &&
            length(this.label[[key]]==length(this.label$label))){
            this.label[[key]] <- this.label[[key]][!na_label]
            }
        }
        for(key in names(this.label$gp)){
            if(length(this.label$gp[[key]])>1 &&
            length(this.label$gp[[key]])==length(this.label$label)){
            this.label$gp[[key]] <- this.label$gp[[key]][!na_label]
            }
        }
        this.label$gp <- do.call(gpar, this.label$gp)
        this.label$label <- this.label$label[!na_label]
        this.label
    }

    grid.lollipop <- function (x1=.5, y1=.5,
                                x2=.5, y2=.75,
                                y3=.04, y4=.02,
                                radius=.8,
                                col=NULL,
                                border=NULL,
                                percent=NULL,
                                edges=100,
                                type=c("circle", "pie", "pin", "pie.stack", "flag"),
                                ratio.yx=1,
                                pin=NULL,
                                scoreMax,
                                scoreType,
                                id=NA,
                                cex=1, lwd=1,
                                dashline.col="gray80",
                                side=c("top", "bottom"),
                                alpha=NULL,
                                shape=shape){
            side <- match.arg(side)
            stopifnot(is.numeric(c(x1, x2, y1, y2, y3, y4, radius, edges)))
            type <- match.arg(type)
            side <- side!="top"
            if(!type %in% c("pie", "pie.stack")){
                this.score <- if(length(percent$score)>0) percent$score else 1
                if(type=="circle"){
                    y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-1)*2*radius*ratio.yx)
                    if(scoreType) y0[4] <- y2+y3+y4
                    if(side) y0 <- 1 - y0
                    grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                            gp=gpar(col=border, lwd=lwd))
                    y0 <- c(y2+y3+y4+this.score*2*radius*ratio.yx, 
                            y2+y3+y4+scoreMax*ratio.yx)
                    if(scoreType) y0[1] <- y2+y3+y4+this.score*2*radius*ratio.yx
                    if(side) y0 <- 1 - y0
                    grid.lines(x=c(x2, x2), 
                            y=y0, 
                            gp=gpar(col=dashline.col, lty=3, lwd=lwd))
                }else{
                    y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-.5)*2*radius*ratio.yx)
                    if(side) y0 <- 1 - y0
                    grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                            gp=gpar(col=border, lwd=lwd))
                }
                
            }else{
                if(type=="pie.stack"){
                    if(percent$stack.factor.first){
                        y0 <- c(y1, y2, y2+y3, y2+y3+y4)
                        if(side) y0 <- 1 - y0
                        grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                                gp=gpar(col=border, lwd=lwd))
                        y0 <- c(y2+y3+y4, y2+y3+y4+scoreMax*ratio.yx)
                        if(side) y0 <- 1 - y0
                        grid.lines(x=c(x2, x2), 
                                y=y0,
                                gp=gpar(col=dashline.col, lty=3, lwd=lwd))
                    }
                }else{
                    y0 <- c(y1, y2, y2+y3, y2+y3+y4)
                    if(side) y0 <- 1 - y0
                    grid.lines(x=c(x1, x1, x2, x2), y=y0, 
                            gp=gpar(col=border, lwd=lwd))
                }
            }
            if(length(pin)>0){
                if(length(border)>0) pin@paths[[2]]@rgb <- rgb2hex(col2rgb(border[1]))
                if(length(col)>0) pin@paths[[1]]@rgb <- rgb2hex(col2rgb(col[1]))
                if(length(col)>1) pin@paths[[3]]@rgb <- rgb2hex(col2rgb(col[2]))
            }
            switch(type,
                circle={
                    if(length(border)==0) border <- "black"
                    if(length(col)==0) col <- "white"
                    if(scoreType){
                        for(i in seq_len(this.score)){
                            y0 <- y2+y3+y4+2*radius*ratio.yx*(i-.5)
                            if(side) y0 <- 1 - y0
                            if(length(shape)==this.score){
                                ## multiple shape for each point
                                this_shape <- shape[i]
                            }else{
                                this_shape <- ifelse(length(shape)>0, shape[1], 'circle')
                            }
                            if(length(border)==this.score){
                                ## multiple border for each point
                                this_border <- border[i]
                            }else{
                                this_border <- border[1]
                            }
                            if(length(col)==this.score){
                                ## multiple color for each point
                                this_col <- col[i]
                            }else{
                                this_col <- col[1]
                            }
                            if(length(lwd)==this.score){
                                ## multiple alpha for each point
                                this_lwd <- lwd[i]
                            }else{
                                this_lwd <- lwd[1]
                            }
                            if(length(alpha)==this.score){
                                ## multiple alpha for each point
                                this_alpha <- alpha[i]
                            }else{
                                this_alpha <- alpha[1]
                            }
                            if(length(this_shape)!=1) this_shape <- 'circle'
                            this_gp <- gpar(col=this_border, fill=this_col, lwd=this_lwd, alpha=this_alpha)
                            switch(this_shape, #"circle", "square", "diamond", "triangle_point_up", "star", or "triangle_point_down"
                                    circle=grid.circle1(x=x2, y=y0,
                                                        r=radius*ratio.yx*cex, 
                                                        gp=this_gp),
                                    square=grid.square(x=x2, y=y0,
                                                        r=radius*ratio.yx*cex, 
                                                        gp=this_gp),
                                    diamond=grid.diamond(x=x2, y=y0,
                                                        r=radius*ratio.yx*cex, 
                                                        gp=this_gp),
                                    triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
                                                        r=radius*ratio.yx*cex, 
                                                        gp=this_gp),
                                    triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
                                                        r=radius*ratio.yx*cex, 
                                                        gp=this_gp),
                                    star=grid.star(x=x2, y=y0,
                                                    r=radius*ratio.yx*cex, 
                                                    gp=this_gp),
                                    grid.circle1(x=x2, y=y0,
                                                r=radius*ratio.yx*cex, 
                                                gp=this_gp))
                            
                        }
                    }else{
                        y0 <- y2+y3+y4+(this.score-.5)*2*radius*ratio.yx
                        if(side) y0 <- 1 - y0
                        switch(shape,
                                circle=grid.circle1(x=x2, y=y0,
                                                    r=radius*ratio.yx*cex, 
                                                    gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                square=grid.square(x=x2, y=y0,
                                                    r=radius*ratio.yx*cex, 
                                                    gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                diamond=grid.diamond(x=x2, y=y0,
                                                    r=radius*ratio.yx*cex, 
                                                    gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
                                                                        r=radius*ratio.yx*cex, 
                                                                        gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
                                                                            r=radius*ratio.yx*cex, 
                                                                            gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                star=grid.star(x=x2, y=y0,
                                                r=radius*ratio.yx*cex, 
                                                gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
                                grid.circle1(x=x2, y=y0,
                                            r=radius*ratio.yx*cex, 
                                            gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)))
                        if(!is.na(id$label)){
                            y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4
                            if(side) y0 <- 1 - y0
                            id$gp$cex <- .75*id$gp$cex
                            id$x <- x2
                            id$y <- y0
                            do.call(grid.text, id)
                        }
                    }
                    },
                pie={
                    y0 <- y2+y3+y4+radius*ratio.yx
                    if(side) y0 <- 1 - y0
                    grid.pie(x=x2, y=y0, 
                                radius = radius*cex, 
                                col = col, 
                                border = border, 
                                percent=percent,
                                edges=edges,
                                lwd=lwd, alpha=alpha)
                    },
                pie.stack={
                    y0 <- y2+y3+y4+(2*percent$stack.factor.order-1)*radius*ratio.yx
                    if(side) y0 <- 1 - y0
                    grid.pie(x=x2, 
                                y=y0, 
                                radius = radius*cex, 
                                col = col, 
                                border = border, 
                                percent=percent[, !colnames(percent) %in% 
                                                    c("stack.factor.order", 
                                                    "stack.factor.first")],
                                edges=edges,
                                lwd=lwd, alpha=alpha)
                    },
                pin={
                    y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
                    if(side) y0 <- 1 - y0
                    grid.picture(picture=pin, x=x2, 
                                    y=y0,
                                    width=2*radius*ratio.yx*cex,
                                    height=3*radius*ratio.yx*cex+y4)
                    if(!is.na(id$label)){
                        y0 <- y2+y3+(this.score-.25)*2*radius*ratio.yx+2*y4/3
                        id$gp$cex <- .5*id$gp$cex
                        id$x <- x2
                        id$y <- y0
                        do.call(grid.text, id)
                    }
                    },
                flag={
                    if(is.na(id$label)){
                    id$label <- " "
                    }
                    this.cex <- id$gp$cex
                    LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))*this.cex
                    y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
                    if(side) y0 <- 1 - y0
                    LINEW <- as.numeric(convertX(stringWidth(paste0("o", id$label, "u")), "npc"))*this.cex
                    LINEW <- LINEW * sign(cos(pi*id$rot/180))
                    LINEH0 <- LINEW*ratio.yx*tan(pi*id$rot/180)
                    grid.polygon(x=c(x2, x2+LINEW, x2+LINEW, x2),
                                y=c(y0, y0+LINEH0, y0+LINEH0+LINEH*1.25, y0+LINEH*1.25),
                                gp=gpar(fill=col, col=border, alpha=alpha))
                    id$x <- x2+LINEW*.5
                    id$y <- y0 + LINEH*.625+LINEH0*.5
                    do.call(grid.text, id)
                },
                grid.pie(x=x2, y=y2+y3+y4+radius*ratio.yx, 
                            radius = radius*cex, 
                            col = col, 
                            border = border, 
                            percent=percent,
                            edges=edges,
                            lwd=lwd, alpha=alpha))
    }

    grid.circle1 <- function(x = 0.5, y = 0.5, r = 0.5, 
                                default.units = "npc", name = NULL, 
                                gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        twopi <- 2 * pi
        ratio.yx <- getYXratio()
        t2xy <- function(t) {
            t2p <- twopi * t + pi/2
            list(x = r * cos(t2p)/ratio.yx, y = r * sin(t2p))
        }
        P <- t2xy(seq.int(0, 1, length.out = 100))
        invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                                gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
        }

        grid.square <- function(x = 0.5, y = 0.5, r = 0.5, 
                                default.units = "npc", name = NULL, 
                                gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        ratio.yx <- getYXratio()
        invisible(grid.rect(unit(x,"npc"), unit(y, "npc"), 
                            width = unit(r*2/ratio.yx, "npc"), 
                            height = unit(r*2, "npc"),
                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
    }

    grid.diamond <- function(x = 0.5, y = 0.5, r = 0.5, 
                                default.units = "npc", name = NULL, 
                                gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        ratio.yx <- getYXratio()
        P <- 
            list(x = c(0, r/ratio.yx, 0, -r/ratio.yx), 
                y = c(-r, 0, r, 0))
        invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                                gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
    }

    grid.triangle_point_up <- function(x = 0.5, y = 0.5, r = 0.5, 
                                default.units = "npc", name = NULL, 
                                gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        ratio.yx <- getYXratio()
        P <- 
            list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx), 
                y = c(-r, -r, r, -r))
        invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                                gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
    }

    grid.triangle_point_down <- function(x = 0.5, y = 0.5, r = 0.5, 
                                        default.units = "npc", name = NULL, 
                                        gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        ratio.yx <- getYXratio()
        P <- 
            list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx), 
                y = c(r, r, -r, r))
        invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                                gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
    }

    plotLegend <- function(legend, this.height, LINEH){
        ypos <- this.height
        pch <- 21
        if(length(legend)>0){
            if(is.list(legend)){
                thisLabels <- legend[["labels"]]
                if(is.null(legend$gp)){
                gp <- legend[!names(legend) %in% formalArgs(legendGrob)]
                class(gp) <- "gpar"
                legend$gp <- gp
                }else{
                gp <- legend$gp
                }
                if(is.null(gp$cex)) gp$cex <- 1
            }else{
            thisLabels <- names(legend)
            gp <- gpar(fill=legend, cex=1)
            legend <- list(
                labels = thisLabels,
                gp = gp
            )
            }
            if(length(thisLabels)>0){ 
            if(is.null(legend$byrow)) legend$byrow <- TRUE
            if(is.null(legend$vgap)) legend$vgap <- unit(.1*gp$cex[1], "lines")
            if(is.null(legend$hgap)) legend$hgap <- unit(.5*gp$cex[1], "lines")
            if(is.null(legend$ncol)){
                legend$ncol <- getColNum(thisLabels, cex=gp$cex)
            }
            ncol <- legend$ncol
            if(is.null(legend$pch)) legend$pch <- pch
            legend <- legend[names(legend) %in% formalArgs(legendGrob)]
            topblank <- ceiling(length(thisLabels) / ncol) * gp$cex[1]
            pushViewport(viewport(x=.5, 
                                    y=ypos+(topblank+.2*gp$cex[1])*LINEH/2, 
                                    width=1,
                                    height=topblank*LINEH,
                                    just="bottom"))
            this.height <- ypos + (topblank+.2*gp$cex[1])*LINEH
            do.call(grid.legend, legend)
            popViewport()
            }
        }
        this.height + LINEH
    }

    getHeight <- function(SNPs, ratio.yx, LINEW, GAP, cex, type, scoreMax,
                        level=c("data", "data&labels")){
        level=match.arg(level)
        stack.factors <- unique(as.character(SNPs$stack.factor))
        stack.factors <- sort(stack.factors)
        if(level=="data"){
            switch(type,
                circle={
                    labels.y <- LINEW + # add gaps for labels
                        6.5*GAP + 
                        scoreMax * LINEW * ratio.yx
                },
                pin={
                    if(length(SNPs$score)>0) {
                        this.scores <- ceiling(SNPs$score)
                    }else {
                        this.scores <- .5
                    }
                    this.scores[is.na(this.scores)] <- .5
                    labels.y <- LINEW + 
                        6.5*GAP + 
                        (this.scores-0.5*cex) * LINEW * ratio.yx
                },
                pie={
                    labels.y <- LINEW*max(ratio.yx, 1.2) + 
                        6.5*GAP + 0.5 * LINEW * ratio.yx * cex
                },
                pie.stack={
                    labels.y <- LINEW + 
                        6.5*GAP + 
                        (scoreMax-0.5*cex) * LINEW * ratio.yx
                },
                flag={
                    labels.y <- LINEW + 
                    6.5*GAP + 
                    scoreMax * LINEW * ratio.yx
                })
            labels.y
        }else{
            if(length(SNPs$label.parameter.rot)>0) {
                labels.rot <- SNPs$label.parameter.rot
            }else{
                labels.rot <- 90
            }
            labels.cex <- 1
            if(length(SNPs$label.parameter.gp)>0){
            if(length(SNPs$label.parameter.gp$cex)>0)
                labels.cex <- SNPs$label.parameter.gp$cex[[1]][1]
            }
            labels.length.rate <- labels.cex * max(cospi((labels.rot-90)/180), 0) * ratio.yx
            stringH <- as.numeric(convertY(stringHeight("W"), "npc"))
            
            switch(type,
                circle={
                    if(length(names(SNPs))>0){
                        maxStrHeight <- 
                            max(as.numeric(
                                convertX(stringWidth(names(SNPs)), "npc")
                            ))+stringH
                    }else{
                        maxStrHeight <- 0
                    }
                    maxStrHeight <- maxStrHeight * labels.length.rate
                    ypos <- LINEW + 6.5*GAP + 
                        scoreMax * LINEW * ratio.yx + maxStrHeight
                },
                pin={
                    if(length(names(SNPs))>0){
                        thisStrHeight <- max(as.numeric(
                            convertX(stringWidth(names(SNPs)), "npc")) ) +
                            stringH
                    }else{
                        thisStrHeight <- 0
                    }
                    thisStrHeight <- thisStrHeight * labels.length.rate
                    if(length(SNPs$score)>0){
                        ypos <- 
                            max(LINEW + 
                                    6.5*GAP + 
                                    (SNPs$score-0.5*cex) * LINEW * ratio.yx + 
                                    thisStrHeight)
                    }else{
                        ypos <- max(LINEW*max(ratio.yx, 1.2) + 
                                        6.5*GAP + thisStrHeight)
                    }
                },
                pie={
                    if(length(names(SNPs))>0){
                        maxStrHeight <- 
                            max(as.numeric(
                            convertX(stringWidth(names(SNPs)), "npc")
                            ))+stringH
                    }else{
                        maxStrHeight <- 0
                    }
                    maxStrHeight <- maxStrHeight * labels.length.rate
                    ypos <- LINEW + 
                        6.5*GAP + maxStrHeight
                },
                pie.stack={
                    if(length(names(SNPs))>0){
                        maxStrHeight <- 
                            max(as.numeric(
                            convertX(stringWidth(names(SNPs)), "npc")
                            ))+stringH
                    }else{
                        maxStrHeight <- 0
                    }
                    maxStrHeight <- maxStrHeight * labels.length.rate
                    ypos <- LINEW + 
                        6.5*GAP + maxStrHeight +
                        (scoreMax-0.5*cex) * LINEW * ratio.yx
                },
                flag={
                    if(length(names(SNPs))>0){
                    maxStrHeight <- 
                        max(as.numeric(
                        convertX(stringWidth(names(SNPs)), "npc")
                        ))+stringH
                    }else{
                    maxStrHeight <- 0
                    }
                    maxStrHeight <- maxStrHeight * labels.length.rate
                    ypos <- LINEW + 6.5*GAP + 
                    scoreMax * LINEW * ratio.yx + maxStrHeight
                }
            )
            ypos
        }
    }


    grid.star <- function(x = 0.5, y = 0.5, r = 0.5, 
                            default.units = "npc", name = NULL, 
                            gp = gpar(), draw = TRUE, vp = NULL){
        fill <- gp$fill
        col <- gp$col
        lwd <- if(length(gp$lwd)>0) gp$lwd else 1
        alpha <- gp$alpha
        if(is.null(fill)) fill <- "white"
        ratio.yx <- getYXratio()
        i <- 1:11
        angle <- 180
        alpha <- 2*pi / 10
        r <- r * (i %% 2 + 1)/2
        omega <- alpha * i + angle * pi /180
        invisible(grid.polygon(unit(r*sin(omega)/ratio.yx+x,"npc"), 
                                unit(r*cos(omega)+y, "npc"), 
                                gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
    }


    plotLollipops <- function(SNPs, feature.height, bottomHeight, baseline, 
                            type, ranges, yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType,
                            LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                            side=c("top", "bottom"), jitter=c("node", "label"),
                            main=TRUE){
        side <- match.arg(side)
        jitter <- match.arg(jitter)
        if(side=="top"){
            pushViewport(viewport(y=bottomHeight,
                                height=1,
                                just="bottom",
                                xscale=c(start(ranges), 
                                        end(ranges)),
                                clip="off"))
        }else{
            pushViewport(viewport(y=bottomHeight+feature.height,
                                height=1,
                                just="top",
                                xscale=c(start(ranges), 
                                        end(ranges)),
                                yscale=c(1, 0),
                                clip="off"))
        }
        
        if(type=="pie.stack" && length(SNPs$stack.factor)>0){
            stopifnot(is.vector(SNPs$stack.factor, mode="character"))
            if(length(SNPs$stack.factor.order)>0 || 
            length(SNPs$stack.factor.first)>0){
                warning("stack.factor.order and stack.factor.first are used by this function!",
                        "The values in these column will be removed.")
            }
            ## condense the SNPs
            stack.factors <- unique(as.character(SNPs$stack.factor))
            stack.factors <- sort(stack.factors)
            stack.factors.order <- seq_along(stack.factors)
            names(stack.factors.order) <- stack.factors
            SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs), 
                            as.character(SNPs$stack.factor))]
            SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
            SNPs$stack.factor.first <- !duplicated(SNPs)
            SNPs.condense <- SNPs
            SNPs.condense$oid <- seq_along(SNPs)
            SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), start(SNPs), end(SNPs))
            SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
            SNPs.condense <- lapply(SNPs.condense, function(.ele){
                .oid <- .ele$oid
                .gr <- .ele[1]
                mcols(.gr) <- NULL
                .gr$oid <- NumericList(.oid)
                .gr
            })
            SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
            SNPs.condense <- sort(SNPs.condense)
            lab.pos.condense <- start(SNPs.condense)
            label.parameter.draw <- rep(TRUE, length(SNPs.condense))
            if(length(SNPs$label.parameter.draw)==length(SNPs)){
            label.parameter.draw <- 
                SNPs$label.parameter.draw[vapply(SNPs.condense$oid, `[`, i=1,
                                                FUN.VALUE = numeric(1L))]
            label.parameter.draw <- vapply(label.parameter.draw, `[`, i=1,
                                            FUN.VALUE = logical(1L))
            }
            lab.pos.condense[label.parameter.draw] <-
            jitterLables(start(SNPs.condense)[label.parameter.draw], 
                                            xscale=c(start(ranges), end(ranges)), 
                                            lineW=LINEW*cex)
            lab.pos.condense[label.parameter.draw] <-
            reAdjustLabels(lab.pos.condense[label.parameter.draw],
                            lineW=LINEW*cex)
            condense.ids <- SNPs.condense$oid
            lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
            lab.pos <- lab.pos[order(unlist(condense.ids))]
        }else{
            lab.pos <- start(SNPs)
            label.parameter.draw <- getDrawLabelParam(SNPs)
            lab.pos[label.parameter.draw] <-
            jitterLables(start(SNPs)[label.parameter.draw], 
                        xscale=c(start(ranges), end(ranges)), 
                        lineW=LINEW*cex)
            lab.pos[label.parameter.draw] <-
            reAdjustLabels(lab.pos[label.parameter.draw], lineW=LINEW*cex)
        }
        if(length(SNPs)>0){
            yaxisat <- NULL
            yaxisLabel <- TRUE
            diameter <- LINEW*ratio.yx
            y2 <- feature.height
            y3 <- 4*GAP*cex
            y4 <- 2.5*GAP*cex
            if(length(yaxis)>1 && is.numeric(yaxis)){
                yaxisat <- yaxis
                if(length(names(yaxis))==length(yaxis)) yaxisLabel <- names(yaxis)
                yaxis <- TRUE
            }
            if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
                if(side=="top"){
                vp <- viewport(x=.5+ifelse(main, -1, 1) *LINEW,
                                y=feature.height + y3 + y4 + scoreMax*diameter/2,
                                width=1,
                                height=(scoreMax+1)*diameter,
                                yscale=c(0, scoreMax0+scoreMax0/scoreMax))
                }else{
                vp <- viewport(x=.5+ifelse(main, -1, 1) *LINEW,
                                y=1-(feature.height + y3 + y4 +
                                        scoreMax*diameter/2),
                                width=1,
                                height=(scoreMax+1)*diameter,
                                yscale=c(scoreMax0+scoreMax0/scoreMax, 0))
                }
            grid.yaxis(at=yaxisat,
                        label=yaxisLabel,
                        main = main,
                        gp=yaxis.gp,
                        vp=vp)
            }
            
            if(length(SNPs$alpha)==length(SNPs)){
            SNPs$alpha[is.na(SNPs$alpha)] <- 0
            if(all(is.numeric(SNPs$alpha))){
                if(any(SNPs$alpha>1)){## convert to 0-1
                SNPs$alpha <- SNPs$alpha/max(SNPs$alpha)
                }
            }else{ ## not correct format.
                SNPs$alpha <- as.numeric(factor(as.character(SNPs$alpha)))
                SNPs$alpha <- (SNPs$alpha+max(SNPs$alpha))/max(SNPs$alpha)/2
            }
            }else{
            SNPs$alpha <- NULL
            }
            if(type=="circle"){
            if(length(SNPs$shape)==length(SNPs)){
                ## shape could only be "circle", "square", "diamond", "triangle_point_up", "triangle_point_down"
                if(!all(unlist(SNPs$shape) %in% c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"))){
                message('shape must be "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"')
                if(is.list(SNPs$shape)){
                    stop("The shape is a list, please confirm the format.",
                        'It must be a list with elements of "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"')
                }
                SNPs$shape <- as.numeric(factor(SNPs$shape))
                SNPs$shape <- rep(c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"), 
                                    max(SNPs$shape))[SNPs$shape]
                }else{
                if(is.list(SNPs$shape)){
                    if(any(lengths(SNPs$shape)[SNPs$score!=0]==0)){
                    stop("The shape is a list, but zero length of shape is detected.")
                    }
                }
                }
                if(scoreType){
                if(is.list(SNPs$shape)){
                    if(!all(lengths(SNPs$shape)[SNPs$score!=0]==
                            SNPs$score[SNPs$score!=0])){
                    warning('Not all the length of shape equal to score.')
                    }
                }
                }
            }else{
                SNPs$shape <- NULL
            }
            }
            for(m in seq_along(SNPs)){
                this.dat <- SNPs[m]
                color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
                border <- 
                    if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
                fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
                alpha <- if(length(this.dat$alpha)>0) this.dat$alpha[[1]] else 1
                lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
                shape <- if(length(this.dat$shape)>0) this.dat$shape[[1]] else "circle"
                rot <- if(length(this.dat$label.rot)>0) this.dat$label.rot else 15
                this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
                this.dashline.col <- 
                if(length(this.dat$dashline.col)>0) this.dat$dashline.col[[1]][1] else dashline.col
                ## control plot the dash line or not
                if(length(names(this.dat))<1) this.dashline.col <- NA
                if(length(SNPs$label.parameter.label)==length(SNPs) && length(SNPs)>0){
                if(this.dat$label.parameter.label[[1]]=="" ||
                    is.na(this.dat$label.parameter.label[[1]])) this.dashline.col <- NA
                }
                if(length(SNPs$label.parameter.pfm)==length(SNPs) && length(SNPs)>0){
                if(is.null(SNPs$label.parameter.pfm[[m]])){
                    this.dashline.col <- NA
                }
                }
                if(length(SNPs$label.parameter.draw)==length(SNPs) && length(SNPs)>0){
                if(!(SNPs$label.parameter.draw[[m]])){
                    this.dashline.col <- NA
                }
                }
                id <- 
                handleLabelParams(this.dat, cex = this.cex, prefix = "node.label.",
                                    label = if(is.character(this.dat$label)) this.dat$label else
                                    if(is.character(this.dat$node.label)) this.dat$node.label else NA,
                                    rot = if(length(this.dat$label.rot)>0) this.dat$label.rot else ifelse(type=="flag", 15, 0),
                                    gp = gpar(cex=this.cex,
                                            col = if(length(this.dat$label.col)>0) this.dat$label.col else "black"),
                                    just = "centre",
                                    hjust = .5,
                                    vjust = .5)
                this.dat.mcols <- mcols(this.dat)
                this.dat.mcols <- cleanDataMcols(this.dat.mcols, type)

                grid.lollipop(x1=convertX(unit(start(this.dat), "native"), "npc", 
                                        valueOnly=TRUE),  
                            y1=baseline,
                            x2=convertX(unit(ifelse(jitter=="node", 
                                                    lab.pos[m], 
                                                    start(this.dat)), 
                                            "native"), "npc", valueOnly=TRUE), 
                            y2=y2,
                            y3=y3, y4=y4, 
                            radius=LINEW/2,
                            col=color,
                            border=border,
                            percent=this.dat.mcols,
                            edges=100,
                            type=type,
                            ratio.yx=ratio.yx,
                            pin=pin,
                            scoreMax=scoreMax * LINEW,
                            scoreType=scoreType,
                            id=id,
                            cex=this.cex, lwd=lwd, dashline.col=this.dashline.col,
                            side=side, alpha=alpha, shape=shape)

            }

            this.height <- getHeight(SNPs, 
                                    ratio.yx, LINEW, GAP, cex, type,
                                    scoreMax=scoreMax,
                                    level="data")
            labels.keep <- getDrawLabelParam(SNPs)
            SNPs <- SNPs[labels.keep]
            lab.pos <- lab.pos[labels.keep]
            if(length(names(SNPs))>0){
                if(type=="pie.stack"){
                    ## unique lab.pos and SNPs
                    idx <- !duplicated(names(SNPs))
                    lab.pos <- lab.pos[idx]
                    SNPs <- SNPs[idx]
                }
                this.label <- 
                handleLabelParams(SNPs, cex = cex,
                                    prefix = "label.parameter.",
                                    x = lab.pos,
                                    label = names(SNPs),
                                    just = ifelse(side=="top", "left", "right"),
                                    hjust = NULL,
                                    vjust = NULL,
                                    rot = 90,
                                    check.overlap = FALSE,
                                    default.units = "native",
                                    gp = gpar(cex=cex),
                                    pfm = NULL,
                                    font = "Helvetica-Bold",
                                    fontface = "bold",
                                    ic.scale = TRUE)
                if(jitter=="label"){
                ## add guide lines
                rased.height <- 4*GAP*cex
                guide.height <- 2.5*GAP*cex
                for(i in seq_along(SNPs)){
                    this.dashline.col <- 
                    if(length(SNPs[i]$dashline.col)>0) 
                        SNPs[i]$dashline.col[[1]][1] else 
                        dashline.col
                    if(length(this.label$label)<1 && length(this.label$pfm)<1){
                    next
                    }
                    if(length(this.label$label)==length(SNPs) && length(SNPs)>0){
                    if(is.na(this.label$label[i])){
                        next
                    }
                    if(this.label$label[i]==""){
                        next
                    }
                    }
                    if(length(this.label$pfm)==length(SNPs) && length(SNPs)>0){
                    if(is.null(this.label$pfm[[i]])){
                        next
                    }
                    }
                    grid.lines(x=c(start(SNPs[i]), this.label$x[i]), 
                            y=c(this.height+feature.height-cex*LINEW, 
                                this.height+feature.height+rased.height),
                            default.units = this.label$default.units,
                            gp=gpar(col=this.dashline.col, lty=3))
                    grid.lines(x=c(this.label$x[i], this.label$x[i]),
                            y=c(this.height+rased.height+feature.height,
                                this.height+rased.height+
                                    guide.height+feature.height),
                            default.units = this.label$default.units,
                            gp=gpar(col=this.dashline.col, lty=3))
                }
                ## add this height
                this.height <- this.height + rased.height + guide.height
                }
                if(length(this.label$pfm)>0){
                if(!requireNamespace("motifStack", quietly = TRUE)){
                    stop("When plot motifs as labels,",
                        " the Bioconductor package 'motifStack' is required!")
                }
                for(idx in seq_along(this.label$pfm)){
                    if(!is.null(this.label$pfm[[idx]])){
                    this_cex <- ifelse(length(cex)==length(this.label$pfm),
                                        cex[[idx]], cex[1])
                    this_just <- ifelse(length(this.label$just)==length(this.label$pfm),
                                        this.label$just[[idx]], this.label$just[1])
                    this_rot <- ifelse(length(this.label$rot)==length(this.label$pfm),
                                        this.label$rot[[idx]], this.label$rot[1])
                    this_font <- ifelse(length(this.label$font)==length(this.label$pfm),
                                        this.label$font[[idx]], this.label$font[1])
                    this_fontface <- ifelse(length(this.label$fontface)==length(this.label$pfm),
                                            this.label$fontface[[idx]], this.label$fontface[1])
                    this_ic.scale <- ifelse(length(this.label$ic.scale)==length(this.label$pfm),
                                            this.label$ic.scale[[idx]], this.label$ic.scale[1])
                    pushViewport(viewport(x=this.label$x[[idx]],
                                            y=this.height + feature.height,  
                                            just = this_just,
                                            width =  convertWidth(
                                            stringWidth(paste(rep("A",
                                                                    ncol(this.label$pfm[[idx]]@mat)),
                                                                collapse = "")), 
                                                        unitTo="npc",
                                                        valueOnly=FALSE),
                                            height = LINEW * this_cex,
                                            angle = this_rot,
                                            default.units = this.label$default.units))
                    motifStack::plotMotifLogoA(pfm = this.label$pfm[[idx]],
                                    font=this_font,
                                    fontface = this_fontface,
                                    ic.scale = this_ic.scale)
                    popViewport()
                    }
                }
                }else{
                if(any(is.na(this.label$label))){
                    ## grid.text will not plot if first element is empty
                    this.label <- filterThisLabel(this.label)
                }
                if(length(this.label$label)>0){
                    grid.text(x=this.label$x, y=this.height + feature.height, 
                            label = this.label$label,  
                            just = this.label$just, 
                            hjust = this.label$hjust,
                            vjust = this.label$vjust,
                            rot=this.label$rot,
                            check.overlap = this.label$check.overlap,
                            default.units = this.label$default.units,
                            gp=this.label$gp)
                }
                }
            }
        }
        popViewport()
    }

    reAdjustLabels <- function(coor, lineW){
        # resort
        coor <- sort(coor)
        bins <- ceiling(1/lineW)
        pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
        pos.bin <- cut(pos, c(-Inf, (0:bins)*lineW, Inf), labels=0:(bins+1), right=FALSE)
        
        ## split the coors by into clusters
        ## give the clusters with more idx more spaces if there are spaces between clusters
        tbl <- table(pos.bin)
        if(all(tbl<2)) return(coor)
        tbl.len <- length(tbl)
        if(tbl.len<3) return(coor)
        loops <- 1000
        loop <- 1
        while(any(tbl==0) && any(tbl>1) && loop < loops){
            tbl.bk <- tbl
            for(i in order(tbl.bk, decreasing=TRUE)){
            if(tbl[i]>1 && tbl.bk[i]==tbl[i]){
                if(i==1){
                if(tbl[2]<tbl[1]){
                    half <- sum(tbl[1:2])/2
                    tbl[2] <- ceiling(half)
                    tbl[1] <- floor(half)
                }
                }else{
                if(i==tbl.len){
                    if(tbl[tbl.len]>tbl[tbl.len-1]){
                    half <- sum(tbl[(tbl.len-1):tbl.len])/2
                    tbl[tbl.len-1] <- ceiling(half)
                    tbl[tbl.len] <- floor(half)
                    }
                }else{
                    if(tbl[i-1]<tbl[i+1]){
                    ## i-1 and i should be balanced
                    half <- sum(tbl[(i-1):i])/2
                    tbl[i-1] <- floor(half)
                    tbl[i] <- ceiling(half)
                    }else{
                    half <- sum(tbl[i:(i+1)])/2
                    tbl[i] <- floor(half)
                    tbl[i+1] <- ceiling(half)
                    }
                }
                }
            }
            }
            loop <- loop + 1
        }
        coef <- unlist(lapply(tbl, function(.ele){
            if(.ele==0) return(0)
            .ele <- seq(from=0, to=1, length.out=.ele+1)
            (.ele[-length(.ele)] + .ele[-1])/2
        }))
        coef <- coef[coef!=0]
        coor <- (rep(as.numeric(names(tbl)), tbl) - 1 + coef) * lineW
        coor <- convertX(unit(coor, "npc"), "native", valueOnly=TRUE)
        coor
    }

    plotFeatures <- function(feature.splited, LINEH, bottomHeight, 
                            label_on_feature=FALSE){
        feature.height <- 0
        for(n in seq_along(feature.splited)){
            this.feature.height <- 
                max(c(feature.splited[[n]]$height/2, 
                    .0001)) + 0.2 * LINEH
            feature.height <- feature.height + this.feature.height
            ##baseline
            grid.lines(x=c(0, 1), y=c(bottomHeight+feature.height, 
                                    bottomHeight+feature.height))
            for(m in seq_along(feature.splited[[n]])){
                this.dat <- feature.splited[[n]][m]
                color <- if(is.list(this.dat$color)) this.dat$color[[1]] else 
                    this.dat$color
                if(length(color)==0) color <- "black"
                fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else 
                    this.dat$fill
                if(length(fill)==0) fill <- "white"
                this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
                lwd <- if(length(this.dat$lwd)>0) this.dat$lwd[[1]][1] else 1
                this.feature.height.m <- 
                    if(length(this.dat$height)>0) 
                        this.dat$height[[1]][1] else 
                            2*this.feature.height
                grid.rect(x=start(this.dat)-.1, y=bottomHeight+feature.height, 
                        width=width(this.dat)-.8, 
                        height=this.feature.height.m,
                        just="left", gp=gpar(col=color, fill=fill, lwd=lwd), 
                        default.units = "native")
                if(!is.null(names(this.dat))){
                if(label_on_feature & !is.na(names(this.dat)[1])){
                    grid.text(x=(start(this.dat)+end(this.dat))/2, 
                            y=bottomHeight+feature.height,
                            just = "centre",
                            label = names(this.dat)[1],
                            gp= gpar(list(cex=this.cex * 
                                            this.feature.height.m/
                                            this.feature.height,
                                            color=color)), 
                            default.units = "native")
                }
                }
            }
            feature.height <- feature.height + this.feature.height
        }
        feature.height
    }

    cleanDataMcols <- function(this.dat.mcols, type){
        this.dat.mcols <- 
            this.dat.mcols[, 
                        !colnames(this.dat.mcols) %in% 
                            c("color", "fill", "lwd", "id", 
                            "cex", "dashline.col", 
                            "id.col", "stack.factor", "SNPsideID",
                            "shape", "alpha"), 
                        drop=FALSE]
        if(type!="pie.stack"){
            this.dat.mcols <- 
            this.dat.mcols[, !colnames(this.dat.mcols) %in% 
                            c("stack.factor.order", 
                                "stack.factor.first"), 
                            drop=FALSE]
        }
        this.dat.mcols <- 
            this.dat.mcols[, !grepl("^label.parameter",
                                    colnames(this.dat.mcols)), 
                        drop=FALSE]
        this.dat.mcols <- 
            this.dat.mcols[, !grepl("^node.label",
                                    colnames(this.dat.mcols)), 
                        drop=FALSE]
        return(this.dat.mcols)
    }

    plot_grid_xaxis <- function(xaxis, gp=gpar(col="black")){
        ## axis, should be in the bottom of transcripts
        if(length(xaxis)==1 && as.logical(xaxis)) {
            grid.xaxis(gp=gp)
        }
        if(length(xaxis)>1 && is.numeric(xaxis)){
            xaxisLabel <- names(xaxis)
            if(length(xaxisLabel)!=length(xaxis)) xaxisLabel <- TRUE
            grid.xaxis(at=xaxis, label=xaxisLabel, gp=gp)
        }
    }

    convertHeight2NPCnum <- function(.ele){
        if(is(.ele, "unit")){
            return(convertHeight(.ele, unitTo="npc", valueOnly=TRUE))
        }else{
            if(is.list(.ele)){
            .ele <- sapply(.ele, function(.e){
                if(is(.e, "unit")){
                .e <- convertHeight(.e, unitTo="npc", valueOnly=TRUE)
                }
                .e[1]
            })
            return(unlist(.ele))
            }else{
            if(is.numeric(.ele)){
                return(.ele)
            }else{
                if(is.integer(.ele)){
                return(.ele)
                }else{
                return(.ele)
                }
            }
            }
        }
    }

    setFeatureLayerID <- function(feature, range){
        feature <- feature[end(feature)>=start(range) & 
                            start(feature)<=end(range)]
        if(length(feature$featureLayerID)!=length(feature)){
            feature$featureLayerID <- rep("1", length(feature))
            feature$featureLayerID <- as.character(feature$featureLayerID)
            start(feature)[start(feature)<start(range)] <- start(range)
            end(feature)[end(feature)>end(range)] <- end(range)
        }
        return(feature)
    }

    plotFeatureLegend <- function(feature, LINEH, range, xaxis, xaxis.gp, label_on_feature=FALSE){
        if(length(xaxis)>1 || as.logical(xaxis[1])){
            xaxisSpace <- 2
            if(is.numeric(xaxis.gp$cex)) xaxisSpace <- 2*xaxis.gp$cex
        }else{
            xaxisSpace <- 0
        }
        if(length(names(feature))>0 & !label_on_feature ){ ## features legend
            feature.s <- feature[!duplicated(names(feature))]
            cex <- if(length(unlist(feature.s$cex))==length(feature.s)) 
            unlist(feature.s$cex) else 1
            ncol <- getColNum(names(feature.s), cex=cex)
            featureLegendSpace <- max(ceiling(length(names(feature.s)) / ncol) * cex + 1 )
            pushViewport(viewport(x=.5, y=featureLegendSpace*LINEH/2, 
                                width=1,
                                height=featureLegendSpace*LINEH,
                                xscale=c(start(range), end(range))))
            color <- if(length(unlist(feature.s$color))==length(feature.s)) 
            unlist(feature.s$color) else "black"
            fill <- if(length(unlist(feature.s$fill))==length(feature.s)) 
            unlist(feature.s$fill) else "black"
            pch <- if(length(unlist(feature.s$pch))==length(feature.s)) 
            unlist(feature.s$pch) else 22
            grid.legend(label=names(feature.s), ncol=ncol,
                        byrow=TRUE, vgap=unit(.2, "lines"),
                        hgap=unit(.5, "lines"),
                        pch=pch,
                        gp=gpar(col=color, fill=fill, cex=cex))
            popViewport()
        }else{
            featureLegendSpace <- 1
        }
        bottomblank <- (xaxisSpace + featureLegendSpace) * LINEH
        return(bottomblank)
    }

    maxStringWidth <- function(labels, spaces="WW", cex){
        max(as.numeric(convertX(stringWidth(paste0(labels, spaces)), "line"))*cex)
    }

    getYXratio <- function(){
        as.numeric(convertHeight(unit(1, 'snpc'), 'npc'))/
            as.numeric(convertWidth(unit(1, "snpc"), "npc"))
    }

    getColNum <- function(labels, spaces="WW", cex){
        ncol <- floor(as.numeric(convertX(unit(1, "npc"), "line")) / 
                        maxStringWidth(labels, spaces=spaces, cex) / 
                as.numeric(convertX(stringWidth("W"), "line")))
        nrow <- ceiling(length(labels) / ncol)
        ncol <- ceiling(length(labels) / nrow)
        ncol
    }

    ######################
    stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList", "list")))
    stopifnot(inherits(features, c("GRanges", "GRangesList", "list")))
    jitter <- match.arg(jitter)
    legendPosition <- checkLegendPosition(legendPosition)
    rescale.old <- rescale
    xaxis.old <- xaxis
    if(any(type!="circle"&jitter=="label")){
      jitter[which(type!="circle"&jitter=="label")] <- "node"
      warning("if jitter set to label, type must be cirle.")
      message("jitter is set to node.")
    }
    SNP.gr.name <- deparse(substitute(SNP.gr))
    if(is(SNP.gr, "GRanges")){
        SNP.gr <- list(SNP.gr)
        if(length(SNP.gr.name)==length(SNP.gr)) {
          names(SNP.gr) <- SNP.gr.name
        }
    }
    len <- length(SNP.gr)
    for(i in seq.int(len)){
        stopifnot(is(SNP.gr[[i]], "GRanges"))
    }
    
    
    ## prepare the feature
    if(inherits(features, c("GRangesList", "list"))){
      for(i in seq_along(features)){
        stopifnot("features must be a GRanges or GRangesList object"=
                    is(features[[i]], "GRanges"))
      }
      features <- features[seq.int(len)]
    }else{
      stopifnot("features must be a GRanges or GRangesList object"=
                  is(features, "GRanges"))
      #features.name <- deparse(substitute(features))
      features <- list(features)[seq.int(len)]
    }
    
    
    TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
    if(any(!type %in% TYPES)){
        stop("Error in match argument: ",
             paste0("'type' should be one of '",  
                    paste(TYPES, collapse="', '"), "'."))
    }
    types <- rep(type, length=len)[seq.int(len)]
    rm(type)
    ############### handle legend ####################
    ## set the legend as a list, 
    ## if all the legend for different tracks is same
    ## set draw legend for last track later
    legend <- handleLegend(legend, len, SNP.gr)
    
    ################ handle ranges #####################
    ## if !missing(ranges) set ranges as feature ranges
    ranges <- handleRanges(ranges, SNP.gr, features, len)
    
    ##cut all SNP.gr by the range
    SNP.gr <- cutSNP(SNP.gr, ranges, len)
    
    dots <- list(...)
    if(length(dots$maxYlimNchar)>0) {
      maxYlimNchar<- dots$maxYlimNchar
    } else {
      if(all(c(TRUE, vapply(yaxis, isTRUE, logical(1L))))){
        maxYlimNchar <- getMaxYlimNchar(SNP.gr, 3.5, types)
      } else {
        if(length(yaxis)>1){
          maxYlimNchar <- max(c(1, nchar(as.character(yaxis))), na.rm = TRUE) +
            4
        }else{
          maxYlimNchar<- 2.5
        }
      }
    }
    
    ################## plot ############################
    ## total height == 1
    height <- 1/sum(lengths(ranges))
    args <- as.list(match.call())
    if(length(args$height0)==0){
      height0 <- 0
    }else{
      height0 <- args$height0
    }
    if(newpage) grid.newpage()
    for(i in seq.int(len)){
      if(length(ranges[[i]])>1){## split into multilayers
        args$newpage <- FALSE
        for(j in rev(seq_along(ranges[[i]]))){
          args$ranges <- ranges[[i]][j]
          args$SNP.gr <- SNP.gr[i]
          args$features <- features[[i]]
          args$type <- types[i]
          args$legend <- legend[[i]]
          args$legendPosition <- legendPosition
          args$height0 <- height0
          args$maxYlimNchar <- maxYlimNchar
          height0 <- do.call(what = lolliplot, args = args)
        }
      }else{## is GRanges with length==1
        type <- match.arg(types[i], TYPES)
        if(type=="pin"){ ## read the pin shape file
          pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
          pin <- readPicture(pinpath)
        }else{
          pin <- NULL
        }
        ## Here we don't know the real height of each tracks
        vp0 <- viewport(x=.5, y=height0 + height*0.5,
                       width=1, height=height)
        pushViewport(vp0)
        LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
        LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
        totalH <- as.numeric(unit(1, "npc"))
        if(LINEH > totalH/20){
          LINEH <- totalH/20
        }
        
        ## GAP the gaps between any elements
        GAP <- .2 * LINEH
        
        SNPs <- SNP.gr[[i]]
        strand(SNPs) <- "*"
        SNPs <- sort(SNPs)
        
        feature <- features[[i]]
        
        ## rescale
        rescale <- rescale.old
        xaxis <- xaxis.old
        if(is.logical(rescale)[1]){
          if(rescale[1]){
            range.tile <- tile(ranges[[i]], n = 5)[[1]]
            if(all(width(range.tile)>2)){
              range.tile.cnt <- countOverlaps(range.tile, SNPs)
              feature.start <- feature.end <- feature
              end(feature.start) <- start(feature.start)
              start(feature.end) <- end(feature.end)
              range.tile.cnt2 <- countOverlaps(range.tile, unique(c(feature.start, feature.end)))
              range.tile.cnt <- range.tile.cnt + range.tile.cnt2
              range.width <- width(ranges[[i]])
              range.tile.width <- log2(range.tile.cnt + 1)
              range.tile.width <- range.tile.width/sum(range.tile.width)
              range.tile.width <- range.width * range.tile.width
              range.tile.width <- cumsum(range.tile.width)
              range.tile.width <- start(ranges[[i]]) + c(0, round(range.tile.width)-1)
              rescale <- data.frame(from.start=start(range.tile), from.end=end(range.tile),
                                    to.start=range.tile.width[-length(range.tile.width)],
                                    to.end=range.tile.width[-1])
              rescale$to.start[-1] <- rescale$to.start[-1] + 1
            }
          }
        }else{
          if(is.numeric(rescale)){## rescale is a percentage, convert it into from.start, from.end, to.start, to.end
            ## reset the scale for exons and introns
            feature.rd <- disjoin(c(feature, ranges[[i]]))
            feature.segment.points <- sort(unique(c(start(feature.rd), end(feature.rd))))
            ## filter the points by viewer port.
            feature.segment.points <- 
              feature.segment.points[feature.segment.points>=start(ranges[[i]]) &
                                       feature.segment.points<=end(ranges[[i]])]
            rescale <- rescale/sum(rescale, na.rm = TRUE)
            rescale[is.na(rescale)] <- 0.01
            if(length(rescale)==length(feature.segment.points)-1){
              rescale.ir <- IRanges(feature.segment.points[-length(feature.segment.points)]+1, 
                                    feature.segment.points[-1])
              start(rescale.ir)[1] <- start(rescale.ir)[1]-1
              rescale.ir.width <- sum(width(rescale.ir))
              rescale.ir.new.width <- cumsum(round(rescale.ir.width*rescale, digits = 0))
              rescale <- data.frame(from.start=start(rescale.ir), 
                                    from.end=end(rescale.ir),
                                    to.start=feature.segment.points[1] + 
                                      c(0, rescale.ir.new.width[-length(rescale.ir.new.width)]),
                                    to.end=feature.segment.points[1] + rescale.ir.new.width)
            }else{
              stop("The length of rescale is not as same as the number of segments (including features and non-features).")
            }
          }else{
            if(is.character(rescale)){
              if(grepl("exon|intron", rescale[1], ignore.case = TRUE)){
                ## reset the scale for exons and introns
                feature.rd <- disjoin(c(reduce(feature), ranges[[i]]), ignore.strand=TRUE)
                feature.rd <- subsetByOverlaps(feature.rd, ranges[[1]], ignore.strand=TRUE)
                feature.segment.exon <- subsetByOverlaps(feature.rd, feature, ignore.strand=TRUE)
                feature.segment.intron <- subsetByOverlaps(feature.rd, feature, invert=TRUE, ignore.strand=TRUE)
                feature.segment.exon$type <- rep("exon", length(feature.segment.exon))
                feature.segment.intron$type <- rep("intron", length(feature.segment.intron))
                feature.rd <- sort(c(feature.segment.exon, feature.segment.intron))
                feature.segment.points <- sort(unique(c(start(feature.segment.exon), end(feature.segment.exon))))
                ## filter the points by viewer port.
                feature.segment.points <- 
                  feature.segment.points[feature.segment.points>=start(ranges[[i]]) &
                                           feature.segment.points<=end(ranges[[i]])]
                ratio.to.full.range <- 9
                if(grepl("exon", rescale[1], ignore.case = TRUE)){#set intron width to 1/10
                  if(grepl("^exon_\\d+$", rescale[1], ignore.case = TRUE)){
                    ratio.to.full.range <- 
                      as.numeric(sub("^exon_(\\d+)$", "\\1", rescale[1], ignore.case = TRUE))
                    ratio.to.full.range <- ratio.to.full.range/(100-ratio.to.full.range)
                  }
                  width.full.range <- sum(width(feature.rd)[feature.rd$type=="exon"])
                  rescale <- width(feature.rd)
                  rescale[feature.rd$type=="intron"] <- 
                    rep(ceiling(width.full.range/ratio.to.full.range/sum(feature.rd$type=="intron")),
                        length(rescale[feature.rd$type=="intron"]))
                }else{#set exon width to 1/10
                  if(grepl("^intron_\\d+$", rescale[1], ignore.case = TRUE)){
                    ratio.to.full.range <- 
                      as.numeric(sub("^intron_(\\d+)$", "\\1", rescale[1], ignore.case = TRUE))
                    ratio.to.full.range <- ratio.to.full.range/(100-ratio.to.full.range)
                  }
                  width.full.range <- sum(width(feature.rd)[feature.rd$type=="intron"])
                  rescale <- width(feature.rd)
                  rescale[feature.rd$type=="exon"] <- 
                    rep(ceiling(width.full.range/ratio.to.full.range/sum(feature.rd$type=="exon")),
                        length(rescale[feature.rd$type=="exon"]))
                }
                rescale <- rescale/sum(rescale)
                if(length(rescale)==length(feature.segment.points)-1){
                  rescale.ir <- IRanges(feature.segment.points[-length(feature.segment.points)]+1, 
                                        feature.segment.points[-1])
                  start(rescale.ir)[1] <- start(rescale.ir)[1]-1
                  rescale.ir.width <- sum(width(rescale.ir))
                  rescale.ir.new.width <- cumsum(round(rescale.ir.width*rescale, digits = 0))
                  rescale <- data.frame(from.start=start(rescale.ir), 
                                        from.end=end(rescale.ir),
                                        to.start=feature.segment.points[1] + 
                                          c(0, rescale.ir.new.width[-length(rescale.ir.new.width)]),
                                        to.end=feature.segment.points[1] + rescale.ir.new.width)
                }else{
                  stop("Something wrong with the auto-scale setting. Please report the bug to https://github.com/jianhong/trackViewer/issues")
                }
              }
            }
          }
        }
        if(is.data.frame(rescale)){
          if(all(c("from.start", "from.end", "to.start", "to.end") %in% colnames(rescale))){
            ## check the from coverage the whole region.
            checkflank <- function(x){
              to <- IRanges(x$to.start, x$to.end)
              xol <- findOverlaps(to, drop.self=TRUE, 
                                  drop.redundant=TRUE,minoverlap=2L)
              if(length(xol)>1){
                stop("There is overlaps of the rescale region for 'to' columns.")
              }
              x <- IRanges(x$from.start, x$from.end)
              xol <- findOverlaps(x, drop.self=TRUE, 
                                  drop.redundant=TRUE,minoverlap=2L)
              if(length(xol)>1){
                stop("There is overlaps of the rescale region for 'from' columns.")
              }
              xgap <- gaps(x, start=min(start(x)), end=max(end(x)))
              if(length(xgap)>0){
                stop("There is gaps of the rescale region for 'from' columns.")
              }
            }
            checkflank(rescale)
            rescale.gr <- function(x){
              if(is(x, "GRanges")){
                x.start <- start(x)
                x.end <- end(x)
                y <- c(x.start, x.end)
                x.cut <- cut(y, breaks=c(rescale$from.start[1], rescale$from.end+1),
                             labels=seq.int(nrow(rescale)), right=FALSE)
                y <- mapply(function(a, b){
                  if(!is.na(b)) {
                    rescale(a, to=c(rescale$to.start[b], rescale$to.end[b]),
                            from=c(rescale$from.start[b], rescale$from.end[b]))
                  }else{
                    a
                  }
                }, y, as.numeric(as.character(x.cut)))
                y <- round(y)
                start(x) <- 1
                end(x) <- y[seq_along(x)+length(x)]
                start(x) <- y[seq_along(x)]
                x
              }else{
                x.cut <- cut(x, breaks=c(rescale$from.start[1], rescale$from.end+1),
                             labels=seq.int(nrow(rescale)), right=FALSE)
                y <- mapply(function(a, b){
                  if(!is.na(b)) {
                    rescale(a, to=c(rescale$to.start[b], rescale$to.end[b]),
                            from=c(rescale$from.start[b], rescale$from.end[b]))
                  }else{
                    a
                  }
                }, x, as.numeric(as.character(x.cut)))
                y <- round(y)
                y
              }
            }
            feature <- rescale.gr(feature)
            SNPs <- rescale.gr(SNPs)
            if(is.logical(xaxis)[1]){
              if(xaxis[1]){
                xaxis <- c(rescale$to.start[1], rescale$to.end)
                names(xaxis) <- c(rescale$from.start[1], rescale$from.end)
              }
            }else{
              xaxis.names <- names(xaxis)
              if(length(xaxis.names)!=length(xaxis)){
                xaxis.names <- as.character(xaxis)
              }
              xaxis <- rescale.gr(xaxis)
              names(xaxis) <- xaxis.names
            }
          }
        }
        
        ## top legend
        if(legendPosition$position=='bottom'){
          ## not supported yet.
          message('not suppport yet.')
        }
        ## convert height to npc number
        feature$height <- convertHeight2NPCnum(feature$height)
        ## multiple transcripts in one gene could be separated by featureLayerID
        feature <- setFeatureLayerID(feature, ranges[[i]])
        feature.splited <- split(feature, feature$featureLayerID)
        
        ## bottomblank, the transcripts legend height
        bottomblank <- plotFeatureLegend(feature, as.numeric(convertY(unit(1, "line"), "npc")),
                                         ranges[[i]], xaxis, xaxis.gp, label_on_feature)
        
        ## get the max score and scoreType
        if(length(SNPs$score)>0){
          SNPs$score <- sapply(SNPs$score, mean) ## fix the bug of score is NumericList
        }
        scoreMax0 <- scoreMax <- 
          if(length(SNPs$score)>0) ceiling(max(c(SNPs$score, 1), na.rm=TRUE)) else 1
        if(type=="pie.stack") scoreMax <- length(unique(SNPs$stack.factor))
        if(!type %in% c("pie", "pie.stack")){
          scoreType <- 
            if(length(SNPs$score)>0) all(floor(SNPs$score)==SNPs$score) else FALSE
          if(length(yaxis)>1 && is.numeric(yaxis)){
            if(length(names(yaxis))!=length(yaxis)){
              names(yaxis) <- yaxis
            }
            scoreMax0 <- max(yaxis, scoreMax0)
            scoreMax <- max(yaxis, scoreMax)
          }
          if(scoreMax>lollipop_style_switch_limit) {
            SNPs$score <- 10*SNPs$score/scoreMax
            scoreMax <- 10*scoreMax0/scoreMax
            scoreType <- FALSE
          }else{
            scoreMax <- scoreMax0
          }
        }else{
          scoreType <- FALSE
        }
        
        ## if the type is caterpillar, there are lollipop in both sides
        ## plot the bottom lollipops first. And push a new viewport
        
        IsCaterpillar <- length(SNPs$SNPsideID) > 0
        if(IsCaterpillar){
          if(any(is.na(SNPs$SNPsideID)) || 
             !all(SNPs$SNPsideID %in% c('top', 'bottom'))){
            warning("Not all SNPsideID is top or bottom")
            IsCaterpillar <- FALSE
          }
        }
        
        if(IsCaterpillar){
          SNPs.top <- SNPs[SNPs$SNPsideID=='top']
          SNPs.bottom <- SNPs[SNPs$SNPsideID=='bottom']
        }else{
          SNPs.top <- SNPs
          SNPs.bottom <- GRanges()
        }
        if(length(SNPs.bottom)<1) IsCaterpillar <- FALSE
        ## viewport of plot region
        if(!IsCaterpillar){
          bottomblank <- bottomblank
        }
        
        ## decide the legend position
        # right side keep 3 lines TODO: figure out the plot region to decide the right side margin
        plotRegionWidth <- unit(1 - (maxYlimNchar+3)*LINEW, 'npc')
        plotRegionX <- unit(maxYlimNchar*LINEW, 'npc')
        if(legendPosition$position %in% c('left', 'right')){
          if(length(legendPosition$width)>0){
            if(is.unit(legendPosition$width)){
              plotRegionWidth <- plotRegionWidth - legendPosition$width
            }else if(is.numeric(legendPosition$width)){
              if(legendPosition$width>0 && legendPosition$width<1){
                legendPosition$width <- unit(legendPosition$width, 'npc')
                plotRegionWidth <- plotRegionWidth - legendPosition$width
              }else{
                stop('legendPosition$width must be a number in the range of 0 ~ 1')
              }
            }else{
              stop('legendPosition$width must be a number of a unit.')
            }
          }else{
            legendPosition$width <- convertX(
              unit(max(legendNchar(legend[[i]]))+3, 'char'),
              'npc')
            plotRegionWidth <- plotRegionWidth - legendPosition$width
          }
          if(legendPosition$position=='left'){
            plotRegionX <- legendPosition$width + unit(2*LINEW, 'npc')
          }else{
            plotRegionX <- unit(maxYlimNchar*LINEW, 'npc')
          }
        }
        vp_track <- viewport(x=plotRegionX, y=bottomblank/2 + .5, 
                             width= plotRegionWidth,
                             height= 1 - bottomblank,
                             xscale=c(start(ranges[[i]]), end(ranges[[i]])),
                             clip="off",
                             just = 'left')
        pushViewport(vp_track)
        ratio.yx <- getYXratio()
        ## plot xaxis
        bottomHeight <- 0
        if(IsCaterpillar){
          ## total height == maxscore + extension + gap + labels
          bottomHeight <- getHeight(SNPs=SNPs.bottom, 
                                    ratio.yx=ratio.yx, 
                                    LINEW=LINEW, 
                                    GAP=GAP, 
                                    cex=cex, 
                                    type=type,
                                    scoreMax=scoreMax,
                                    level="data&labels")
          #bottomHeight <- bottomHeight/len
          vp <- viewport(y=bottomHeight, just="bottom",
                         xscale=c(start(ranges[[i]]), end(ranges[[i]])))
          pushViewport(vp)
          xaxis.gp$col <- "gray"
          plot_grid_xaxis(xaxis, gp=xaxis.gp)
          popViewport()
        }else{
          plot_grid_xaxis(xaxis, gp=xaxis.gp)
        }
        
        ## the baseline, the center of the first transcript
        baseline <- 
          max(c(feature.splited[[1]]$height/2, 
                .0001)) + 0.2 * LINEH
        baselineN <- 
          max(c(feature.splited[[length(feature.splited)]]$height/2, 
                .0001)) + 0.2 * LINEH
        
        ##plot features
        feature.height <- plotFeatures(feature.splited, LINEH, bottomHeight, label_on_feature)
        
        if(length(SNPs.bottom)>0){
          plotLollipops(SNPs.bottom, feature.height, bottomHeight, baselineN, 
                        type, ranges[[i]], yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType, 
                        LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                        side="bottom", jitter=jitter)
        }
        feature.height <- feature.height + 2*GAP
        if(length(SNPs.top)>0){
          plotLollipops(SNPs.top, feature.height, bottomHeight, baseline, 
                        type, ranges[[i]], yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType, 
                        LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                        side="top", jitter=jitter)
        }
        
        ## top legend
        this.height <- getHeight(SNPs.top, 
                                 ratio.yx, LINEW, GAP, cex, type,
                                 scoreMax=scoreMax,
                                 level="data&labels")
        this.height0 <- this.height <- this.height + bottomHeight + feature.height
        if(legendPosition$position=='top'){
          this.height <- plotLegend(legend[[i]], this.height, LINEH)
          if('alpha' %in% names(legend[[i]])){
            legend[[i]]$alpha <- NULL
            if('pch' %in% names(legend[[i]])){
              legend[[i]]$pch <- NA
            }
            plotLegend(legend[[i]], this.height0, LINEH)
          }
        }
        popViewport()## vp_track
        
        this.height <- bottomblank + 
          this.height * (1 - bottomblank)
        
        ## ylab
        vp_ylab <- viewport(x=ifelse(legendPosition$position!='left', 0, 
                                     convertX(legendPosition$width, 'npc',
                                              valueOnly = TRUE)),
                            y=this.height*0.5, 
                            width=plotRegionWidth +
                              unit(2*maxYlimNchar*LINEW, 'npc'),
                            height=this.height,
                            just='left')
        pushViewport(vp_ylab)
        if(is.logical(ylab)){
          if(ylab && length(names(SNP.gr))>0){
            grid.text(names(SNP.gr)[i], x = LINEW, 
                      y = .5, rot = 90, gp=ylab.gp)
          }
        }
        if(is.character(ylab)){
          if(length(ylab)==1) ylab <- rep(ylab, len)
          grid.text(ylab[i], x = LINEW,
                    y = .5, rot = 90, gp=ylab.gp)
        }
        popViewport()#vp_ylab
        
        ## left or right legend
        if(legendPosition$position %in% c('left', 'right')){
          legendRegionWidth <- legendPosition$width
          if(legendPosition$position=='left'){
            legendRegionX <- unit(0, 'npc')
          }else{
            legendRegionX <- unit(1, 'npc') - legendPosition$width
          }
          vp_legend <- viewport(x=legendRegionX, y=.5, 
                                width= legendRegionWidth,
                                height= 1,
                                just = 'left')
          pushViewport(vp_legend)
          plotLegend(legend[[i]], 0, LINEH)
          popViewport()## vp_legend
        }
        
        popViewport()#vp0
        height0 <-  height0 + this.height*height
      }
    }
    return(invisible(height0))
}


set.seed(23)

# Palettete ---------------------------------------------------------------------
my_pal <- c(paletteer_d("RColorBrewer::Set1"),paletteer_d("RColorBrewer::Set3"))
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#D5BAFF", 
                   "#FFC8E1", "#C2F0FC","#C5FAD5", "#FFDAC1", "#E2F0CB", "#C3B1E1", "#F8B195", 
                   "#F67280","#F6E2B3", "#B5EAD7","#FFABAB", "#D0F4DE", "#A0CED9", "#FFCBF2")


# Load Data ---------------------------------------------------------------------
## Variants
variants_table <- readRDS(clean_table)

## Protein gff
protein_info <- rtracklayer::import(prot_file)

## Exon tsv
exon_info <- read_tsv(exon_file)

## Participant Metadata
metadata_info <- read_tsv(p_metadata_file)


# Visual Plots -----------------------------------------------------------------
## Distribution of MAF_variants ------
pdf(paste0(gene_name, "_MAF_Distribution.pdf"))
ggplot(variants_table, aes(x = log10(MAF_variant + 1e-6))) +
      geom_histogram(bins = 50, fill = "black") +
          labs(title = paste0(gene_name, " MAF Distribution"), x = "MAF (log10) + 1e-6",
                     y = "Number of variants") +
      theme_bw()
dev.off()


## MAF vs IMPACT -------
pdf(paste0(gene_name,"_MAFvsIMPACT.pdf"), width=14, height=8)
ggplot(variants_table, aes(x = log10(MAF_variant + 1e-6), y = IMPACT_annotation, 
                  color = CLIN_SIG_annotation, shape = SYMBOL_annotation)) +
      geom_jitter(width = 0.1, height = 0.1, size = 3, alpha = 0.8) +
      scale_x_continuous(name = "log10(MAF + 1e-6)") +
      scale_y_discrete(name = "IMPACT") +
      scale_shape_manual(values = c(16, 17)) +  
      theme_minimal(base_size = 14) +
      theme(legend.position = "right",plot.title = element_text(face="bold", size=16)) +
      labs(title = "MAF vs. IMPACT", color = "ClinVar",shape = "Gene")
dev.off()


## IMPACT Frequency -------
pdf(paste0(gene_name,"_IMPACT_Frequency.pdf"))
print(ggplot(variants_table, aes(x = IMPACT_annotation)) + 
      geom_bar(fill = "#0072B2", color = "black", width = 0.7) +
      geom_text(stat="count", aes(label= after_stat(count)),vjust=-0.5, color="black", size=5) +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Frequency of ", gene_name, " by IMPACT category"),
        x = "IMPACT Category", y = "Number of variants") +
      theme(plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1)))
dev.off()

## Density Variants -------
# Add a check to ensure there are at least two points in each group before calculating density
valid_groups <- variants_table %>% 
  group_by(IMPACT_annotation) %>% 
  filter(n() > 1)  # Filter groups with more than one point

if (nrow(valid_groups) > 0) {
  densities <- valid_groups %>% 
    group_by(IMPACT_annotation) %>% 
    group_map(~ density(.x$POS_variant)$y)

  # Filter out NA values from densities
  max_density <- max(unlist(densities), na.rm = TRUE)
  if (is.infinite(max_density)) {
    max_density <- 0  # Set to 0 if no valid densities are found
  }
} else {
  cat("Warning: Not enough data points for density calculation. Skipping.")
  max_density <- 0
}

# Calculate gene length and adaptive bandwidth
gene_start <- min(c(variants_table$POS_variant, exon_info$ExonStart), na.rm = TRUE)
gene_end <- max(c(variants_table$POS_variant, exon_info$ExonEnd), na.rm = TRUE)
gene_length <- gene_end - gene_start

adaptive_bw <- case_when(  # Define adaptive bandwidth based on gene length
  gene_length < 10000 ~ max(gene_length / 20, 100),    # For genes <10kb, bw = gene_length/20, minimum 100bp
  gene_length < 100000 ~ gene_length / 50,             # For genes 10kb-100kb, bw = gene_length/50
  TRUE ~ gene_length / 100                             # For genes >100kb, bw = gene_length/100
)

adaptive_bw <- round(adaptive_bw / 100) * 100

cat(paste0("Gene length: ", format(gene_length, big.mark = ","), " bp\n"))
cat(paste0("Adaptive bandwidth: ", format(adaptive_bw, big.mark = ","), " bp\n"))

pdf(paste0(gene_name,"_Density_variants.pdf"), width=14, height=8)
ggplot(variants_table, aes(x = POS_variant, color = IMPACT_annotation, fill = IMPACT_annotation)) +
    geom_density(bw=adaptive_bw,alpha = 0.3) +
    geom_rect(data = exon_info,
              aes(xmin = ExonStart, xmax = ExonEnd, ymin = -max_density*0.03, ymax = max_density*0.03),
              fill = "#00441B", alpha = 1, inherit.aes = FALSE) +
    geom_text(data = exon_info,
              aes(x = (ExonStart + ExonEnd) / 2, y = -max_density*0.03 * 3, label = seq_len(nrow(exon_info))),
              size = 3.5, inherit.aes = FALSE) +
    scale_color_manual(values =paletteer_d("RColorBrewer::RdYlGn")[c(1,3,5,10)]) +
    scale_fill_manual(values = paletteer_d("RColorBrewer::RdYlGn")[c(1,3,5,10)]) +
    facet_wrap(~ CHROM_variant, scales = "free_x") +
    labs(title = paste0("Density distribution of ",gene_name ," variants in genomic positions by impact annotation"),
        x = "Genomic Position (exons in dark green)", y = paste0("Density bw = ",adaptive_bw," bp"), color = "Impact", fill = "Impact") +
    theme_minimal()
dev.off()

    
# Lollipop Plot Function -----------------------------------------------------
# Function to create trackViewer lollipop plot
create_lollipop_plot <- function(df_data, gene_name, subset_label, protein_info, metadata_info) {
  
  # Prepare features
  feature2plot <- c("Domain", "Motif", "Region", "Zinc finger")
  features <- protein_info[protein_info$type %in% feature2plot, ]
  
  if(length(features) == 0) {
    cat("No protein features available for Lollipop plot\n")
    return()
  }
  
  # Prepare colors and labels for features
  features$Note <- trimws(features$Note)
  note_colors <- setNames(my_pal[1:length(unique(features$Note))], unique(features$Note))
  features$fill <- note_colors[features$Note]
  names(features) <- features$Note
  
  # Create variants GRanges
  sample.gr <- GRanges(seqnames = paste0(unique(seqnames(features))), 
                       ranges = IRanges(start = df_data$Protein_pos_start, width = 1, 
                                       names = df_data$LabelVarPlot))
  
  # Add participant Info 
  df_data <- df_data %>%
    rowwise() %>%
    mutate(samples_all = paste(Het_samples, Hom_samples, Hemi_samples, sep = ","),
           samples_all = gsub(",+", ",", samples_all),     
           samples_all = gsub("^,|,$", "", samples_all),   
           samples_all = trimws(samples_all),
           samples_vec = list(unique(unlist(strsplit(samples_all, ",")))),
           participants_vec = list(
             unique(na.omit(metadata_info$participant_id[match(samples_vec, metadata_info$plate_key)]))),
           participants_ID = paste(unique(participants_vec), collapse = ","),
           participants_num = length(unique(participants_vec))
    ) %>%
    ungroup()  %>%
    select(-samples_vec, -participants_vec)
  
  # Prepare variant properties
  sample.gr$score <- df_data$participants_num  
  sample.gr$Consequence_annotation <- df_data$Consequence_annotation
  sample.gr$color <- setNames(pastel_colors, unique(sample.gr$Consequence_annotation))[sample.gr$Consequence_annotation]
  legends <- list(labels=unique(sample.gr$Consequence_annotation), fill=unique(sample.gr$color))
  sample.gr$shape <- ifelse(grepl("^(rs|COSV)", names(sample.gr)), "circle", "diamond")
  
  # Prepare labels for plotting
  sample.gr.rot <- sample.gr
  sample.gr.rot$label.parameter.gp <- gpar(fontsize=7)
  sample.gr.rot$label.parameter.label <- sub(",.*", "", names(sample.gr))
  sample.gr.rot$label.parameter.label[which(sample.gr.rot$score < 2)] <- NA
  #sample.gr.rot$label.parameter.rot <- 45
  
  
  # Create filename and plot
  filename <- paste0(gene_name, "_Lollipop_", subset_label, ".pdf")

  pdf(filename, height=9)
  lolliplot(sample.gr.rot, features, legend=legends, ylab="Num. of Participants",
            yaxis.gp = gpar(fontsize=15), xaxis.gp = gpar(fontsize=15))
  grid.text(paste0(subset_label, " Variants of ", gene_name), 
                   x=.5, y=.98, gp=gpar(cex=1.5, fontface="bold"))
  dev.off()
  
  cat("Generated:", filename, "\n")
}

# Configuration for variant subsets -------------------------------------------
# Define which combinations to analyze
variant_subsets <- list(
  # High priority combinations
  list(
    name = "HIGH_All",
    filters = list(IMPACT_annotation = "HIGH"),
    description = "High_Impact_All"
  ),
  list(
    name = "HIGH_Canonical",
    filters = list(IMPACT_annotation = "HIGH", CANONICAL_annotation = "YES"),
    description = "High_Impact_Canonical"
  ),
  list(
    name = "HIGH_NonCanonical", 
    filters = list(IMPACT_annotation = "HIGH"),
    special_filters = list(CANONICAL_annotation = "not_YES"),  # This will include NA and other values
    description = "High_Impact_NonCanonical"
  ),
  list(
    name = "MODERATE_All",
    filters = list(IMPACT_annotation = "MODERATE"),
    description = "Moderate_Impact_All"
  ),
  list(
    name = "MODERATE_Canonical",
    filters = list(IMPACT_annotation = "MODERATE", CANONICAL_annotation = "YES"), 
    description = "Moderate_Impact_Canonical"
  )#,
  # Lower priority combinations (uncomment if needed)
  # list(
  #   name = "MODERATE_NonCanonical",
  #   filters = list(IMPACT_annotation = "MODERATE"),
  #   special_filters = list(CANONICAL_annotation = "not_YES"),
  #   description = "Moderate_Impact_NonCanonical"
  # ),
  # list(
  #   name = "LOW_All", 
  #   filters = list(IMPACT_annotation = "LOW"),
  #   description = "Low_Impact_All"
  # )
)

# Generate Lollipop Plots ------------------------------------------------------
cat("Generating lollipop plots for prioritized variant subsets...\n")

for (subset_config in variant_subsets) {
  
  # Apply standard filters
  df_subset <- variants_table
  for (column in names(subset_config$filters)) {
    if (column %in% colnames(variants_table)) {
      df_subset <- df_subset %>% 
        filter(!!sym(column) == subset_config$filters[[column]])
    } else {
      cat("Warning: Column", column, "not found in data. Skipping subset", subset_config$name, "\n")
      next
    }
  }
  
  # Apply special filters (for more complex filtering like "not YES")
  if (!is.null(subset_config$special_filters)) {
    for (column in names(subset_config$special_filters)) {
      if (column %in% colnames(variants_table)) {
        filter_value <- subset_config$special_filters[[column]]
        if (filter_value == "not_YES") {
          # Include everything except "YES" (includes NA, "NO", etc.)
          df_subset <- df_subset %>% 
            filter(is.na(!!sym(column)) | !!sym(column) != "YES")
        }
        # Add more special filter types here if needed
      } else {
        cat("Warning: Column", column, "not found in data. Skipping special filter\n")
      }
    }
  }
  
  # Additional base filters for lollipop plots
  df_subset <- df_subset %>% 
    filter(!is.na(Protein_pos_start), !is.na(NS_variant), NS_variant > 0)
  
  # Check if we have data to plot
  if (nrow(df_subset) == 0) {
    cat("No variants found for subset:", subset_config$description, "\n")
    next
  }
  

  # Generate lollipop plot
  cat("Processing subset:", subset_config$description, "(", nrow(df_subset), "variants )\n")
  create_lollipop_plot(df_subset, gene_name, subset_config$description, protein_info, metadata_info)
}

cat("Plots generated successfully for gene:", gene_name, "\n")
