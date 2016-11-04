#!/usr/bin/env Rscript
#
# Program: avgProf-Y.r
# Purpose: Plot average profiles across all genes from a given genomic region window
#
#          The generated avgProf plot is the inverse of that created by ngs.plot.r
#          which creates average profiles across the studied region for all genes
#
# Author:  Andy Saurin (andrew.saurin@univ-amu.fr)
# Created: 27 Sept 2014.
#

library(Biobase, quietly=TRUE, verbose=FALSE, warn.conflicts = FALSE)

cmd.help <- function() {
    cat("\n")
    cat("## Mandatory parameters:\n")
    cat("    -I   Result zip file created by ngs.plot\n")
    cat("    -O   Output name\n")
    cat("\n")
    cat("### Avg. profiles parameters:\n")
    cat("    -GR The genomic region as (min_val,max_val) pair to use to calculate average profiles, from 1-100 (default=10,30) (overridden by a non-zero value of -SW)\n")
    cat("    -SW Scan across the genomic region (1-100) using this window size (integer from 1-100). (default=0 (disabled))\n")
    cat("        If provided, the script will output 100/n PDF files where n is the integer provided.\n")
    cat("    -RM The running mean window size. (default=0 (auto))\n")
    cat("    -SE The amount of extra smoothing to apply to the plot (none=0, v.light=1, light=2 (default), medium=3, heavy=4)\n")
    cat("\n")
    cat("### Plot parameters:\n")
    cat("    -WD Image width(default=20)\n")
    cat("    -HG Image height(default=8)\n")
    cat("    -YAS Y-axis scale: auto(default) or min_val,max_val(custom scale)\n")
    cat("    -LEG Draw legend? 1(default) or 0\n")
    cat("    -BOX Draw box around plot? 1 or 0(default)\n")
    cat("    -XL Draw X-axis and label? 1(default) or 0\n")
    cat("    -LWD Line width(default=5)\n")
    cat("    -ZL When using BAM pairs, draw an abline through 0 on the Y-axis? 1(default) or 0 \n")
    cat("\n")
    cat("### Misc parameters:\n")
    cat("    -OM Output the matrix of running means as tab-delimited file? 1=yes 0=no (default=0)\n")
    cat("\n")
    cat("\n")
}

# Parse command line arguments and extract them into an associate array.
# Check if the required arguments are all satisfied.

parseArgs <- function(args, manditories) {
    if(length(args) %% 2 == 1 || length(args) == 0) {
        cat('Unpaired argument and value.\n')
        return(NULL)
    }
    n.i <- seq(1, length(args), by=2)
    v.i <- seq(2, length(args), by=2)
    args.name <- args[n.i]
    args.value <- args[v.i]

    # Check if required argument values are supplied.
    miss_tag <- F
    man.bool <- manditories %in% args.name
    if(!all(man.bool)){
        cat(paste('Missing argument: ', paste(manditories[!man.bool],
                                              collapse=','), '.', sep='')
           )
        miss_tag <- T
    }
    if(miss_tag){
        res <- NULL
    }else{
        res <- args.value
        names(res) <- args.name
    }
    res
}

PlotVars <- function(args.tbl, existing.vl=vector('character')) {

	updated.vl <- list()

	#### Genomic Region to scan ####
    if('-GR' %in% names(args.tbl)) {
        gr <- args.tbl['-GR']
        if(gr != 'auto') {
            gr <- unlist(strsplit(gr, ','))
            if(length(gr) == 2) {
                s.col <- as.numeric(gr[1])
                e.col <- as.numeric(gr[2])
                stopifnot(s.col < e.col)
            } else {
                stop("-GR must be a pair of numerics separated by ','\n")
            }
            updated.vl$s.col <- s.col
            updated.vl$e.col <- e.col
        } else {
            updated.vl$s.col <- 10
            updated.vl$e.col <- 30

        }
    } else if(!'s.col' %in% existing.vl) {
            updated.vl$s.col <- 10
            updated.vl$e.col <- 30
    }


    #### Scanning window. ####
    if('-SW' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-SW']) > 0 )
        if ( as.integer(args.tbl['-SW']) > 100 ) {
        	stop("-SW cannot be greater than 100")
        }
        updated.vl$winsize <- as.integer(args.tbl['-SW'])
    } else if(!'winsize' %in% existing.vl) {
        updated.vl$winsize <- 0
    }

    #### AvgProf Running Mean window. ####
    if('-RM' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-RM']) >= 0  )
        updated.vl$running.mean <- as.integer(args.tbl['-RM'])
    } else if(!'running.mean' %in% existing.vl) {
        updated.vl$running.mean <- 0
	}

    #### AvgProf smoothing. ####
    if('-SE' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-SE']) >= 0  )
        if ( as.integer(args.tbl['-SE']) > 4 ) {
        	stop("-SW cannot be greater than 3 (options: 0 (none), 1 (v.light), 2 (light=default), 3 (medium), 4 (heavy)")
        }
        updated.vl$smooth.extent <- as.integer(args.tbl['-SE'])
    } else if(!'smooth.extent' %in% existing.vl) {
        updated.vl$smooth.extent <- 2
    }


    #### Plot width. ####
    if('-WD' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-WD']) > 0)
        updated.vl$plot.width <- as.integer(args.tbl['-WD'])
    } else if(!'plot.width' %in% existing.vl) {
        updated.vl$plot.width <- 20
    }

    #### Plot height. ####
    if('-HG' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-HG']) > 0)
        updated.vl$plot.height <- as.integer(args.tbl['-HG'])
    } else if(!'plot.height' %in% existing.vl) {
        updated.vl$plot.height <- 8
    }

    #### Plot line width. ####
    if('-LWD' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LWD']) > 0)
        updated.vl$lwd <- as.integer(args.tbl['-LWD'])
    } else if(!'lwd' %in% existing.vl) {
        updated.vl$lwd <- 5
    }

	#### Y-axis scale ####
    if('-YAS' %in% names(args.tbl)) {
        ystr <- args.tbl['-YAS']
        if(ystr != 'auto') {
            yp <- unlist(strsplit(ystr, ','))
            if(length(yp) == 2) {
                y.min <- as.numeric(yp[1])
                y.max <- as.numeric(yp[2])
                stopifnot(y.min < y.max)
            } else {
                stop("-YAS must be 'auto' or a pair of numerics separated by ','\n")
            }
            updated.vl$yscale <- c(y.min, y.max)
        } else {
            updated.vl$yscale <- 'auto'
        }
    } else if(!'yscale' %in% existing.vl) {
        updated.vl$yscale <- 'auto'
    }

    #### Plot legend. ####
    if('-LEG' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LEG']) != 0 || as.integer(args.tbl['-LEG']) != 1)
        updated.vl$legend <- as.integer(args.tbl['-LEG'])
    } else if(!'legend' %in% existing.vl) {
        updated.vl$legend <- 1
    }

    #### Plot box. ####
    if('-BOX' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-BOX']) != 0 || as.integer(args.tbl['-BOX']) != 1)
        updated.vl$box <- ifelse( (args.tbl['-BOX']==1), TRUE, FALSE)
    } else if(!'box' %in% existing.vl) {
        updated.vl$box <- FALSE
    }

    #### X-axis label. ####
    if('-XL' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-XL']) != 0 || as.integer(args.tbl['-XL']) != 1)
        updated.vl$xlab <- ifelse( (args.tbl['-XL']==1), TRUE, FALSE)
    } else if(!'xlab' %in% existing.vl) {
        updated.vl$xlab <- TRUE
    }

    #### Draw Y-axis 0 abline. ####
    if('-ZL' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-ZL']) != 0 || as.integer(args.tbl['-ZL']) != 1)
        updated.vl$abline.y <- ifelse( (args.tbl['-ZL']==1), TRUE, FALSE)
    } else if(!'abline.y' %in% existing.vl) {
        updated.vl$abline.y <- TRUE
    }

    #### Output running mean matrix. ####
    if('-OM' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-OM']) != 0 || as.integer(args.tbl['-OM']) != 1)
        updated.vl$matrix.out <- as.integer(args.tbl['-OM'])
    } else if(!'matrix.out' %in% existing.vl) {
        updated.vl$matrix.out <- 0
    }

	return(updated.vl)
}

# Read command arguments.
args <- commandArgs(T)
if(length(args) < 1) {
    cmd.help()
    stop("No arguments provided.\n")
}


# Parse arguments.
args.tbl <- parseArgs(args, c('-I', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}
iname <- args.tbl['-I']  # iname: input zip file name
oname <- args.tbl['-O']  # oname: output file root name


# Load plotting parameters and data; update plotting parameters.
ifolder <- sub('.zip$', '', basename(iname))
if(ifolder == basename(iname)) {
    stop("Input filename must end with .zip\n")
}
load(unz(iname, file.path(ifolder, 'heatmap.RData')))
updated.vl <- PlotVars(args.tbl, ls())

# Load required libraries.
cat("Loading R libraries")
if(!suppressMessages(require(caTools, warn.conflicts=F))) {
    install.packages('caTools')
    if(!suppressMessages(require(caTools, warn.conflicts=F))) {
        stop('Loading package caTools failed!')
    }
}
cat('.')
if(!suppressMessages(require(matrixStats, warn.conflicts=F))) {
    install.packages('matrixStats')
    if(!suppressMessages(require(matrixStats, warn.conflicts=F))) {
        stop('Loading package matrixStats failed!')
    }
}
cat('.')
if(!suppressMessages(require(iterators, warn.conflicts=F))) {
    install.packages('iterators')
    if(!suppressMessages(require(iterators, warn.conflicts=F))) {
        stop('Loading package matrixStats failed!')
    }
}
cat('.')
if(!suppressMessages(require(foreach, warn.conflicts=F))) {
    install.packages('foreach')
    if(!suppressMessages(require(foreach, warn.conflicts=F))) {
        stop('Loading package foreach failed!')
    }
}
cat('Done\n')


smooth.mean <- function (x, k=10, iter=10, alg="C", endrule="mean" ) {
	# Calculate the smoothed mean of values using a running window
	#
	# Args:
	#   x: the data matrix
	#   k: running mean window size
	#   iter: how many iteractions of the function to run
	#   alg: running mean algorithm to use
	#   endrule: running mean endrule to use

	i=0
#	k=k*1.5
	repeat {
		if (i == iter) {
			break
		}
		x <- as.matrix(runmean(x, k=k, alg=alg, endrule=endrule))
		i <- i+1
		#if ( i%%k ) { cat(".") }
	}
	return(x)
}

plotmat <- function(regcovMat, title2plot=NULL, yscale, legend=1, box=FALSE, lwd=5, xlable=TRUE, bam.pair=FALSE, abline.y=TRUE) {
    # Plot the average profiles
    #
    # Args:
    #	regvoMat: the smoothed scores matrix
    #	title2plot: The plot lables
    #	yscale: The Y-axis scale
    #   legend: plot the lables (1) or not (0)
    #   box: box around the plot? (bool)
    #   lwd: The line width of the plots
    #   xlable: plot the X-axis scale? (bool)
    #   bam.pair: is this a bam pair enrichment profile? (bool)
    #   abline.y: draw an abline at 0 when using bam pairs? (bool)

    # Choose colors.
    ncurve <- ncol(regcovMat)
    if(ncurve <= 8) {
        suppressMessages(require(RColorBrewer, warn.conflicts=F))
        col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
        col2use <- col2use[1:ncurve]
    } else {
        col2use <- rainbow(ncurve)
    }
    #col2use <- col2alpha(col2use, 0.8)

	if ( length(yscale) != 2 ) {
		max <- ceiling(max(regcovMat))
		min <- floor(min(regcovMat))
		yscale <- c(min,max)
	}

	xrange <- matrix(seq(1:nrow(regcovMat)))

	par(cex=2, mar=c(5, 5, 4, 2))

	if ( xlable == TRUE ) {
		xaxt <- 's'
		xlab <- "Gene Number"
	} else {
		xaxt <- 'n'
		xlab <- NULL
	}

	xscale <- c(1,nrow(regcovMat))

	if ( bam.pair == TRUE ) {
		ylab <- "log2(Fold Change)"
	} else {
		ylab <- "RPM"
	}
    matplot(xrange, regcovMat,
            xaxt=xaxt, type="l", col=col2use, ylim=yscale, xlim=xscale,
            lty="solid", lwd=lwd, frame.plot=box, ann=F, xaxs="i",  yaxs="i")
    axis(2,lwd=lwd, lwd.ticks=lwd)
    if ( xlable == TRUE ) {
    	axis(1,lwd=lwd, lwd.ticks=lwd, at=xscale)
    }
    title(ylab=ylab, xlab=xlab)

    if ( bam.pair == TRUE && abline.y == TRUE ) {
	    abline(0,0, col="grey", lwd=lwd, lty=2)
    }

    if ( !is.null(title2plot) && legend == 1) {
    	legend("topright", as.matrix(title2plot), text.col=col2use)
    }

}

iterate_files <- function (iname, hm.num, s.col, e.col, running.mean, smooth.extent) {
	# For a given dirname of ngs.plot data calculate running means of hm.num heatmaps
	#
	# Args:
	#	iname: the zip file of ngs.plot data
	#	hm.num: the number of heatmap data files in the directory
	#	s.col: the first column in the data frame from which to take values
	#	e.col: the last column in the data frame from to take values
	#   running.mean: the running mean smoothing window size
	#   smooth.extent: a factor to increase smoothing

	#s.col <- 17
	#e.col <- 24

	ifolder <- sub('.zip$', '', basename(iname))
	if(ifolder == basename(iname)) {
    	stop("Input filename must end with .zip\n")
	}

	i <- 1

	if (is.null(hm.num)) {
		hm.num <- 1
	}

	rmeans <- NULL

	col.names <- matrix(0,ncol=hm.num,nrow=1)

	repeat {

		hm <- read.delim( unz(iname, file.path(ifolder, sprintf( "hm%d.txt", i ))), sep="\t", header=T )

		#fpath <- sprintf( "%s/hm%d.txt", ifolder, i )

		#hm <- read.delim( fpath, sep="\t", header=T)

		#extract just the data columns
		hm <- hm[,5:ncol(hm)]

		hm.scores <- as.matrix(hm[,s.col:e.col])
		hm.means <- rowMedians(hm.scores[,1:ncol(hm.scores)])

		## the smoothing (running mean) window
		if ( running.mean == 0 ) {
			sw <- ceiling( nrow(hm.scores) / 180 )
		} else {
			sw <- running.mean
		}

		## how many iterations of the running mean calculation to perform to extend smoothing
		if ( smooth.extent == 0 ) {
			iter <- 1
		} else {
			smooth <- ( (smooth.extent)^2 /2 ) * 3
			iter <- ceiling( sw * smooth )
		}

		hm.sm <- smooth.mean(hm.means, sw, iter)
		col <- sprintf( "V%d",i )

		#cat( sprintf( "\tPerforming %d iterations using a window width of %d on each iteration. ", iter, sw ))
		cat(".")

		if ( !is.matrix(rmeans) ) {
			rmeans <- matrix(0, ncol=hm.num, nrow=nrow(hm.sm))
		}
		rmeans[,i] <- hm.sm
		col.names[,i] <- sprintf( "hm%d.txt", i )

		if (i == hm.num) {
			#cat( "Done\n" )
			colnames(rmeans) <- col.names

			return(rmeans)
		}

		i <- i+1
	}
}

with(updated.vl, {

updated.vl

	if ( winsize > 0 ) {
		s.col <- 1
		e.col <- winsize
	} else {
		cat( sprintf("Average profiles from genomic region %d-%d", s.col, e.col) )
		cat("\n")
	}

	hm.num <- length(ctg.tbl$title)

	cat( sprintf("Calculating running means on heatmaps from %s ", basename(iname) ) )

	maxcol <- 101

	repeat {

		if ( e.col > maxcol ) {
			break
		}

		regcovMat <- iterate_files(iname, length(ctg.tbl$title), s.col, e.col, running.mean, smooth.extent)

		name <- sprintf("%s.avgProfY_X%d-%d", oname, s.col, e.col)

		### write the avg profile to disk
		pdf( sprintf("%s.pdf", name), width=plot.width, height=plot.height )

		plotmat(regcovMat, ctg.tbl$title, yscale, legend, box, lwd, xlab, bam.pair, abline.y)
		out.dev <- dev.off()

		###Write regcovMat to file?
		if ( matrix.out == 1 ) {
			write.table(regcovMat, file=sprintf("%s.txt", name), sep="\t", row.names=F, col.names=T, quote=F)
		}

		###Scanning across the genomic region?
		if ( winsize == 0 ) #no
			break

		s.col <- e.col+1
		e.col <- ifelse( (s.col+winsize)-1 >= maxcol, maxcol, (s.col+winsize)-1 )

		if ( e.col-s.col == 1 || e.col == maxcol) {
			break
		}

	}

	cat( "Done\n" )

})

