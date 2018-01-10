# Define function to calculate tissue specificity
tissueScore <- function(x){
        x_norm <- (log2(x+1))/sum(log2(x+1))
        ntissues <- length(x_norm)
        H <- function(p){
                p <- p[p>0]
                -sum(p*log(p))
        }

        JSD <- function(p1, p2){
                H((p1+p2)/2) - ((H(p1) + H(p2))/2)
        }

        scores <- as.numeric(rep(NA,ntissues))
        names(scores) <- names(x_norm)
        for(i in 1:ntissues){
                et <- rep(0,ntissues)
                et[i] <- 1
                scores[[i]] <- 1-sqrt(JSD(x_norm, et))
        }
        maxScore <- scores[which.max(scores)]
        return(c(Tissue=names(maxScore)[1], Score=maxScore, maxFPKM=x[which.max(x)]))

}


#this function converts a vector of values("z") to a vector of color
#levels. One must define the number of colors. The limits of the color
#scale("zlim") or the break points for the color changes("breaks") can 
#also be defined. when breaks and zlim are defined, breaks overrides zlim.
# From: http://menugget.blogspot.co.uk/2011/09/converting-values-to-color-levels.html
val2col<-function(z, zlim, col = heat.colors(12), breaks){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 CUT <- cut(z, breaks=breaks)
 colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
 return(colorlevels)
}

# Useful functions for scatterplot pairs
upPan <- function(...){
	points(..., col = "darkblue", pch=".")
	abline(a=0, b=1, col = "red")
}
lowPan <- function(x, y, ...) {
	text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
	signif(cor(x, y), 2), cex = 1.5)
}

fancy_scientific <- function(l) { 
     # turn in to character string in scientific notation 
     l <- format(l,  scientific = TRUE) 
     # quote the part before the exponent to keep all the digits 
     l <- gsub("^(.*)e", "'\\1'e", l) 
     # turn the 'e+' into plotmath format 
     l <- gsub("e", "%*%10^", l) 
     l <- gsub("10\\^\\+", "10^", l)
     # return this as an expression 
     parse(text=l) 
} 

