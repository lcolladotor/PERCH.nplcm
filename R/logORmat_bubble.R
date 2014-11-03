#' Bubble plot of log odds ratio (OR) matrix
#'
#' Plot pairwise log odds ratios for bronze-standard measurements, in controls
#' and case, respectively. The radius of each circle is proportional to
#' \code{1/var(hat{logOR})}. The number in the center of the circle is the \code{hat{logOR}}.
#' The other number above the circle is the \code{s.e.{hat{logOR}}}. We suggest
#' you create a pdf with large width and height if the number of pathogens is large.
#'
#' @param MBS.case Case bronze-standard measurements.
#' @param MBS.ctrl Control bronze-standard measurements.
#' @param pathogen_name The string of pathogen names corresponding to each row/column
#' in the bubble plot.
#'
#' @return A bubble plot. Pairs in cases on the topright; controls at the bottom left.
#' The pathogen names are listed on the left.
#' @export
#'
logORmat_bubble = function(MBS.case,MBS.ctrl,pathogen_name){
  # dat = datobs
  J <- ncol(MBS.case)
  Y   <- c(rep(1,nrow(MBS.case)),rep(0,nrow(MBS.ctrl)))
  Nd  <-  sum(Y)
  MBS <-  as.matrix(rbind(MBS.case,MBS.ctrl))
  Nu <-  length(Y)-Nd

  logORmat = matrix(NA,nrow=J,ncol=J)
  logORmat.se = matrix(NA,nrow=J,ncol=J)

  for (j2 in 1:(J-1)){ #case (j2,j1); ctrl (j1,j2).
    for (j1 in (j2+1):J){

      ## cases:
      x = MBS.case[,j2]
      y = MBS.case[,j1]

      fit = glm(y~x,family = binomial(link="logit"))

      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      ##controls:
      x = MBS.ctrl[,j2]
      y = MBS.ctrl[,j1]

      fit = glm(y~x,family = binomial(link="logit"))

      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j1,j2] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j1,j2] = round(summary(fit)$coef["x",2],3)
      }
    }
  }

  #cell.num = logORmat/logORmat.se
  tmp = logORmat
  tmp[abs(logORmat.se)>10]=NA
  tmp[abs(tmp)<0] = NA
  cell.num.std.logOR = tmp

  tmp2 = logORmat.se
  tmp2[abs(logORmat.se)>10]=NA
  tmp2[abs(tmp2)<0] = NA
  cell.num.prec = 1/tmp2^2

  colors <-rev(brewer.pal(10,"PuOr"))
  pal <- colorRampPalette(colors)

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


   circle.cor = function(cor, cor.se, axes = FALSE, xlab = "",
                        ylab = "", asp = 1,
                        title="",
                        #                         title = "log odds ratio (no adjustment);
                        #                         circle area proportional to 1/var(log OR);
                        #                         color reflects magnitude of log OR in each cell;
                        #                         numbers are log ORs",
                        ...) {
    n = nrow(cor)
    par(mar = c(0, 0, 5, 0), bg = "white")
    plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
         ylab = "", asp = 1, type = "n")
    ##add grid
    segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
             0.5 + 0:n, col = "gray")
    segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                        n), col = "gray")
    ##define circles' background color.
    ##black for positive correlation coefficient and white for negative
    bg = matrix(0,nrow=n,ncol=n)

    for (i in 1:n){
      for (j in 1:n){
        if (!is.na(cor[i,j])){
          bg[i,j] = val2col(z=cor[i,j],c(-10,10),col=pal(64))
        }
      }
    }
    #bg = pal(64)
    #bg[cor > 0] = "black"
    #bg[cor <= 0] = "white"
    ##plot n*n circles using vector language, suggested by Yihui Xie
    symbols(rep(1:n, each = n), rep(n:1, n), add = TRUE, inches = F,
            circles = as.vector(sqrt(abs(cor.se)/max(cor.se[!is.na(cor.se)]))/2),
            bg = as.vector(bg))
    text(rep(0, n), 1:n, n:1, col = "black")
    text(1:n, rep(n + 1), 1:n, col = "black")

    cor.txt<- round(t(cor)[,n:1],1)
    cor.se.txt <-round(t(1/sqrt(cor.se))[,n:1],1)
    for (i in 1:n){
      for (j in 1:n){
        text(i,j,cor.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue"))
        text(i,j+0.25,cor.se.txt[i,j],col="black")
      }
    }
    segments(0.5,.5+n,.5+n,0.5,col="black",lty=3,lwd=3)
    mtext(title,3,cex=1,line=1)
  }

  circle.cor(cell.num.std.logOR,cell.num.prec)

  for (s in rev(1:length(pathogen_name))){
    mtext(paste0(s,":",pathogen_name[s]),side=2,las=2,
          at=par("usr")[1]+0.03*(length(pathogen_name)-s)*diff(par("usr")[1:2])+5,
          cex=1,line=-6)
  }
  mtext("cases",4,cex=2,line=-6)
  mtext("controls",1,cex=2,line=-1)
  mtext(paste(eda_options$X_names,eda_options$X_values,sep="="),3,cex=2,line=-0.5)

}
