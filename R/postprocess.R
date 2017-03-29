#' @import data.table
#' @import magrittr
#' @importFrom stringr str_count
#' @importFrom cowplot plot_grid
#' @importFrom Rcpp sourceCpp
NULL

################################################################################

## reading

PI <- c(0.91243863, 0.08756137)
rblist <- function(L) rbindlist(L,use.names=TRUE) #do.call("rbind",L)
##' Read abcGWAS output for a (list of) directory(ies) in a defined format
##'
##' @title reader
##' @param d directory or list of directories
##' @param ... optional parameters passed to inner functions.  Most
##'     useful for test runs is setting n, for example setting n=10
##'     will use only the first 10 files in each directory
##' @return a data.table summarising the passed results and coverage
##'     so far.  This has class c("ABC","data.table") and some
##'     convenience functions make use of this additional S3 class.
##' @export
##' @author Chris Wallace
reader <- function(d,...) {
    message("reading ABC passes")
    x <- lapply(d,readd,...)%>%rblist()
    message("reading reference data")
    ntlist <- lapply(d,readR,...)
    nt1 <- ntlist[[1]]$nt1
    nt2 <- ntlist[[1]]$nt2
    thr <- ntlist[[1]]$thr
    snps <- ntlist[[1]]$snps
    if(length(ntlist)>1) {
        for(i in 2:length(ntlist)) {
            if(!identical(snps,ntlist[[i]]$snps))
                stop("snps mismatch between d[[1]] and d[[",i,"]]")
            nt1 <- nt1 + ntlist[[i]]$nt1
            nt2 <- nt2 + ntlist[[i]]$nt2
        }
    }
    ref <- if(all(is.na(x$Var2))) {
               onedt(nt1,snps)
           } else {
               rbind(onedt(nt1,snps),twodt(nt2,snps))
           }
    vname <- grep("Var",names(ref),value=TRUE)
    x <- merge(x,ref,by=c(vname,"ncv"),all=TRUE)
    x[,ssthr:=cut(sumsq, c(0,thr), include.lowest=TRUE)]
    attributes(x) <- c(attributes(x),list(thresholds=thr))
    class(x) <- c("ABC",class(x))
    x
}
readd <- function(d,n=NULL) {
    files <- list.files(d,full=TRUE,pattern="RData$") %>% sub(".RData",".out",.)    
    if(!is.null(n) && n<length(files))
        files <- files[1:n]
    message(d, ":\t", length(files), " files found")
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    ret <- lapply(seq_along(files),function(i) {
        #setTxtProgressBar(pb, i)
        x <- fread(files[i])
        x$ncv <- str_count(x$CV," ")+1
        x
    })%>%rbindlist()
    close(pb)
    maxcv <- max(ret$ncv)
    if(maxcv==1) {
        setnames(ret,"CV","Var1")
        ret[,Var2:=NA]
    } else {
        ret[,Var1:=sub(" .*","",CV)]
        ret[ncv==2,Var2:=sub(".* ","",CV)]
        ret[ncv==1,Var2:=NA]
        ret[,CV:=NULL]
    }
    ret[,c("Var1","Var2"):=list(as.character(Var1),as.character(Var2))]
    ret[Var1<Var2, c("Var1","Var2"):=list(Var2,Var1)]
    ret
}
onedt <- function(nt,snps) {
    names(nt) <- snps
    mt <- data.table(Var1=names(nt),Var2=as.character(NA),count=nt,ncv=1,prior=PI[1]/length(snps))
    mt[count>0,weight:=prior/count]
    mt
}
twodt <- function(nt,snps) {
    dimnames(nt) <- list(snps,snps)
    mt2 <- as.data.table(melt(nt,value.name="count"))
    mt2[,c("Var1","Var2"):=list(as.character(Var1),as.character(Var2))]
    mt2[Var1<Var2, c("Var1","Var2"):=list(Var2,Var1)]
    mt2 <- unique(mt2,by=c("Var1","Var2"))
    mt2[,ncv:=2]
    N <- length(snps)
    mt2[,prior:=PI[2]/(N*(N-1)/2)]
    mt2[count>0,weight:=prior/count]
    mt2
}
##' Given a directory of abcGWAS output, optionally shrink the .RData files
##'
##' Shrunk files will retain only a subset of the information, and
##' will be copied to a sub-directory, .shrink, of d.  Doing this once
##' will make subsequent reading of data faster.
##' @title preshrink.
##' @param d directory containing abcGWAS output
##' @param newsnps don't use (here for test purposes only, will be
##'     deleted)
##' @return 0 if success
##' @export
##' @author Chris Wallace
preshrink <- function(d,newsnps=NULL) {
    files <-  list.files(d,pattern="RData")
    nd <- file.path(d,".shrink")
    if(!file.exists(nd))
        dir.create(nd)
    newfiles <- list.files(nd,pattern="RData")
    files <- setdiff(files,newfiles) # only shrink once
    if(!length(files))
        return(0)
    message("shrinking files in ",d)
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for(i in seq_along(files)) {
        setTxtProgressBar(pb, i)
        f <- file.path(d,files[i])
        nf <- file.path(nd,files[i])
        (objs <- load(f))
        if("snps.all" %in% objs)
            snps <- snps.all
        if(!is.null(newsnps))
            snps <- newsnps
        oname <- grep("numtested",objs,value=TRUE)
        save(list=c(oname,"snps","thr"),file=nf)
    }
    close(pb)
    0
}
fcombine <- function(fuse,fbak) {
    buse <- basename(fuse)
    bbak <- basename(fbak)
    usebak <- which(!(bbak %in% buse))
    c(fuse,fbak[usebak])
}
readR <- function(d,n=NULL) {
    files <- list.files(d,full=TRUE,pattern="RData")
    nd <- file.path(d,".shrink")
    if(file.exists(nd)) {
        newfiles <- list.files(nd,full=TRUE,pattern="RData")
        files <- fcombine(newfiles,files)
    }
    if(!is.null(n) && n<length(files))
        files <- files[1:n]
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    (objs <- load(files[1]))
   ## 1CV
    if("numtested" %in% objs) {
        nt1 <- numtested
        nt2 <- matrix(0,length(snps),length(snps))
        for(i in 2:length(files)) {
           setTxtProgressBar(pb, i)
           load(files[i])
            nt1 <- nt1 + numtested
        }
#        return(list(ref=onedt(nt,snps),thr=thr))
    } else {
    ## 2CV
        ##    print(objs)
        nt1 <- numeric(length(snps))
        if("numtested_1CV" %in% objs)
            nt1 <- numtested_1CV
        nt2 <- numtested_2CV
        for(i in 2:length(files)) {
            setTxtProgressBar(pb, i)
            (load(files[i]))
            if("numtested_1CV" %in% objs)
                nt1 <- nt1 + numtested_1CV
            nt2 <- nt2 + numtested_2CV
        }
    }
    close(pb)    
    return(list(nt1=nt1,nt2=nt2,snps=snps,thr=thr))
}

################################################################################

## qc
##' print a brief summary of an abcGWAS run
##'
##' @title summary
##' @param x object of class ABC
##' @return no return value
##' @export
##' @author Chris Wallace
summary.ABC <- function(x) {
    message("ABC object containing results for ",sum(x$count>0)," samples")
    print(head(x))
    coverage(x)
}
##' print coverage summary, where coverage is defined by the number of times a model was visted
##'
##' @title coverage of abcGWAS run
##' @inheritParams summary.ABC
##' @return no return value
##' @export
##' @author Chris Wallace
coverage <- function(x) {
    ux1 <- unique(x[ncv==1,],by="Var1")
    message("1CV model coverage:")
    print(summary(ux1$count))
    ux2 <- unique(x[ncv==2,],by=c("Var1","Var2"))
    if(nrow(COV$ux2)) {
        message("2CV model coverage:")
        print(summary(ux2$count))
    }
    invisible(list(ux1=ux1,ux2=ux2))
}

################################################################################

## plotting
##' diagnostic plot to decide whether posterior estimates have converged
##'
##' @title plot.convergence
##' @inheritParams summary.ABC
##' @param ithr which threshold to use. Optional.  Default is second element in thr.
##' @return plot
##' @export
##' @author Chris Wallace
plot.convergence <- function(x,ithr=2) {
    x <- x[sample(1:nrow(x)),]
    x2 <- x[as.numeric(ssthr)==ithr & count>0,]
    splits <- quantile(1:nrow(x2),seq(0.1,1,by=0.1))
    po <- lapply(splits, function(n) {
        x3 <- x2[1:n,.(po=sum(weight)),by=c("Var1","Var2")]
        x3$po <- x3$po/sum(x3$po)
        x3 <- stack.ABCpost(x3)
        x3 <- x3[,.(po=sum(po)),by="Var"]
        x3$split <- n
        return(x3)
    })%>%rblist()
    po[,maxpo:=max(po),by="Var"]
    po <- po[maxpo>quantile(maxpo,0.9) | maxpo>0.01,]
    po <- po[order(split),]
    ggplot(po,aes(x=split,y=po,col=Var,group=Var)) + geom_path()
}
##' @title plotting
##' @inheritParams summary.ABC
##' @param what aspect of abcGWAS run to plot
##' * coverage call plot.coverage (default)
##' * snp plot mppi summary by threshold
##' * ncv plot ncv posterior by threshold
##' * model plot model posterior by threshold
##' @return the plot
##' @export
##' @author Chris Wallace
plot.ABC <- function(x,what=c("coverage","nsnp","ncv","model")){
    what <- match.arg(what)
    switch(what,
           "coverage"=plot.coverage(x),
           "nsnp"=plot(post(x,by="nsnp")),
           "ncv"=plot(post(x,by="ncv")),
           "model"=plot(post(x,by="model")),
           stop("option not recognised: ",what))
}
##' Plot posterior summary
##'
##' @title plot posterior inference
##' @param y object of class postABC
##' @param po.lower lower threshold - do not display values with posterior probability or (MPPI) < po.lower
##' @return a plot
##' @export
##' @author Chris Wallace
plot.ABCpost <- function(y,po.lower=0) {
    y <- y[order(y$ssthr),]
    yvar <- intersect(names(y),c("mppi","po"))
    colvar <- intersect(names(y),c("CV","Var","ncv"))
    ncols <- length(unique(y[po>po.lower,][[colvar]]))
    y[[colvar]] <- as.factor(y[[colvar]])
    p <- ggplot(y[po>po.lower,], aes_string(x="ssthr",y=yvar,col=colvar,group=colvar)) + geom_path()
    if(ncols>9)
        p <- p + theme(legend.position="none")
    return(p)
}
plot.gamma1 <- function(x) {
    y <- post(x,by="snp")
    use <- unique(y[mppi>0.01,]$Var)
    ggplot(x[Var1 %in% use,],aes(x=gamma,fill=Var1)) + geom_histogram() + facet_wrap(~Var1)
}

##' plot histogram of coverage
##'
##' @title coverage of abcGWAS run
##' @inheritParams summary.ABC
##' @return no return value
##' @export
##' @author Chris Wallace
plot.coverage <- function(x) {
    COV <- coverage(x)
    mx <- max(COV$ux1$count)
    if(nrow(COV$ux2)) {
        mx <- max(c(mx,ux2$count))
    }
    p1 <- ggplot(COV$ux1,aes(x=count)) + geom_histogram(binwidth=500) + geom_vline(xintercept=0,col="red",linetype="dashed") + ggtitle("Number of one CV samples per model") 
    if(!nrow(ux2))
        return(p1)
    p2 <- ggplot(COV$ux2,aes(x=count)) + geom_histogram(binwidth=500) + geom_vline(xintercept=0,col="red",linetype="dashed") + ggtitle("Number of two CV samples per model")
    plot_grid(p1,p2,ncol=1,align="h")
}



################################################################################

## posterior inference
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @inheritParams summary.ABC
##' @param by choose to make posterior inference on
##' * ncv the number of causal variants
##' * model the choice of W (Var1 and Var2)
##' * snp generate the marginal posterior probability of inclusion for each SNP
##' @return object of class ABCpost
##' @export
##' @author Chris Wallace
post <- function(x, by=c("ncv","model","snp")) {
    by <- match.arg(by)
    dtby <- switch(by,
                   ncv="ncv",
                   c("Var1","Var2"))
    y <- x[count>0,.(pass=sum(weight),prior=prior[1]),by=c("ssthr",dtby)]
    y <- y[order(ssthr),]
    y[,pass:=cumsum(pass),by=dtby]
    y[,po:=pass/sum(pass),by="ssthr"]
    y <- y[order(po,decreasing=TRUE),]
    if(by %in% c("ncv","model")) {
        class(y) <- c("ABCpost",class(y))
        if(by=="ncv")
            y[,prior:=PI[y$ncv]]
        return(y)
    }
    if(all(is.na(y$Var2))) {
        setnames(y,c("Var1"),c("Var"))
        y[,prior:=NULL]
        y[,pass:=NULL]
        y[,Var2:=NULL]
        class(y) <- c("mppi","ABCpost",class(y))
        return(y)
    }
    po <- stack.ABCpost(y)
    po <- po[,.(po=sum(po)),by=c("Var","ssthr")]
    class(po) <- c("mppi","ABCpost",class(po))
    po[order(po,decreasing=TRUE),]    
}

stack.ABCpost <- function(y) {
    po1 <- y
    po2 <- y[!is.na(Var2),]
    setnames(po1,"Var1","Var")
    po1[,Var2:=NULL]
    setnames(po2,"Var2","Var")
    po2[,Var1:=NULL]
    po <- rbind(po1,po2)
}
##' Generate 99% credible set of models
##'
##' .. content for \details{} ..
##' @title
##' @inheritParams summary.ABC
##' @param ... other arguments passed to methods, specifically size which sets the desired credible set size
##' @return
##' @export
##' @author Chris Wallace
credset <- function(x,...) UseMethod("credset",x)
    
##' @export
##' @author Chris Wallace
credset.ABC <- function(x,size=0.99) {
    po <- post(x,by="model")
    credset(po)
}
##' @export
##' @author Chris Wallace
credset.postABC <- function(x,size=0.99) {
    x[,cumpp:=cumsum(po)]
    w <- which(x$cumpp>size)[1]
    x[1:w,]
}
