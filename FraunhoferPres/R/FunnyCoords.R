
setwd("/home/tim/git/FraunhoferPres/FraunhoferPres")

# time (APCTDL) can be an absolute coordinate, but anything that maps
# to a time dimension can itself become a coordinate on which to 
# project other values.

# example of a few projections, where age or time have been swapped out with 
# standardized lifetable functions.


library(HMDHFDplus)

# define username and password in console
ltf <- readHMDweb("USA","fltper_1x1",username=us,password=pw)
fert <- readHFDweb("USA","asfrRR",username=us,password=pw)

# a prospective age calculator, given e(x)
ex2age <- function(ex,Age = (1:length(ex)-1), ex.at = seq(50,5,by=-5)){
	# splines make it easier to interpolate to exact values.
	out <- splinefun(Age~ex)(ex.at)
	names(out) <- ex.at
	out
}
# a quantile age calculator, given l(x)
lxquantile2age <- function(lx, Age = (1:length(lx)-1), p = seq(.9,.1,by=-.1), closeout = 2){
	# includes closing out the lifetable at open age + 2, whatever. Only OK
	# if your lifetable goes to 110+ (could also assume exponential or something)
	# same use of spline, but we force monotonic.
	out <- splinefun(c(Age,max(Age)+closeout)~c(lx,0),method="monoH.FC")(p)
	names(out) <- p
	out
}

Age <- 0:110
yr1 <- 1950 ; yr2 <- 2010

# a look at selected period prospective ages
ex1 <- ltf$ex[ltf$Year == yr1]
ex2 <- ltf$ex[ltf$Year == yr2]

prosp <- cbind(ex = seq(50,5,by=-5), 
		age1 = round(ex2age(ex1),1),
		age2 = round(ex2age(ex2),1))
head(prosp)
#   ex age1 age2
#50 50 24.1 32.3
#45 45 29.4 37.5
#40 40 34.8 42.8
#35 35 40.3 48.2
#30 30 45.9 53.8
#25 25 51.9 59.6
# or a look at period quantile ages
lx1 <- ltf$lx[ltf$Year == yr1] / 1e5
lx2 <- ltf$lx[ltf$Year == yr2] / 1e5

quant <- cbind(prob = seq(.9,.1,by=-.1), 
		age1 = round(lxquantile2age(lx1),1),
		age2 = round(lxquantile2age(lx2),1))
head(quant)
#    prob age1 age2
#0.9  0.9 48.1 62.5
#0.8  0.8 60.8 72.1
#0.7  0.7 67.5 77.6
#0.6  0.6 72.2 81.6
#0.5  0.5 75.8 84.7
#0.4  0.4 79.0 87.5

# you didn't think I'd not try to make a surface out
# of this stuff, did you?
# [period assumption sinner here]
library(reshape2)
# an age-period matrix of e(x)
EX <- acast(ltf,Age~Year,value.var = "ex")
LX <- acast(ltf,Age~Year,value.var = "lx")
# a big matrix of prospective ages, sort of:
# we're using thresholds rather than comparing
# with a standard, but in the end such comparisons
# will be possible.
Sand    <- apply(EX,2,ex2age,ex.at=5:50)
QT    <- apply(LX/1e5,2,lxquantile2age,p=seq(.9,.1,by=-.01))
# define a color ramp, because colors
colRamp <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"),space="Lab")

# a custom cohort line function, 
# not necessarily robust or user-friendly.
# made for the sake of this single plot...
LexLines <- function(EX,N=5,...){
	Sand <- apply(EX,2,ex2age,ex.at=5:50)
	# make it hmmm.
	ages <- as.integer(rownames(EX))
	aN <- ages[ages %% N == 0]
	aN <- aN[aN < max(Sand)]
	aN <- aN[aN > min(Sand)]
	years <- as.integer(colnames(EX)) 
	yN <-  years[years %% N == 0]
	
	EXN <- EX[as.character(aN),as.character(yN)]
	for (m in 1:(nrow(EXN)-1)){
		for (n in 1:(ncol(EXN)-1)){
			segments(yN[n],EXN[m,n],yN[n+1],EXN[m+1,n+1], xpd=FALSE, ...)
		}
	}
	# left side:
	offl <- min(yN) - min(years)
	if (offl > 0){
		EXleft <- EX[as.character(aN+offl),1]
		for (m in 1:(nrow(EXN)-1)){
			segments(min(years),EXleft[m],yN[1],EXN[m+1,1], xpd=FALSE, ...)
		}
	}
	# right side
	offr <- max(years) -  max(yN)
	if(offr > 0){
		EXright <- EX[as.character(aN-offl),ncol(EX)]
		for (m in 1:(nrow(EXN)-1)){
			segments(max(yN),EXN[m,ncol(EXN)], max(years), EXright[m+1], xpd=FALSE, ...)
		}
	}
	cohL <- min(yN) - aN
	
	for (i in 1:(length(cohL)-1)){
		rise <- EXN[i+1,1] - EXN[i,1]
		run <- N
		slope <- rise / run
		angle <- atan(slope) * 180 / pi
		text(min(yN), EXN[i,1]+1, cohL[i], col = "yellow", cex = .8, srt = angle)	
	}
	
	cohR <- max(yN) - aN
	
	for (i in 1:(length(cohR)-1)){
		rise <- EXN[i+1,ncol(EXN)] - EXN[i,ncol(EXN)]
		run <- N
		slope <- rise / run
		angle <- atan(slope) * 180 / pi
		text(max(yN), EXN[i,ncol(EXN)]+1, cohR[i], col = "yellow", cex = .8, srt = angle)	
	}
	
}

getwd()
# the complicated plot, will describe in blog.
#png("Figures/prospectiveAgeLexis.png",width=1000,height=600)
pdf("Figures/prospectiveAgeLexis.pdf",width=10,height=6.1)
filled.contour(as.integer(colnames(Sand)),as.integer(rownames(Sand)),t(Sand),
		asp=1,color=colRamp, xlab = "Year",
		ylab = "remaining life expectancy",key.title = title("Age"),
		# but quick tip: since filled.contour() spits back a device with crazy coords,
		# if you want to plot more in the figure, you need to do it by passing arguments
		# to plot.axes (non-intuitively). It's just like panel.first = list(bla bla)
		plot.axes = { LexLines(EX,N=5,col="#FFFF0080");
			contour(as.integer(colnames(Sand)),as.integer(rownames(Sand)),t(Sand), 
					drawlabels = TRUE, axes = FALSE, 
					frame.plot = FALSE, add = TRUE,
					labcex = 1);
			axis(1); axis(2);
			# grids are helpful in this case to match prospective ages...
			abline(h=seq(5,50,by=5),col="#FFFFFF30");
			abline(v=seq(1935,2015,by=5),col="#FFFFFF30")})
dev.off()
# one could do this for quantiles as well, maybe some other time.


ltf <- readHMDweb("SWE","fltper_1x1",username=us,password=pw)
fert <- readHFDweb("SWE","asfrRR",username=us,password=pw)
FX       <- acast(fert,Age~Year,value.var="ASFR")
LX <- acast(ltf[ltf$Year %in% fert$Year,],Age~Year,value.var = "lx") / 1e5
#Quant    <- apply(LX,2,lxquantile2age,p=seq(0,1,by=.01))
#
#filled.contour(as.integer(colnames(Quant)),as.numeric(rownames(Quant))*100,t(Quant),
#		asp=1,color=colRamp, xlab = "Year",
#		ylab = "percent surviving until at least age",key.title = title("Age"),
#		plot.axes = {contour(as.integer(colnames(Quant)),as.numeric(rownames(Quant))*100,t(Quant), 
#					drawlabels = TRUE, axes = FALSE, 
#					frame.plot = FALSE, add = TRUE,
#					labcex = 1); 
#			axis(1); axis(2)})

# but what about a function to assign fert to the quantiles?
# a quantile age calculator, given l(x)
age2fert <- function(age=12.5:55.5, fertvec, ageout=12.5:55.5){
	# includes closing out the lifetable at open age + 2, whatever. Only OK
	# if your lifetable goes to 110+ (could also assume exponential or something)
	# same use of spline, but we force monotonic.
	out <- splinefun(c(0,fertvec,0)~c(age[1]-1,age,age[length(age)]+1),method="monoH.FC")(ageout)
	names(out) <- ageout
	out[out<0] <- 0
	out
}
p <- seq(.5,1,by=.005)
Quant    <- apply(LX,2,lxquantile2age,p=p)
FX       <- acast(fert,Age~Year,value.var="ASFR")
QuantF   <- Quant * 0
for (yr in 1:ncol(Quant)){
	QuantF[,yr] <- age2fert(age=12.5:55.5,fertvec=FX[,yr],ageout=lxquantile2age(LX[,yr],p=p))
}


pdf("Figures/FertQuant.pdf",width=10,height=6)
par(mai=c(.8,.8,.5,.5))
filled.contour(as.integer(colnames(QuantF)),as.numeric(rownames(QuantF))*100,t(QuantF),
		color=colRamp, xlab = "Year",
		ylab = "percent surviving",key.title = title("Fertility
rate"),plot.axes = { contour(as.integer(colnames(QuantF)),as.numeric(rownames(QuantF))*100,t(QuantF), 
			drawlabels = TRUE, axes = FALSE, 
			frame.plot = FALSE, add = TRUE,
			labcex = 1, col = "#00000040");
	axis(1); axis(2)})
dev.off()



pdf("Figures/FertAPC.pdf",width=10,height=6)
par(mai=c(.8,.8,.5,.5), xaxs='i',yaxs='i')
filled.contour(as.integer(colnames(FX)),as.numeric(rownames(FX)),t(FX),
		asp=1,color=colRamp, xlab = "Year",
		ylab = "Age",key.title = title("Fertility
						rate"),ylim=c(0,70),
		plot.axes = { contour(as.integer(colnames(FX)),as.integer(rownames(FX)),t(FX), 
							drawlabels = TRUE, axes = FALSE, 
							frame.plot = FALSE, add = TRUE,
							labcex = 1, col = "#00000040");
					axis(1); axis(2)})
dev.off()


#############################3
# for blog:

png("Figures/FertQuant.png",width=1000,height=600)
par(mai=c(.8,.8,.5,.5))
filled.contour(as.integer(colnames(QuantF)),as.numeric(rownames(QuantF))*100,t(QuantF),
		color=colRamp, xlab = "Year",
		ylab = "percent surviving",key.title = title("Fertility
						rate"),plot.axes = { contour(as.integer(colnames(QuantF)),as.numeric(rownames(QuantF))*100,t(QuantF), 
					drawlabels = TRUE, axes = FALSE, 
					frame.plot = FALSE, add = TRUE,
					labcex = 1, col = "#00000040");
			axis(1); axis(2)})
dev.off()

png("Figures/FertAPC.png",width=1000,height=600)
par(mai=c(.8,.8,.5,.5), xaxs='i',yaxs='i')
filled.contour(as.integer(colnames(FX)),as.numeric(rownames(FX)),t(FX),
		asp=1,color=colRamp, xlab = "Year",
		ylab = "Age",key.title = title("Fertility
						rate"),ylim=c(0,70),
		plot.axes = { contour(as.integer(colnames(FX)),as.integer(rownames(FX)),t(FX), 
					drawlabels = TRUE, axes = FALSE, 
					frame.plot = FALSE, add = TRUE,
					labcex = 1, col = "#00000040");
			axis(1); axis(2)})
dev.off()
