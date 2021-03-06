### R code from vignette source 'flup.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: flup.rnw:156-162
###################################################
options( width=90,
#        prompt=" ", continue=" ", # Makes it easier for students to
                                   # copy from the final pdf document
# Makes margins smaller, turns labels and removes the unnecessary box 
         SweaveHooks=list( fig=function()
         par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1,bty="n") ) )


###################################################
### code chunk number 2: flup.rnw:239-241
###################################################
library(Epi)
print( sessionInfo(), l=F )


###################################################
### code chunk number 3: flup.rnw:247-254
###################################################
data( nickel )
nicL <- Lexis( entry = list( per=agein+dob,
                             age=agein,
                             tfh=agein-age1st ),
                exit = list( age=ageout ),
         exit.status = ( icd %in% c(162,163) )*1,
                data = nickel )


###################################################
### code chunk number 4: flup.rnw:264-267
###################################################
str( nickel )
str( nicL )
head( nicL )


###################################################
### code chunk number 5: flup.rnw:276-277
###################################################
summary( nicL )


###################################################
### code chunk number 6: nicL1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot( nicL )


###################################################
### code chunk number 7: nicL2
###################################################
getOption("SweaveHooks")[["fig"]]()
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
plot( nicL, 1:2, lwd=1, col=c("blue","red")[(nicL$exp>0)+1],
      grid=TRUE, lty.grid=1, col.grid=gray(0.7),
      xlim=1900+c(0,90), xaxs="i",
      ylim=  10+c(0,90), yaxs="i", las=1 )
points( nicL, 1:2, pch=c(NA,3)[nicL$lex.Xst+1],
        col="lightgray", lwd=3, cex=1.5 )
points( nicL, 1:2, pch=c(NA,3)[nicL$lex.Xst+1],
        col=c("blue","red")[(nicL$exp>0)+1], lwd=1, cex=1.5 )


###################################################
### code chunk number 8: flup.rnw:331-334
###################################################
nicS1 <- splitLexis( nicL, "age", breaks=seq(0,100,10) )
summary( nicL )
summary( nicS1 )


###################################################
### code chunk number 9: flup.rnw:342-343
###################################################
round( subset( nicS1, id %in% 8:10 ), 2 )


###################################################
### code chunk number 10: flup.rnw:349-351
###################################################
nicS2 <- splitLexis( nicS1, "tfh", breaks=c(0,1,5,10,20,30,100) )
round( subset( nicS2, id %in% 8:10 ), 2 )


###################################################
### code chunk number 11: flup.rnw:356-361
###################################################
library( popEpi )
nicM <- splitMulti( nicL, age = seq(0,100,10), 
                          tfh = c(0,1,5,10,20,30,100) )
summary( nicS2 )
summary( nicM )


###################################################
### code chunk number 12: flup.rnw:365-368
###################################################
identical( nicS2, nicM )
class( nicS2 )
class( nicM )


###################################################
### code chunk number 13: flup.rnw:388-396
###################################################
timeBand( nicM, "age", "middle" )[1:20]
# For nice printing and column labelling use the data.frame() function:
data.frame( nicS2[,c("lex.id","per","age","tfh","lex.dur")],
            mid.age=timeBand( nicS2, "age", "middle" ),
              mid.t=timeBand( nicS2, "tfh", "middle" ),
             left.t=timeBand( nicS2, "tfh", "left"   ),
            right.t=timeBand( nicS2, "tfh", "right"  ),
             fact.t=timeBand( nicS2, "tfh", "factor" ) )[1:20,]


###################################################
### code chunk number 14: flup.rnw:416-417
###################################################
summary( (nicS2$age-nicS2$tfh) - nicS2$age1st ) 


###################################################
### code chunk number 15: flup.rnw:422-424
###################################################
summary( timeBand( nicS2, "age", "middle" ) -
         timeBand( nicS2, "tfh", "middle" ) - nicS2$age1st )


###################################################
### code chunk number 16: flup.rnw:446-454
###################################################
subset( nicL, id %in% 8:10 )
agehi <- nicL$age1st + 50 / nicL$exposure
nicC <- cutLexis( data = nicL, 
                   cut = agehi, 
             timescale = "age",
             new.state = 2, 
      precursor.states = 0 )
subset( nicC, id %in% 8:10 )


###################################################
### code chunk number 17: flup.rnw:461-469
###################################################
subset( nicS2, id %in% 8:10 )
agehi <- nicS2$age1st + 50 / nicS2$exposure
nicS2C <- cutLexis( data = nicS2, 
                     cut = agehi, 
               timescale = "age",
               new.state = 2, 
        precursor.states = 0 )
subset( nicS2C, id %in% 8:10 )


###################################################
### code chunk number 18: flup.rnw:528-530
###################################################
( a.kn <- with( subset( nicM, lex.Xst==1 ), quantile( age+lex.dur, (1:5-0.5)/5 ) ) )
( t.kn <- with( subset( nicM, lex.Xst==1 ), quantile( tfh+lex.dur, (1:5-0.5)/5 ) ) )


###################################################
### code chunk number 19: flup.rnw:543-548
###################################################
ma <- glm( (lex.Xst==1) ~ Ns(age,knots=a.kn),
           family = poisson,
           offset = log(lex.dur),
             data = nicM )
summary( ma )


###################################################
### code chunk number 20: pr-a
###################################################
getOption("SweaveHooks")[["fig"]]()
nd <- data.frame( age=40:85, lex.dur=1000 )
pr.a <- ci.pred( ma, newdata = nd )
matplot( nd$age, pr.a,
         type="l", lty=1, col=1, lwd=c(3,1,1),
         log="y", xlab="Age (years)",
         ylab="Lunng cancer mortality per 1000 PY")


###################################################
### code chunk number 21: flup.rnw:583-585
###################################################
mat <- update( ma, . ~ . + Ns(tfh,knots=t.kn) )
summary( mat )


###################################################
### code chunk number 22: flup.rnw:595-596
###################################################
summary( nickel$age1st )


###################################################
### code chunk number 23: flup.rnw:600-606
###################################################
nd <- data.frame( expand.grid( age=c(20:90,NA), age1st=seq(15,45,10) ) ) 
nd <- transform( nd, tfh = ifelse( age > age1st, age-age1st, NA ),
                 lex.dur = 1000 )
# makes no sense to have age < age1st
nd <- transform( nd, age = ifelse( age > age1st, age, NA ) )
head( nd )


###################################################
### code chunk number 24: pr-at
###################################################
getOption("SweaveHooks")[["fig"]]()
pr.at <- ci.pred( mat, newdata = nd )
matplot( nd$age, pr.at,
         type="l", lty=1, col=1, lwd=c(3,1,1),
         log="y", xlab="Age (years)",
         ylab="Lunng cancer mortality per 1000 PY")


###################################################
### code chunk number 25: flup.rnw:630-631
###################################################
anova( ma, mat, test="Chisq" )


###################################################
### code chunk number 26: flup.rnw:642-646
###################################################
( f.kn <- with( subset( nicM, lex.Xst==1 ), quantile( age1st, (1:5-0.5)/5 ) ) )
maf <- update( ma, . ~ . + Ns(age1st,knots=f.kn) )
summary( maf )
anova( maf, ma, mat, test="Chisq" )


###################################################
### code chunk number 27: pr-at-af
###################################################
getOption("SweaveHooks")[["fig"]]()
pr.af <- ci.pred( maf, newdata = nd )
matplot( nd$age, pr.at,
         type="l", lty=1, col=1, lwd=c(3,1,1),
         log="y", xlab="Age (years)",
         ylab="Lunng cancer mortality per 1000 PY")
matlines( nd$age, pr.af,
         type="l", lty=1, col=2, lwd=c(3,0,0) )


###################################################
### code chunk number 28: flup.rnw:674-685
###################################################
maft <- update( mat, . ~ . + Ns(age1st,knots=f.kn) )
summary( maft )
mft <- update( maft, . ~ . - Ns(age,knots=a.kn) )
mf  <- update( maf , . ~ . - Ns(age,knots=a.kn) )
mt  <- update( mat , . ~ . - Ns(age,knots=a.kn) )
allp <- anova( maft, mat, ma, maf, mf, mft, mt, mat,
               maf, maft, mft,
               test="Chisq" )
mall <- as.matrix( allp )
cbind( mod = c("maft","mat","ma","maf","mf","mft","mt","mat","maf","maft","mft"),
       round( allp[,1:5], 3 ) )


###################################################
### code chunk number 29: flup.rnw:700-709
###################################################
data( nickel )
nicL <- Lexis( entry = list( per = agein+dob,
                             age = agein,
                             tfh = agein-age1st ),
                exit = list( age = ageout ),
         exit.status = ( icd > 0 ) + ( icd %in% c(162,163) ),
                data = nickel )
summary( nicL )
subset( nicL, id %in% 8:10 )


###################################################
### code chunk number 30: flup.rnw:714-723
###################################################
nicL <- Lexis( entry = list( per = agein+dob,
                             age = agein,
                             tfh = agein-age1st ),
                exit = list( age = ageout ),
         exit.status = ( icd > 0 ) + ( icd %in% c(162,163) ),
                data = nickel,
              states = c("Alive","D.oth","D.lung") )
summary( nicL )
str( nicL )


###################################################
### code chunk number 31: flup.rnw:735-743
###################################################
nicL$agehi <- nicL$age1st + 50 / nicL$exposure
nicC <- cutLexis( data = nicL,
                   cut = nicL$agehi,
             timescale = "age",
             new.state = "HiExp",
      precursor.states = "Alive" )
subset( nicC, id %in% 8:10 )
summary( nicC, scale=1000 )


###################################################
### code chunk number 32: flup.rnw:762-771
###################################################
nicC <- cutLexis( data = nicL,
                   cut = nicL$agehi,
             timescale = "age",
             new.state = "HiExp",
             new.scale = "tfe",
          split.states = TRUE,
      precursor.states = "Alive" )
subset( nicC, id %in% 8:10 )
summary( nicC, scale=1000, timeScales=TRUE )


###################################################
### code chunk number 33: nic-box
###################################################
getOption("SweaveHooks")[["fig"]]()
boxes( nicC, boxpos = list(x=c(10,10,80,80,80,80),
                           y=c(75,25,87,63,13,37)),
            scale.Y = 1000,
            show.BE = TRUE )


