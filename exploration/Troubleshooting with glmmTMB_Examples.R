# An example of overly complex model
rm(list = ls())
zinbm0 <- glmmTMB(count~spp + (1|site), zi=~spp, Salamanders, family=nbinom2)
fixef(zinbm0)

ff <- fixef(zinbm0)$zi
signif(plogis(c(sppGP=unname(ff[1]),ff[-1]+ff[1])), 2)

Salamanders <- transform(Salamanders, GP=as.numeric(spp=="GP"))
zinbm0_A <- update(zinbm0, ziformula=~GP)

zinbm0_B <- update(zinbm0, ziformula=~(1|spp))
fixef(zinbm0_B)[["zi"]]
VarCorr(zinbm0_B)

zinbm1 <- glmmTMB(count~spp + (1|site), zi=~mined, Salamanders, family=nbinom2)
fixef(zinbm1)[["zi"]]

## at present we need to specify the parameter by number; for
##  extreme cases need to specify the parameter range
## (not sure why the upper bound needs to be so high ... ?)
cc <- confint(zinbm1,method="uniroot",parm=9, parm.range=c(-20,20))
print(cc)
