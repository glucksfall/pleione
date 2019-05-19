{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "mawi <- function(alpha,m,n,eps1_,eps2_,x,y)                {\n",
    "eqctr <- 0.5 + (eps2_-eps1_)/2 \n",
    "eqleng <- eps1_ + eps2_\n",
    "\n",
    "wxy <- 0\n",
    "pihxxy <- 0\n",
    "pihxyy <- 0\n",
    "\n",
    "for (i in 1:m)\n",
    "    for (j in 1:n)\n",
    "        wxy <- wxy + trunc(0.5*(sign(x[i] - y[j]) + 1))\n",
    "\n",
    "for (i in 1:m)\n",
    "    for (j1 in 1:(n-1))\n",
    "        for (j2 in (j1+1):n)\n",
    "            pihxyy <- pihxyy + trunc(0.5*(sign(x[i] - max(y[j1],y[j2])) + 1))\n",
    "\n",
    "for (i1 in 1:(m-1))\n",
    "    for (i2 in (i1+1):m)\n",
    "        for (j in 1:n)\n",
    "            pihxxy <- pihxxy + trunc(0.5*(sign(min(x[i1],x[i2]) - y[j]) + 1))\n",
    "\n",
    "wxy <- wxy / (m*n)\n",
    "pihxxy <- pihxxy*2 / (m*(m-1)*n)\n",
    "pihxyy <- pihxyy*2 / (n*(n-1)*m)\n",
    "sigmah <- sqrt((wxy-(m+n-1)*wxy**2+(m-1)*pihxxy+(n-1)*pihxyy)/(m*n))\n",
    "\n",
    "crit <- sqrt(qchisq(alpha,1,(eqleng/2/sigmah)**2))\n",
    "\n",
    "if (abs((wxy-eqctr)/sigmah) >= crit) rej <- 0\n",
    "if (abs((wxy-eqctr)/sigmah) < crit)  rej <- 1\n",
    "\n",
    "if (is.na(sigmah) || is.na(crit)) rej <- 0\n",
    "\n",
    "cat(\" alpha =\",alpha,\"  m =\",m,\"  n =\",n,\"  eps1_ =\",eps1_,\"  eps2_ =\",eps2_,\n",
    " \"\\n\",\"W+ =\",wxy,\"  SIGMAH =\",sigmah,\"  CRIT =\",crit,\"  REJ =\",rej)\n",
    "}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      " alpha = 0.05   m = 12   n = 11   eps1_ = 0.1382   eps2_ = 0.2602 \n",
      " W+ = 0.4242424   SIGMAH = 0.1132147   CRIT = 0.2865093   REJ = 0"
     ]
    }
   ],
   "source": [
    "x <- c(10.3,11.3,2.0,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)\n",
    "y <- c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3.0,6.0,3.5,18.7)\n",
    "mawi(0.05,12,11,0.1382,0.2602,x,y)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "R",
   "language": "R",
   "name": "cran"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "3.5.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
