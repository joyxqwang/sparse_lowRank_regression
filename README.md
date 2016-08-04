# sparse_lowRank_regression
This function optimizes the following problem:

    \min_{W,Qi} ||X'*W-Y||_F^2 + gamma1*(\sum_i^numG||W*Qi||_Sp^p)^k + gamma2*||W||_{2,q}

Format of input:
    
    n: number of samples
    d: number of SNPs
    c: number of QTs
    X: dim*n SNP data
    Y: n*(cT) phenotype matrix
    numG: number of groups of tasks

Simply run the code in matlab as below:

    [outW, outQ, outObj, outNumIter] = L2QTrL21Squaretrain(X, Y, numG);
