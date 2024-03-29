        seed = -1
        seqfile = ../single_copy.cds.concatenated.phy
        treefile = ../mcmc.tree
        outfile = out_usedata2

        ndata = 1
        usedata = 2    * 0: no data; 1:seq like; 2:normal approximation
        clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
        RootAge = '<7'  * constraint on root age, used if no fossil for root.

        model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha = 0.5   * alpha for gamma rates at sites
        ncatG = 5    * No. categories in discrete gamma

        cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

        BDparas = 1 1 0   * birth, death, sampling
        kappa_gamma = 6 2      * gamma prior for kappa
        alpha_gamma = 1 1      * gamma prior for alpha

        rgene_gamma = 1 58.8     * gamma prior for rate for genes
        sigma2_gamma = 1 4.5    * gamma prior for sigma^2     (for clock=2 or 3)

        finetune = 1: 0.06  0.5  0.006  0.12 0.4  * times, rates, mixing, paras, RateParas

        print = 1
        burnin = 50000
        sampfreq = 100
        nsample = 2000000
