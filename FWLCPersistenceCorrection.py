def FWLCPersistenceCorrection(Length, P_inf=51.51, alpha=2.78):
    """
    Persistence length correction 190121 Steve Yeo
    Based off of Seol et al. Biop. Jour. 2007
    Short DNA Length/50nm ~ 6~20
    Defaults: P_inf = 51.51 (+/- 0.31 nm), alpha = 2.78 (+/- 0.12)
    """
    Pfwlc = P_inf / (1 + (alpha * P_inf) / Length)
    return Pfwlc
