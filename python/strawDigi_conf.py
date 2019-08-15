import shipunit as u

class StrawtubesMisalign:
    # TurnOn the module of misalign or no
    misalign = True
    # rand type : "None", "Gaus", "Unif", may added more, need modify strawtubeDigi
    # "None" = all tubes have same maximum sagging, no distribution of max. sag
    # "Gaus" = the tubes have different max sagging with gaus distribution
    # "Unif" = the different max sagging is uniformly distributed in a range given range, not the same as "None"
    randType = "Unif"

    # maximum sagging at the middle of the tube
    # for using distribution, these are mean value /mpv /something like that
    maxTubeSagging = 0.0*u.cm
    maxWireSagging = 0.7*u.cm

    # uniform distribtion, the delta is half of the range
    tubeUnifDelta = 0.0*u.cm
    wireUnifDelta = 0.0*u.cm

    # Gaussian distribtion, the sigma
    tubeGausSigma = 0.*u.cm
    wireGausSigma = 0.*u.cm

    # debug or not
    debug = False

# Using which way to calculate drift time
class DriftTimeCalculate:
    defaultDriftTime = True   # true for using drift velocity (default)


