import uncertainties


def cfloat(string):
    """
    Create "certain float" (`ufloat.n`) from string by converting it first
    to `uncertainties.ufloat` and then taking the nominal value.

    :param string: string to be converted
    :type string: str
    :return: float with no uncertainty
    :rtype: float
    """
    return uncertainties.ufloat_fromstr(representation=string).n
