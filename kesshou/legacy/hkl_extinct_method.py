def extinct(self, rule='hkl:', point_group=PG['1']):
    """
    Removes from dataframe all reflections which are in a specified
    *domain*, but do not meet the *condition*.
    The *rules* have a format "domain:condition",
    where the domain specifies a part of reciprocal space which
    we plan to extinct, while the condition describes all the reflections,
    which should **not** be extinct despite belonging to domain.

    The rule should be written using a format present in
    International Tables of Crystallography,
    which is easily accessible using the web page of
    `Department of Chemistry, University College London
    <http://img.chem.ucl.ac.uk/sgp/large/sgp.htm>`_
    In the current release the *domains* are hard-coded and thus not all
    theoretically possible domains can be specified.
    However, all domains specified in the International Tables
    should be accessible.

    :param rule: A string with *domain* and *condition* separated by colon.
    :type rule: str
    :param point_group: Point group containing information about present
        symmetry operations, necessary to correctly apply the extinction.
    :type point_group: kesshou.symmetry.Group
    """

    def _interpret_rule():
        try:
            _split_rule = [*rule.split(':', 1), '']
        except AttributeError:
            _split_rule = rule[0:2]
        _domain_string = _split_rule[0].strip(' ')
        _condition_string = _split_rule[1].strip(' ')
        return _domain_string, _condition_string

    domain_string, condition_string = _interpret_rule()
    dom = self._domain(domain_string)
    con = self._condition(condition_string)
    # obtain a part of domain which does NOT meet the condition
    dom = dom.loc[~dom.isin(con).all(1)][['h', 'k', 'l']]
    # make a set with all hkls which should be extinct
    extinct = set()
    for op in point_group.operations:
        extinct = extinct | {tuple(row) for row in op.transform(dom.to_numpy())}
    # remove all reflections whose hkls are in extinct set
    for h, k, l in list(extinct):
        self.table = self.table[~((self.table['h'] == h) &
                                  (self.table['k'] == k) &
                                  (self.table['l'] == l))]
    # reset the indices for other methods to use
    self.table.reset_index(drop=True, inplace=True)


def _domain(self, address='hkl'):
    """
    This method limits the reflection data to the ones "living" in address.

    :param address: Address which the reflections must have to be preserved
    :type address: str
    :return: Dataframe containing only reflections with given address
    :rtype: pd.DataFrame
    """

    address = address.lower().replace(' ', '').replace('_', '')
    df = self.table

    # no zeroes in address
    if address == 'hkl':
        return df
    if address in ('hhl', 'kkl'):
        return df.loc[df['h'] == df['k']]
    if address in ('hkh', 'lkl'):
        return df.loc[df['h'] == df['l']]
    if address in ('hkk', 'hll'):
        return df.loc[df['k'] == df['l']]

    # one zero in address
    if address == 'hk0':
        return df.loc[df['l'] == 0]
    if address == 'h0l':
        return df.loc[df['k'] == 0]
    if address == '0kl':
        return df.loc[df['h'] == 0]
    if address in ('hh0', 'kk0'):
        return df.loc[(df['h'] == df['k']) & (df['l'] == 0)]
    if address in ('h0h', 'l0l'):
        return df.loc[(df['h'] == df['l']) & (df['k'] == 0)]
    if address in ('0kk', '0ll'):
        return df.loc[(df['k'] == df['l']) & (df['h'] == 0)]

    # two zeroes in address
    if address == 'h00':
        return df.loc[(df['k'] == 0) & (df['l'] == 0)]
    if address == '0k0':
        return df.loc[(df['h'] == 0) & (df['l'] == 0)]
    if address == '00l':
        return df.loc[(df['h'] == 0) & (df['k'] == 0)]

    # three zeroes in address
    if address == '000':

        return df.loc[(df['h'] == 0) & (df['k'] == 0) & (df['l'] == 0)]
    # raise exception if the address is unknown
    raise ValueError('Unknown domain address have been supplied')

def _condition(self, equation=''):
    """
    This method limits the reflection data based on truth of equation

    :param equation: Equation which must be met for data to be preserved
    :type equation: str
    :return: Dataframe containing only reflections which meet the equation
    :rtype: pd.DataFrame
    """
    equation = equation.lower().replace(' ', '').replace('_', '')
    df = self.table
    # no variables return dataframe with no rows
    if equation == '':
        return df.iloc[0:0]
    # one variable
    if equation == 'h=2n':
        return df.loc[is2n(df['h'])]
    if equation == 'k=2n':
        return df.loc[is2n(df['k'])]
    if equation == 'l=2n':
        return df.loc[is2n(df['l'])]
    if equation == 'h=3n':
        return df.loc[is3n(df['h'])]
    if equation == 'k=3n':
        return df.loc[is3n(df['k'])]
    if equation == 'l=3n':
        return df.loc[is3n(df['l'])]
    if equation == 'h=4n':
        return df.loc[is4n(df['h'])]
    if equation == 'k=4n':
        return df.loc[is4n(df['k'])]
    if equation == 'l=4n':
        return df.loc[is4n(df['l'])]
    if equation == 'h=6n':
        return df.loc[is6n(df['h'])]
    if equation == 'k=6n':
        return df.loc[is6n(df['k'])]
    if equation == 'l=6n':
        return df.loc[is6n(df['l'])]
    # sum of variables
    if equation in ('h+k=2n', 'k+h=2n'):
        return df.loc[is2n(df['h'] + df['k'])]
    if equation in ('h+l=2n', 'l+h=2n'):
        return df.loc[is2n(df['h'] + df['l'])]
    if equation in ('k+l=2n', 'l+k=2n'):
        return df.loc[is2n(df['k'] + df['l'])]
    if equation in ('2h+k=4n', 'k+2h=4n'):
        return df.loc[is4n(df['h'] + df['h'] + df['k'])]
    if equation in ('h+2k=4n', '2k+h=4n'):
        return df.loc[is4n(df['h'] + df['k'] + df['k'])]
    if equation in ('2h+l=4n', 'l+2h=4n'):
        return df.loc[is4n(df['h'] + df['h'] + df['l'])]
    if equation in ('h+2l=4n', '2l+h=4n'):
        return df.loc[is4n(df['h'] + df['l'] + df['l'])]
    if equation in ('2k+l=4n', 'l+2k=4n'):
        return df.loc[is4n(df['k'] + df['k'] + df['l'])]
    if equation in ('k+2l=4n', '2l+k=4n'):
        return df.loc[is4n(df['k'] + df['l'] + df['l'])]
    if equation in ('h+k+l=2n', 'h+l+k=2n', 'k+h+l=2n',
                    'k+l+h=2n', 'l+h+k=2n', 'l+k+h=2n'):
        return df.loc[is2n(df['h'] + df['k'] + df['l'])]
    # mixed sum and multiple variables
    if equation in ('h,k=2n,h+k=4n', 'h+k=4n,h,k=2n',
                    'k,h=2n,h+k=4n', 'h+k=4n,k,h=2n',
                    'h,k=2n,k+h=4n', 'k+h=4n,h,k=2n',
                    'k,h=2n,k+h=4n', 'k+h=4n,k,h=2n'):
        return df.loc[(is2n(df['h'])) & (is2n(df['k'])) &
                      (is4n(df['h'] + df['k']))]
    if equation in ('h,l=2n,h+l=4n', 'h+l=4n,h,l=2n',
                    'l,h=2n,h+l=4n', 'h+l=4n,l,h=2n',
                    'h,l=2n,l+h=4n', 'l+h=4n,h,l=2n',
                    'l,h=2n,l+h=4n', 'l+h=4n,l,h=2n'):
        return df.loc[(is2n(df['h'])) & (is2n(df['l'])) &
                      (is4n(df['h'] + df['l']))]
    if equation in ('k,l=2n,k+l=4n', 'k+l=4n,k,l=2n',
                    'l,k=2n,k+l=4n', 'k+l=4n,l,k=2n',
                    'k,l=2n,l+k=4n', 'l+k=4n,k,l=2n',
                    'l,k=2n,l+k=4n', 'l+k=4n,l,k=2n'):
        return df.loc[(is2n(df['k'])) & (is2n(df['l'])) &
                      (is4n(df['k'] + df['l']))]
    # multiple variables
    if equation in ('h,k=2n', 'k,h=2n'):
        return df.loc[(is2n(df['h'])) & (is2n(df['k']))]
    if equation in ('h,l=2n', 'l,h=2n'):
        return df.loc[(is2n(df['h'])) & (is2n(df['l']))]
    if equation in ('k,l=2n', 'l,k=2n'):
        return df.loc[(is2n(df['k'])) & (is2n(df['l']))]
    if equation in ('h,k,l=2n', 'h,l,k=2n', 'k,h,l=2n',
                    'k,l,h=2n', 'l,h,k=2n', 'l,k,h=2n', ):
        return df.loc[(is2n(df['h'])) & (is2n(df['k'])) & (is2n(df['l']))]
    # multiple sums of variables
    if equation in ('h+k,h+l,k+l=2n', 'k+h,h+l,k+l=2n', 'h+k,l+h,k+l=2n',
                    'h+k,h+l,l+k=2n', 'k+h,l+h,k+l=2n', 'k+h,l+h,l+k=2n'):
        return df.loc[(is2n(df['h'] + df['k'])) &
                      (is2n(df['h'] + df['l'])) &
                      (is2n(df['k'] + df['l']))]
    # raise exception if the equation is unknown
    raise ValueError('Unknown condition equation have been supplied')