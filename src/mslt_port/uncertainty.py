"""Add support for uncertainty analyses by varying input data."""

import pandas as pd
import numpy as np
import scipy.stats as st


class Normal:
    """
    Draw samples from a normal distribution, whose standard deviation is a
    percentage of the mean.

    :param sd_pcnt: The standard deviation, represented as a percentage of the
        mean; valid values are between ``0.1`` and ``20``.
    """

    _sd_min = 0.1
    _sd_max = 20

    def __init__(self, sd_pcnt):
        if sd_pcnt < self._sd_min or sd_pcnt > self._sd_max:
            raise ValueError('SD% of {} is not between {}% and {}%'.format(
                sd_pcnt, self._sd_min, self._sd_max))
        self.sd_frac = sd_pcnt / 100.0

    def uncorrelated_samples(self, mean_values, num_samples, prng):
        """
        Draw **uncorrelated** values at random.

        :param mean_values: The mean values.
        :param num_samples: The number of samples to return.
        :param prng: A ``numpy.random.RandomState`` instance.
        :returns: An array of shape ``(n, v)`` for ``n`` samples and ``v``
            mean values.
        """
        mean = mean_values.values
        sd = np.absolute(mean) * self.sd_frac
        size = (num_samples, ) + mean.shape
        return prng.normal(mean[np.newaxis, ...],
                           sd[np.newaxis, ...],
                           size=size)

    def correlated_samples(self, mean_values, samples):
        """
        Draw **correlated** values according to the provided percentiles.

        :param mean_values: The mean values.
        :param samples: The percentiles.
        :returns: An array of shape ``(n, v)`` for ``n`` percentile samples
            and ``v`` mean values.
        """
        mean = mean_values.values
        zero_mask = mean == 0.0
        sd = np.absolute(mean) * self.sd_frac
        if np.any(zero_mask):
            # Only draw samples where the SD is non-zero.
            values = np.zeros(samples.shape + mean.shape)
            nonz_mask = ~ zero_mask
            rv = st.norm(loc=mean[nonz_mask], scale=sd[nonz_mask])
            values[:, nonz_mask] = rv.ppf(samples[..., np.newaxis])
            return values
        else:
            rv = st.norm(loc=mean, scale=sd)
            return rv.ppf(samples[..., np.newaxis])


class Beta:
    """
    Draw samples from a Beta distribution, whose standard deviation is a
    percentage of the mean.

    :param sd_pcnt: The standard deviation, represented as a percentage of the
        mean; valid values are between ``0.1`` and ``20``.
    """

    _sd_min = 0.1
    _sd_max = 20

    def __init__(self, sd_pcnt):
        if sd_pcnt < self._sd_min or sd_pcnt > self._sd_max:
            raise ValueError('SD% of {} is not between {}% and {}%'.format(
                sd_pcnt, self._sd_min, self._sd_max))
        self.sd_frac = sd_pcnt / 100.0

    def shape_parameters(self, mean_values):
        """Calculate the shape parameters alpha and beta."""
        mean = mean_values.values
        sd = np.absolute(mean) * self.sd_frac
        varn = np.square(sd)
        nu = mean * (1 - mean) / varn - 1
        alpha = mean * nu
        beta = (1 - mean) * nu
        return (alpha, beta)


    def uncorrelated_samples(self, mean_values, num_samples, prng):
        """
        Draw **uncorrelated** values at random.

        :param mean_values: The mean values.
        :param num_samples: The number of samples to return.
        :param prng: A ``numpy.random.RandomState`` instance.
        :returns: An array of shape ``(n, v)`` for ``n`` samples and ``v``
            mean values.
        """
        size = (num_samples, ) + mean_values.shape
        alpha, beta = self.shape_parameters(mean_values)
        return prng.beta(a=alpha[np.newaxis, ...],
                         b=beta[np.newaxis, ...],
                         size=size)

    def correlated_samples(self, mean_values, samples):
        """
        Draw **correlated** values according to the provided percentiles.

        :param mean_values: The mean values.
        :param samples: The percentiles.
        :returns: An array of shape ``(n, v)`` for ``n`` percentile samples
            and ``v`` mean values.
        """
        zero_mask = mean_values.values == 0.0
        if np.any(zero_mask):
            # Only draw samples where the SD is non-zero.
            values = np.zeros(samples.shape + mean_values.shape)
            nonz_mask = ~ zero_mask
            alpha, beta = self.shape_parameters(mean_values[nonz_mask])
            rv = st.beta(a=alpha, b=beta)
            values[:, nonz_mask] = rv.ppf(samples[..., np.newaxis])
            return values
        else:
            alpha, beta = self.shape_parameters(mean_values)
            rv = st.beta(a=alpha, b=beta)
            return rv.ppf(samples[..., np.newaxis])


def test(num=10):
    """Select samples using both methods."""
    prng = np.random.RandomState(seed=49430)

    sd_pcnt = 5
    n = Normal(sd_pcnt)
    print('Normal')

    x = pd.Series(range(1, 6))
    print('x = {}'.format(x.values))
    y = n.uncorrelated_samples(x, num, prng)
    print(y)
    print()

    samples = prng.random_sample(size=num)
    print('percentiles = {}'.format(samples))
    y = n.correlated_samples(x, samples)
    print(y)
    print()

    b = Beta(sd_pcnt)
    print('Beta')

    x = pd.Series([0.1, 0.3, 0.5])
    print('x = {}'.format(x.values))
    y = b.uncorrelated_samples(x, num, prng)
    print(y)
    print()

    samples = prng.random_sample(size=num)
    print('percentiles = {}'.format(samples))
    y = b.correlated_samples(x, samples)
    print(y)
    print()

    # Construct a table by adding each draw column in turn.
    df = pd.DataFrame(x, columns=['x'])
    for ix, col_values in enumerate(y):
        col_name = 'draw_{}'.format(ix)
        df[col_name] = col_values
    print(df.to_string())
    print()

    # Construct a table by concatenating all draw columns in one action.
    df2 = pd.DataFrame(x, columns=['x'])
    df_draws = pd.DataFrame(np.transpose(y))
    df_draws.columns = ['draw_{}'.format(ix) for ix in range(len(y))]
    df2 = pd.concat([df2, df_draws], axis=1)
    print(df2.to_string())
    print()
