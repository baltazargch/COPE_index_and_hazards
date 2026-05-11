#!/usr/bin/env python3

import numpy as np
import numba as nu

@nu.jit(nopython=True, parallel=False)
def __bioclimvars_nu(prec, tmin, tmax, tmp, wettest_quater, driest_quater, wettest_month, driest_month, warmest_quater,
                   coldest_quater, warmest_month, coldest_month, tmp_reshaped, quaters):

    bio1 = np.mean(tmp)
    bio2 = np.mean(tmax - tmin)
    bio4 = np.std(tmp) * 100.
    bio5 = tmax[warmest_month]
    bio6 = tmin[coldest_month]
    bio7 = bio5 - bio6
    bio3 = (bio2 / bio7) * 100.
    bio8 = np.mean(tmp_reshaped[wettest_quater, :])
    bio9 = np.mean(tmp[quaters[driest_quater]])
    bio10 = np.mean(tmp[quaters[warmest_quater]])
    bio11 = np.mean(tmp[quaters[coldest_quater]])
    bio12 = np.sum(prec)
    bio13 = prec[wettest_month]
    bio14 = prec[driest_month]

    prec_mean = np.mean(prec)
    bio15 = 0.
    if prec_mean != 0:
        bio15 = np.std(prec) / prec_mean

    bio16 = np.sum(prec[quaters[wettest_quater]])
    bio17 = np.sum(prec[quaters[driest_quater]])
    bio18 = np.sum(prec[quaters[warmest_quater]])
    bio19 = np.sum(prec[quaters[coldest_quater]])

    return (bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10,
            bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)


def biovars(pre, tmp, tmn, tmx):
    '''
    pre: numpy array of monthly precipitations
    tmp: numpy array of average monthly temperatures
    tmn: numpy array of minimum monthly temperatures
    tmx: numpy array of maximum monthly temperatures

    All arrays must have length 12. Otherwise ValueError exception is raised

    return: tuple with 19 values of bioclimatic variables BIO1-BIO19
    '''

    quaters = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])


    pre_quart_means = np.mean(pre.reshape(4, 3), axis=1)
    wettest_quater = np.argmax(pre_quart_means)
    driest_quater = np.argmin(pre_quart_means)
    wettest_month = np.argmax(pre)
    driest_month = np.argmin(pre)

    tmp_reshaped = tmp.reshape(4, 3)
    tmp_quart_means = np.mean(tmp.reshape(4, 3), axis=1)
    warmest_quater = np.argmax(tmp_quart_means)
    coldest_quater = np.argmin(tmp_quart_means)
    warmest_month = np.argmax(tmp)
    coldest_month = np.argmin(tmp)

    res_bio = __bioclimvars_nu(pre, tmn, tmx, tmp, wettest_quater, driest_quater, wettest_month, driest_month,
                             warmest_quater, coldest_quater, warmest_month, coldest_month, tmp_reshaped, quaters)

    return res_bio

if __name__ == "__main__":
    print ("Test")
    pre = np.array([9., 2., 3.2, 8.,
                    9., 2., 3.3, 14.2,
                    9., 3.7, 2.2, 1.4])
    tmn = np.array([-2.5, -3.1, 0.3, 2.4,
                    6.2, 19.3, 22.4, 20.1,
                    12.1, 7.2, 0.1, -4.1])
    tmx = tmn + 3.9
    tmp = (tmn + tmx) / 2.

    bio = biovars(pre, tmp, tmn, tmx )
    print (bio)
