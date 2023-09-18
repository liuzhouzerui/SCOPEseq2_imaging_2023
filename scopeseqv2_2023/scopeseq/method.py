import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.stats import pearsonr


def register(position, reference, d_th, doublet_removal=False):
    """
    find the nearest position in another dataset
    :param position: data points [X,Y]
    :param reference: reference data points [X,Y]
    :param d_th: the distance threshold of near
    :param doublet_removal: if match to multiple reference set, then decided as not registered.
    :return:
        d_position: the nearest ref data point num, size=position
        d_inrange: whether in the distance range to the nearest ref data point, size=position
    """
    d_tmp = distance.cdist(position, reference, 'euclidean')
    d_position = np.argmin(d_tmp, axis=1)
    d_value = np.min(d_tmp, axis=1)
    if doublet_removal:
        double = (d_tmp < d_th).sum(axis=1) > 1
        d_value[double] = d_th
    d_inrange = (d_value < d_th)
    return d_position, d_inrange


def find_rotation_matrix(x, y):
    """
    :param x: 2d vector
    :param y: 2d vector
    :return: the rotation matrix from y to x
    """
    x = x / np.sqrt(np.sum(x**2))
    y = y / np.sqrt(np.sum(y**2))

    cos_xy = x[0]*y[0]+x[1]*y[1]
    sin_xy = x[0]*y[1]-x[1]*y[0]

    rotation_matrix = np.array([[cos_xy, - sin_xy], [sin_xy, cos_xy]])

    return rotation_matrix


def rotate_position(position, target_start, target_end, initial_start, initial_end):
    target_vector = target_end - target_start
    initial_vector = initial_end - initial_start

    # get rotation matrix
    rotation_matrix = find_rotation_matrix(target_vector, initial_vector)

    rotated_position = target_start + np.dot(np.array(position - initial_start),
                                             rotation_matrix) * target_vector / np.dot(initial_vector,
                                                                                       rotation_matrix)

    # return rotated position
    return rotated_position


def quantile_linear_transformation(x, round):
    if round in [1, 2, 3, 4, 5, 8]:
        q1 = 0.25
        q2 = 0.75
    if round in [6, 7]:
        q1 = 0.33
        q2 = 0.83
    marker_1 = x.quantile(q1)
    marker_2 = x.quantile(q2)
    x_new = (x-marker_1)/(marker_2-marker_1)*14883+1500
    x_new[x_new < 1500] = 1500
    x_new[x_new > 16383] = 16383
    return x_new.values, marker_1, marker_2


def assign_obc(x, barcode_ref, no_signal_th=None, mode='all'):
    """
    assign optical bead barcode using bead-by-bead method. obc num is 0 based
    :param x: probe intensity vector
    :param barcode_ref: reference obc pool
    :param no_signal_th: threshold for no signal
    :param mode: 'all' for iteration, 'max' for only first round
    :return: obc vector. or -1 for unmatched obc
    """
    if min(x) == -1:
        return -1, -1
    if (no_signal_th is not None) and (max(x) < no_signal_th):
        return -1, -1
    x_sorted = np.sort(x)
    x_dif = np.diff(x_sorted)/x_sorted[0:(len(x)-1)]
    x_th = x_sorted[np.argmax(x_dif)]
    barcode = (x > x_th) * 1
    i = 1
    if mode == 'all':
        while i < (len(x)-1):
            if ''.join(barcode.astype("str")) in barcode_ref.values:
                return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0], i
            else:
                x_dif[np.argmax(x_dif)] = 0
                x_th = x_sorted[np.argmax(x_dif)]
                barcode = (x > x_th) * 1
                i = i + 1
        if ''.join(barcode.astype("str")) in barcode_ref.values:
            return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0], i
        else:
            return -1, i
    if mode == 'max':
        if ''.join(barcode.astype("str")) in barcode_ref.values:
            return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0], i
        else:
            return -1, i


def assign_obc_stage(x, barcode_ref):
    """
    assign optical bead barcode using cycle-by-cycle method
    :param x:
    :param barcode_ref:
    :return:
    """
    if min(x) == -1:
        return -1
    if ''.join(x.astype("str")) in barcode_ref.values:
        return np.where(barcode_ref.values == ''.join(x.astype("str")))[0][0]
    else:
        return -1


def obc_correlation(x, barcode_ref):
    if min(x) == -1:
        return pd.Series(-1, index=barcode_ref.index)
    else:
        correlation = barcode_ref.apply(lambda y: pearsonr(x, y)[0], axis=1)
        return correlation


def patch_index_replace(obc, replace):
    print('map image patch to sequencing i7 index...')

    # check params
    if not isinstance(replace, dict):
        raise TypeError(f'replace should be dict type.')

    patch = obc.apply(lambda x: int(x.split('_')[0]))
    s = obc.apply(lambda x: x.split('_')[1])
    q = obc.apply(lambda x: x.split('_')[2])

    seq_index = patch.copy()
    for k, v in replace.items():
        seq_index[patch == k] = v

    obc_new = seq_index.astype('str') + '_' + s + '_' + q
    return obc_new


def assign_cluster(data, thres):
    data['image_cluster'] = ''
    c = 0
    for i in thres.keys():
        pos = data[i].values > thres[i]
        neg = data[i].values <= thres[i]
        if c > 0:
            data.loc[pos, 'image_cluster'] = data.loc[pos, 'image_cluster'] + '_pos'
            data.loc[neg, 'image_cluster'] = data.loc[neg, 'image_cluster'] + '_neg'
        else:
            data.loc[pos, 'image_cluster'] = 'pos'
            data.loc[neg, 'image_cluster'] = 'neg'
        c += 1

    return data
