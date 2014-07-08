"""
Dictionary utils

"""
import numpy


def combine_dict(keydict, valuedict):
    """
    Combine two dictionaries when the values of keydict are iterables. The
    combined dictionary has the same keys as keydict and the its values are sets
    containing all the values associated to keydict values in valuedict.

    .. digraph:: keydict
        :alt: key1 -> [v1, v2, .., vN]

        "keydict" -> "key1" -> "[v1, v2, .., vN]";

    .. digraph:: valuedict
        :alt: v1 -> [u1, u2, .., uN]
              v2 -> [t1, t2, .., tN]

        "valuedict" -> "v1" -> "[u1, u2, .., uN]";
        "valuedict" -> "v2" -> "[t1, t2, .., tN]";

    Resulting dictionary will be

    .. digraph:: combined
        :alt: key1->{u1, u2, .., uN}

        "combined" -> "key1" -> "{u1, u2, .., uN, t1, t2, .., tN}";

    :param dict keydict: dictionary whose keys are the same as the returned
        dictionary
    :param dict valuedict: dictionary whose values are the same as the returned
        dictionary

    :return dict: combined dictionary
    """

    comb_dict = dict((key, set()) for key in keydict)

    for key, values in keydict.items():
        for value in values:
            try:
                comb_dict[key].update(valuedict[value])
            except KeyError:
                #in case a value isn't in valuedict keys, silently pass
                pass

    return comb_dict


def combine_dict_one_value(keydict, valuedict):
    """
    Combine two dictionaries by the value of the keydict is used as a key in
    valuedict and the resulting dictionary is composed of keydict keys and
    valuedict values.

    Same as :func:`comb_dict`, but each value in keydict is a single element
    that is key in valuedict.

    :param dict keydict: dictionary whose keys are the same as the returned
        dictionary
    :param dict valuedict: dictionary whose values are the same as the returned
        dictionary

    :return dict: combined dictionary
    """

    comb_dict = {}

    for key, value in keydict.items():
        comb_dict[key] = valuedict[value]

    return comb_dict


def link_ids(id_map, black_list=None):
    """
    Given a dictionary whose values (iterables) can be linked back to other
    keys, it returns a dictionary in which the keys are the original keys and
    the values are sets of keys to which they can be linked.

    .. digraph:: input
        :alt: key1->[v1, v2]
              key2->[v3, v4]
              key3->[v2, v4]

        "id_map" -> "key1" -> "[v1, v2]";
        "id_map" -> "key2" -> "[v3, v4]";
        "id_map" -> "key3" -> "[v2, v4]";

    Becomes:

    .. digraph:: output
        :alt: key1->[key1, key3]
              key2->[key3]
              key3->[key1, key2]

        "linked" -> "key1" -> "[key1, key3]"
        "linked" -> "key2" -> "[key3]"
        "linked" -> "key3" -> "[key1, key2]"

    :param dict id_map: dictionary of keys to link
    :param iterable black_list: iterable of values to skip in making the links

    :return dict: linked dictionary
    """
    id_to_id = {}

    for s_id, s_cps in id_map.iteritems():
        id_links = set()
        for e_id, e_cps in id_map.iteritems():
            if e_id == s_id:
                continue
            for s_id2 in s_cps:
                if not black_list is None:
                    if s_id2 in black_list:
                        continue
                if s_id2 in e_cps:
                    id_links.add(e_id)
                    break
        id_to_id[s_id] = id_links
    return id_to_id


def reverse_mapping(map_dict):
    """
    Given a dictionary in the form: key->[v1, v2, .., vN], returns a dictionary
    in the form: v1->[key1, key2, .., keyN]

    :param dict map_dict: dictionary to reverse

    :return dict: reversed dictionary
    """
    rev_map = {}

    for key, value_ids in map_dict.iteritems():
        for value_id in value_ids:
            try:
                rev_map[value_id].add(key)
            except KeyError:
                rev_map[value_id] = set([key])

    return rev_map


def find_id_in_dict(s_id, s_dict):
    """
    Finds a value 's_id' in a dictionary in which the values are iterables.
    Returns a list of keys that contain the value.

    :param dict s_id: element to look for in the dictionary's values
    :param object d: dictionary to search in

    :return list: list of keys in which d was found
    """
    f_list = []
    for k_id, v_ids in s_dict.items():
        if s_id in v_ids:
            f_list.append(k_id)
    return f_list


def split_dictionary_by_value(value_dict, threshold, aggr_func=numpy.median,
                              key_filter=None):
    """
    Splits a dictionary, whose values are iterables, based on a threshold:

        * one in which the result of aggr_func is lower than the threshold
          (first)
        * one in which the result of aggr_func is equal or greater than the
          threshold (second)

    :param dict valuedict: dictionary to be splitted
    :param number threshold: must be comparable to threshold
    :param func aggr_func: function used to aggregate the dictionary values
    :param iterable key_filter: if specified, only these key will be in the
        resulting dictionary

    :return: two dictionaries
    """
    lower_dict = {}
    higher_dict = {}

    if key_filter is None:
        key_filter = value_dict.iterkeys()

    for key in key_filter:
        values = value_dict[key]
        if aggr_func(values) < threshold:
            lower_dict[key] = values
        else:
            higher_dict[key] = values

    return lower_dict, higher_dict


def apply_func_to_values(dictionary, func):
    """
    .. versionadded:: 0.1.12

    Assuming a dictionary whose values are iterables, *func* is applied to each
    element of the iterable, retuning a *set* of all transformed elements.

    Arguments:
        dictionary (dict): dictionary whose values are iterables
        func (func): function to apply to the dictionary values

    Returns:
        dict: dictionary with transformed values
    """
    return dict(
        (key, set(func(value) for value in values))
        for key, values in dictionary.iteritems()
    )


def filter_ratios_by_numbers(ratios, min_num):
    """
    Returns from a dictionary only the items for which the length of the
    iterables that is the value of the item, is equal or greater of min_num.

    :param dict ratios: dictionary key->list
    :param int min_num: minimum number of elements in the value iterable

    :return dict: filtered dictionary
    """
    return dict(
        (key, values) for key, values in ratios.iteritems()
        if len(values) >= min_num
    )


def filter_nan(ratios):
    """
    Returns a dictionary with the NaN values taken out
    """

    return dict(
        (key, [ratio for ratio in ratios[key] if not numpy.isnan(ratio)])
        for key in ratios
    )
