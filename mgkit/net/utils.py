"""
Utility functions for the network package
"""
import requests


def url_open(url, data=None, headers=None, agent=None, get=True, stream=False):
    """
    .. versionchanged:: 0.3.4
        now uses *requests*

    Arguments:
        url (str): url to request
        data (dict): parameters to pass to the request
        headers (dict): any additional headers
        agent (str): user agent to use
        get (bool): True if the request is a GET, False for POST
        stream (bool): returns an iterator to stream over

    :param str url: url to request
    :param str data: data to add to the request
    :param bool compress: if the response should be compressed
    :param str agent: if supplied, the 'User-Agent' header we'll be added to
        the request
    :return: the response handle
    """
    if (agent is not None) and (headers is None):
        headers = {'user-agent': agent}
    if get:
        request = requests.get(url, params=data, headers=headers, stream=stream)
    else:
        request = requests.post(url, params=data, headers=headers, stream=stream)

    if stream:
        return request.iter_lines()
    else:
        return request


def url_read(url, data=None, agent=None, headers=None, get=True):
    """
    .. versionchanged:: 0.3.4
        now uses *requests*, removed *compressed* and added *headers*, *get*

    Opens an URL and reads the

    Wrapper of :func:`url_open` which reads the full response

    Arguments:
        url (str): url to request
        data (dict or None): data to add to the request
        headers (dict or None): additional headers
        agent (str or None): if supplied, the 'User-Agent' header we'll be
            added to the request
        get (bool): uses a GET operation if True, POST if False

    Returns:
        str: the response data
    """
    request = url_open(url, data=data, agent=agent, headers=headers, get=get, stream=False)
    return request.text
