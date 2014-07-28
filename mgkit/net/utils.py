"""
Utility functions for the network package
"""

import urllib2
import gzip
import cStringIO
import sys


def url_open(url, data=None, compress=True, agent=None):
    """
    Utility function that compresses the request and add the user agent header
    if supplied.

    :param str url: url to request
    :param str data: data to add to the request
    :param bool compress: if the response should be compressed
    :param str agent: if supplied, the 'User-Agent' header we'll be added to
        the request
    :return: the response handle
    """
    if data is None:
        request = urllib2.Request(url)
    else:
        request = urllib2.Request(url, data)

    if agent is not None:
        request.add_header(
            'User-Agent', "Python {0}.{1} - {2}".format(
                sys.version_info[0],
                sys.version_info[1],
                agent
            )
        )

    if compress:
        request.add_header('Accept-encoding', 'gzip')

    response = urllib2.urlopen(request)

    if response.info().get('Content-Encoding') == 'gzip':
        data_handle = cStringIO.StringIO(response.read())
        data_handle = gzip.GzipFile(fileobj=data_handle, mode='rb')
    else:
        data_handle = response

    return data_handle


def url_read(url, data=None, compress=True, agent=None):
    """
    Utility function that compresses the request and add the user agent header
    if supplied.

    Wrapper of :func:`url_open` which reads the full response

    :param str url: url to request
    :param str data: data to add to the request
    :param bool compress: if the response should be compressed
    :param str agent: if supplied, the 'User-Agent' header we'll be added to
        the request
    :return: the response data
    """
    return url_open(url, data, compress, agent).read()
