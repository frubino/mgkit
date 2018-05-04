import pytest
import requests

try:
    requests.get('http://www.google.com')
    conn_ko = False
except requests.exceptions.ConnectionError:
    conn_ko = True


skipif = pytest.mark.skipif(conn_ko, reason='No connection available')
