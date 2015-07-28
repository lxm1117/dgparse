import pytest

import collections
import os

from contextlib import contextmanager


@pytest.fixture(scope='module')
def fixture_data():
    'Load fixture data, relative to tests dir'
    @contextmanager
    def fixture_data_fn(path, mode='rb'):
        test_root = os.path.dirname(__file__)
        handle = open(os.path.join(test_root, path), mode)
        yield handle
        handle.close()
    return fixture_data_fn


@pytest.fixture(scope='module')
def utils(fixture_data):
    attrs_dict = {
        'fixture_data': fixture_data,
    }
    return collections.namedtuple('utils', attrs_dict.keys())(**attrs_dict)
