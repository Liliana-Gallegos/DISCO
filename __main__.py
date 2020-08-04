#!/usr/bin/python
from __future__ import absolute_import
import os, sys

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ == '':
    path = os.path.dirname(os.path.dirname(__file__))
    sys.path.insert(0, path)

from disco import DISCO  # noqa

if __name__ == '__main__':
    sys.exit(DISCO.main())
