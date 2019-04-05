import hy
import os

from django.urls import path

from . import views
from .toplevel import initialize

## urls.py is evaluated once upon initialization, so any functions that need to be
## called from the top-level are called here.
initialize()

urlpatterns = [
    path('', views.query, name='query'),
]
