import hy
import os

from django.urls import path

from . import views
from .template import build_templates

# Top-level initialization.

print("Building templates...")
for (template, _) in build_templates():
    print("Wrote {}.".format(template))
print()

urlpatterns = [
    path('', views.query, name='query'),
]
