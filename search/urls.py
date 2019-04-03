import hy
import os

from django.urls import path

from . import views
from .query_page import build_query_page

# Top-level initialization.

print("Building templates...")
query_page_path = os.path.join(os.path.dirname(__file__), "templates/search/index.html")
with open(query_page_path, "w+") as query_page:
    query_page.write(build_query_page())
    print("Wrote {}".format(query_page_path))
print()

urlpatterns = [
    path('', views.query, name='query'),
]
