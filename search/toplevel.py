"""Top-level initialization for search."""

from .template import build_templates

def initialize():
    """Builds Django templates and writes them to the relevant directory on
    disk.
    """
    print("Building templates...")
    for (template, _) in build_templates():
        print("Wrote {}.".format(template))
    print()
