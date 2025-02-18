#!/usr/bin/env python3
# %module
# % description: Export the current GRASS location/mapset as a STAC-compliant Catalog JSON file including Collections, Items, and Assets.
# % keywords: export, stac, GRASS, mapset, catalog, collection, items, assets
# %end
#
# %option
# % key: output
# % type: string
# % description: Path for the output JSON file
# % required: yes
# %end

import sys


def main():
    import grass.script as gs

    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
