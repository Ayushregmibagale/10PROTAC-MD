# Exact public subset of a recorded execution snapshot

The manifest retains every executed source/config hash. Separately licensed
Rosetta helper source is omitted from Git; obtain its exact pinned files with
`python baseline_engine/fetch_rosetta_helpers.py` in this directory in an
appropriately licensed environment. The official helpers' source and hashes are
recorded in vendor/SOURCE.json. No PyRosetta binary/database is included.
This snapshot preserves the precise five-compound validation code even as the
main working tree advances. Complete raw artifacts remain in the local project.
