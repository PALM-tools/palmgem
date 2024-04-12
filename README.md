PALM-GeM: Geospatial Data Merging and preprocessing into PALM
=============================================================

## General TODO:
* finish documentation
* Add filtration routines
* Merge with codes from our palm_inputs
* Merge slanted faces (in production branch, with own docs, config in branch, will be merged latter)
* test that all works fine
* finish usercases - bergen/brno/prague
* In user cases add some shell script that would copy example configs into config, import to import etc, so the folders are not messy
* push to gitlab, send note to gitlab
* finish manuscript

A versatile tool for processing geospatial data from various origins and
preparing static drivers for the PALM modelling system.

## Opensource PALM Static Driver preprocessor
Short paragraph about tool
* All mayor cities in Europe are there
* Note about how user can prepare static driver for every bigger european city and with other PALM tools (WRF interface, PROMET, PALM itself) can run PALM simulation with reasonable good resolution.

## Installation / Requirements
Detailed installation can be found here: [Installation](docs/install.md)
TODO: Jarda please write installation and requirements for PostgreSQL a PostGIS.\
TODO: python version + version of libraries

## Documentation
General description of used modules, function etc.. [docs](docs/general.md) \
Detailed list of all available configuration are listed here: [configuration list](docs/configuration_docs.md) \
Running preprocessor python script and its configuration is documented here: [preprocesor run](docs/run_preprocessor.md) \
The main run of static driver generator with its configuration is documented here: [static driver generator](docs/run_palm_static_driver.md)
[Tutorial](docs/visuallization.md) how to connection PostgreSQL database with QGIS application and enable interactive visualization.

## Contributing / Support
If you wish to use your own dataset with more detailed parameterization, please contact us. We will gladly help with implementation.

## Change log
Current version developed to process open-source datase (UrbanAtlas, OpenStreetMaps).
2024.03.19 - This is an initial commit, more details will be added later.

## Example
We have prepared several example testcases \
[Brno testcase](examples/brno/README.md) \
[Tromso testcase](examples/tromso/README.md)

## Benchmark
| Domain | Size [km x km] | Grid size | Time [h] |
|:------------------|:---------------|:---------------|:----------|
| Brno | 0.640 x 0.640 | 64 x 64 | 12 [s] |
| Brno | 1.280 x 1.280 | 128 x 128 | 23 [s] |
| Brno | 2.560 x 2.560 | 256 x 256 | 117 [s] |
| Brno | 5.120 x 5.120 | 512 x 512 | 14.3 [min] |
| Brno | 10.240 x 10.240 | 1024 x 1024 |  [min] |
| Brno | 20.480 x 20.480 | 2048 x 2048 |  [min] |
| Brno | 40.960 x 40.960 | 4096 x 4096 |  [min] |
| Prague |  x  |  x  |  [min] |
| Bergen | 0.640 x 0.64 | 64 x 64 |  88 [s] |
| Bergen | 1.280 x 1.280 | 128 x 128 | 100  [s] |
| Bergen | 2.560 x 2.560 | 256 x 256 | 155 [s] |
| Bergen | 5.120 x 5.120 | 512 x 512 | 11 [min] |
| Bergen | 10.240 x 10.240 | 1024 x 1024 |  [min] |
| Bergen | 20.480 x 20.480 | 2048 x 2048 | 20 [h] |
| Berlin | 0.320 x 0.320 | 64 x 64 | 13 [s] |
| Berlin | 0.640 x 0.640 | 128 x 128 | 134 [s] |
| Berlin | 1.280 x 1.280 | 256 x 256 | 123 [s] |
| Berlin | 2.560 x 2.560 | 512 x 512 | 30 [min] |
| Berlin | 4.000 x 4.000 | 800 x 800 | 67 [min] |

## Example static driver

## Project status & Future versions
Currently, we are able to produce PALM's static driver for most of the bigger cities in EU. We have developed (in testing branch: ....) an extension for cut cell topography (link to PALM CCT or article). The cut cell tool is under development, but the version is working in most cases. There is an extension under development that would process finer geospatial datasets into PALM finer parametrization (vegetation pars, building pars, building surface pars, etc.). This extension would also process individual trees location into LAD (leaf area density) variable. Examples of this extensions can be found: cite our work. The extension will be soon part of main branch in PALM-GEM.

## License


## Authors
List of authors with contacts.

## Cite us