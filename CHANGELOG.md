# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project should adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2023-08-13

### Changed

- Updated computation of maximum traffic density to utilize Matlab perforamnce improvements since original release: no longer requires compiling executable if using Matlab version R2020b or newer.

### Fixed

- Radar coverage data files corrected (`cellcoverage.mat` and `cellcoveragelowh.mat` on Zenodo): previously, some radars were erroneously included or excluded. 240 radar sites remained, 10 were removed, and 36 were added.

## [1.0] - 2020-12-15

### Added

- Initial public release

[1.1]: https://github.com/mit-ll/traffic-density-database/releases/tag/v1.1
[1.0]: https://github.com/mit-ll/traffic-density-database/releases/tag/v1.0
