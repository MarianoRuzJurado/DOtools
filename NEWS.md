# DOtools 0.4.0

All notable changes to this project will be documented in this file.

---

## [0.4.0] - 2025-06-09

### Added
* `DO.SplitBarGSEA`: Bar plot function for GSEA results (e.g., from Metascape).
* `DO.scVI`: Wrapper function to run scVI integration on a Seurat object.

### Changed
* Minor adjustments across various functions for improved consistency.

---

## [0.3.0] - 2025-05-28

### Added
* `DO.Import`: Simplifies Seurat object creation.
* `DO.CellBender`: Enables running CellBender on aligner-generated HDF5 files.

### Changed
* Minor adjustments and internal improvements.

---

## [0.2.0] - 2025-05-26

### Added
* `DO.CellBender`: Functionality to run CellBender (initial version).
* `DO.BarcodeRanks`: New function for barcode rank visualization.

### Changed
* Minor adjustments for stability and usability.

---

## [0.1.1] - 2025-05-26

### Changed
* Minor adjustments to function behavior and formatting.

---

## [0.1.0] - 2025-04-09

### Added
* `DO.CellCompositions`: Integrates the scanpro tool to detect changes in cell composition.
* New connected stacked bar plot in scanpro output highlighting significant composition changes.
* `DO.Dotplot`: Added pseudobulk option for y-axis.

---

## [0.0.2] - 2025-04-03

### Fixed
* Bug in `DO.Box.plot`: Meta.data features could not be specified.
* Bug in `DO.Vln.plot`: Meta.data features could not be specified.

### Added
* New argument in `DO.Box.plot`: Enables user-defined x-axis sorting.

---

## [0.0.1] - 2025-03-27

### Added
* Initial commit: core structure and foundational functions.
