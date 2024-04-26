Latest changes
==============

Release 1.3.1 -- 2024/04/26
---------------------------

- Bug fix: `--heatmap-midpoint` default value changed from None to 0.5 to
  prevent `TypeError` if this value is not specified

Release 1.3.0 -- 2024/04/25
---------------------------

- Exposed heatmap colors and the midpoint of 3-color scales as commandline
  arguments so users can run wild with their creative preferences

- Joblib version 1.4.0 was officially released, so bumped dependency. Versions
  below this one are no longer supported in PhamClust

Release 1.2.0 -- 2024/03/06
---------------------------

- Added support for joblib 1.4.dev0 (installed from GitHub), unlocking
  ``Parallel(return_as="generator_unordered")``. This restores pre-joblib speed
  while maintaining the very low memory footprint enabled by joblib.

- Matrices can now be sorted in tree order; call ``.reorder()`` on your matrix
  without supplying a node-order and it will be put in tree order.

- PDF heatmaps traded for SVG format instead.

Release 1.1.0 -- 2024/02/08
---------------------------

- Replaced internal parallelization with joblib.Parallel for drastically
  reduced memory footprint. This results in approximately 10-20% reduction in
  speed for the time being - see `In development`.

- Refactored `__main__.py` to separate the setup steps from phamclust logic
