Latest changes
==============

In development
--------------

- Add distogram to heatmaps

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
