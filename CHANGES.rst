Latest changes
==============

In development
--------------

- Awaiting next joblib release, to unlock
  `Parallel(return_as="generator_unordered")`. This should restore previously
  established speed while maintaining low memory footprint.

- Add distogram to heatmaps

- Allow sorting of clusters and subclusters by distance-tree order rather than
  descending size order

Release 1.1.0 -- 2023/02/08
---------------------------

- Replaced internal parallelization with joblib.Parallel for drastically
  reduced memory footprint. This results in approximately 10-20% reduction in
  speed for the time being - see `In development`.

- Refactored `__main__.py` to separate the setup steps from phamclust logic