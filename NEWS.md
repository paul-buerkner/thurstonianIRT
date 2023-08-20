# thurstonianIRT 0.12.3

* Avoid a CRAN check warning.

# thurstonianIRT 0.12.2

* Fix line-length issues in Mplus code generation. (#34)

# thurstonianIRT 0.12.1

* Fix data simulations via `sim_TIRT_data` for non-triplet formats.
* Ensure compatibility with both rstan versions 2.21 and 2.26.

# thurstonianIRT 0.12.0

* Switch to the latest `rstantools` toolchain.
* Add a `gof` method to compute goodness-of-fit measures
for lavaan-fitted models thanks to Angus Hughes. (#19)


# thurstonianIRT 0.11.1

* Fix tests failing after an update of `lavaan`. (#25)


# thurstonianIRT 0.11.0

## Bug Fixes

* Fix usage of `gamma` parameters in `sim_TIRT_data`
thanks to @IanEisenberg. (#13)
* Prevent impossible rankings from being sampled in 
`sim_TIRT_data` thanks to Susanne Frick.

## New Features

* Support predictions of trait scores for new persons. (#12, #15)
* Specify multiple blocks at once via `set_blocks_from_df`
thanks to @awhug. (#11)


# thurstonianIRT 0.10.0

## Bug Fixes

* Correctly handle standard errors of trait scores as returned by Mplus.
* Fix the `triplets` example data to work with lavaan and Mplus without
convergence issues. (#3)

## New Features

* Support family `gaussian` when using lavaan as the model fitting engine.

## Other Changes

* Improve the documentation and presentation of the package across the board 
thanks to Russell S. Pierce, Thomas J. Faulkenberry, and Daniel S. Katz.
(#2, #5, #7, #8, #9)


# thurstonianIRT 0.9.0
  
* Initial release version
