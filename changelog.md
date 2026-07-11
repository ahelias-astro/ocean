# Changelog

## Maintained by Adrien Hélias, Western University, Canada

### v1.0.1 (July 11, 2026)

Bug fixes:
- The function slepian2 was indexing into the eigenvector result (evs) in a way that was returning a scalar instead of a 1D vector, as initially intended. This could produce much smaller values than expected. It is now fixed.
- Expanded the masking of the light curve values to time and magnitude errors in the three sigma filter feature to prevent code errors.

### v1.0.0 (May 16, 2026)

Initial release of *ocean*.

