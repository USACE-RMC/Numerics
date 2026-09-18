# Independent generalized-family affine-overflow tail references

Generated on 2026-09-08 with **R 4.4.3 (2025-02-28 ucrt)**, platform
`x86_64-w64-mingw32`, locale `C`, using `generalized-affine-overflow.R`.
The standalone script uses base R only and never loads Numerics. It emits the
runtime metadata and four reference rows; it does not alter other frozen fixtures.
R emitted four startup locale-setting warnings but exited successfully with all
numerical checks passing.

All cases use location zero and scale `1E-200`. At observation `1E200` and shape
`-20`, the ratio `x/alpha` overflows binary64, but the defining Hosking latent value
is finite:

```
logSupport = log(20) + log(abs(x)) - log(alpha)
             + log1p(exp(-(log(20) + log(abs(x)) - log(alpha))))
z = -logSupport/kappa = 46.201488473558619
```

Reflecting both x and shape negates z and gives the identical opposite-tail and
density limits. The first omitted ordinary ratio is never formed. Normal tails
and density use R `pnorm(..., log.p=TRUE)` and `dnorm(..., log=TRUE)`. Logistic
values use their defining log formulas, independently cross-checked against
R `plogis` and `dlogis`. The density Jacobian is `kappa*z-log(alpha)`.

| Family | x | kappa | Tail | Log tail | Log PDF | PDF |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| GNO | -1E200 | 20 | lower | -1072.0411870643629 | -1531.7204579917529 | 0 |
| GNO | 1E200 | -20 | upper | -1072.0411870643629 | -1531.7204579917529 | 0 |
| GLO | -1E200 | 20 | lower | -46.201488473558619 | -509.71423934592178 | 4.3044582966584941E-222 |
| GLO | 1E200 | -20 | upper | -46.201488473558619 | -509.71423934592178 | 4.3044582966584941E-222 |

The GNO ordinary PDF and tail genuinely underflow; their logarithms remain finite.
Both GLO ordinary tail and PDF remain representable. Tests retain a `5E-12`
absolute allowance for log values and `3E-12` relative allowance for ordinary
values, with exact zero required for underflow. References are R binary64
arithmetic, not arbitrary-precision claims. Existing compensated near-support
transforms remain separately covered by the accepted adjacent-boundary fixtures.
