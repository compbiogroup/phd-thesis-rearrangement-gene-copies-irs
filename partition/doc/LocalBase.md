# `src/LocalBase.hs`

Generic helpers used by every other module. The module has no export list, so everything here is exported.

## Types and classes

- **`Dist`**: a newtype over `Int` for distances. It isn't used by the current algorithms.
- **`Ori = LR | RL`**: orientation, left-to-right or right-to-left.
- **`class Orientable o`**:
  - `getOri` gives the orientation of a value.
  - `invOri` inverts the value.
  - Instances: `Gene` (negative means `RL`; `invOri` negates), `GenesIRs` and the newtypes built on it (orientation of a whole chromosome, see [Genomes](Genomes.md)).
- **`canonicOri`**: returns the value in `LR` orientation, applying `invOri` if needed. For genes this is `abs`. It's used whenever genes are compared regardless of sign.

## List helpers

| Function | Behaviour |
|----------|-----------|
| `unique` | Removes duplicates (via `HashSet`, so order isn't kept) |
| `lPairs l` | Consecutive pairs `zip l (tail l)` |
| `interleavelists l1 l2` | `a1 b1 a2 b2 …`; used to put `|` separators or breakpoint IRs between blocks |
| `evens` / `odds` | Elements at even / odd positions (0-based) |
| `replace i f xs` | Applies `f` to the element at position `i` (0-based); `indexError` if out of range |
| `rotateL k l` | `drop k l ++ take k l` |
| `minWith f` / `maxWith f` | `Maybe` minimum/maximum by key. `maxWith` uses `maximumBy`, which returns the **last** maximum on ties |
| `sortWith f` | Stable sort by key, so equal keys keep their order |
| `localHead` | `head` that fails with `patternError` instead of a partial-function warning |
| `skipFirst` | `drop 1` |
| `fpow n f x` | Applies `f` `n` times |
| `flipMap` | Swaps keys and values of a `Map`; the values must be distinct |
| `xOr` | Boolean exclusive or |

## Error messages

The algorithms call `error` with one of these strings, so a crash message tells which kind of failure happened:

- `patternError`: a case that should be impossible.
- `indexError`: an index out of bounds.
- `inputError`: invalid input.
- `logicError`: the algorithm reached a state its logic should prevent.
- `undefinedError`: a case that isn't implemented.

## Debugging

`traceValue` / `traceValueS prefix` print a value via `Debug.Trace` and return it unchanged.
