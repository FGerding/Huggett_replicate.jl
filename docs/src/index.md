# My Replication of Asset Prices in a Huggett economy
# This is joint work with Matteo Cremonini

> This replication study was part of my evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at Bocconi in Winter 2021/2022

In this replication study, our main contribution has been to code a procedure which plots and computes the excess demand function of a Huggett economy, for any interest rate level and for any borrowing constraint level.

The package exports all values needed to replicate the paper into the user's Environment, while the file Graphs.jl plots the two figures of the paper.

```@docs
AssetGrid
crra
initguess
next_guess
household
lambda
eD
```

end
