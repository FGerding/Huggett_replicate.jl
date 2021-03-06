# Hugget_Rep.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://FGerding.github.io/Huggett_replicate.jl/dev/)

## Installation

To install this package, please clone [this repository](https://github.com/FGerding/Huggett_replicate.jl) to your local files.To do this, type in your terminal:

    git clone https://github.com/FGerding/Huggett_replicate.jl Huggett_replicate

## Load this package

Start up Julia. Then, in package mode, activate a new environment and install all necessary dependencies for Hugget_replicate by typing

    activate .
    instantiate 

Then you are able to use this package by typing (back in Julia mode)
    
    using Hugget_replicate


## Use the package

The package exports all values needed to replicate the paper into the user's Environment, while the file Graphs.jl plots the two figures of the paper.
In order to reproduce Figures 1 and 2 of the paper you should:

1) Run the Huggett.jl file which loads all the values necessary for the replication and the necessary functions.
2) Run: import Pkg; Pkg.add("Huggett")
3) Run the graph.jl file which gives Figures 1 and 2!


## External documentation

The online documentation can be found [here](https://FGerding.github.io/Huggett_replicate.jl/dev/)

The paper we replicated can be found [here](https://www.sciencedirect.com/science/article/abs/pii/S002205311100055X)
