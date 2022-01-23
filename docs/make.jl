push!(LOAD_PATH,"../src/")
using Documenter, Huggett

makedocs(modules = [Huggett], sitename = "Huggett.jl")

deploydocs(repo = "github.com/FGerding/Huggett.jl.git", devbranch = "main")
