module Huggett


using BenchmarkTools
using LinearAlgebra
using Plots
using GeometryTypes
using QuantEcon

export AssetGrid, crra, initguess, next_guess, household, lambda, eD


mutable struct interest ; max::Float64; min::Float64; r::Float64 ; end 

export interest, β, σ, borrow, itmx_m, tol_v, tol_m, itmx_v, vmin, eps_clear, a_max, na, ns, ϵ_h, ϵ_l, y, Π, r_store, a_eye, rmax, rmin, r, V, TV, pol, PP, excess, VI, income, q, r_seq, eD_seq_A0, eD_seq_A1, eD_seq_A2, sig_seq, q_star, q_star_0
    β = 0.993362
	σ = 1.5       
	borrow = -2
	itmx_m = 100
	tol_v = 1e-4
	tol_m = 1e-3
	itmx_v = 5000
	vmin = -1e10
	eps_clear = 1e-3
	a_max = 24
	na = 500
	ns = 2
	ϵ_h = 1
	ϵ_l = 0.1
	y = [ϵ_h, ϵ_l]
	Π = [0.925 0.075; 0.5 0.5]
	r_store = zeros(itmx_m,1)
	a_eye = Matrix{Int}(I, na, na)      
	rmax = -0.01
	rmin = -0.02
    r = interest(rmax,rmin,rmax)
    a_eye = Matrix{Int}(I, na, na)      
	V = zeros(ns*na,1)
    TV = V
    pol = zeros(ns*na)
	PP = zeros(ns*na,ns*na)
	excess = []
	VI = zeros(ns,na)
	income = vcat(repeat([y[1]],na,1),repeat([y[2]],na,1))
    q = sort(unique(Iterators.flatten((range(1,1.5,length=20),
    range(1.5,4,length=6)))))
    r_seq = sort(unique(Iterators.flatten((range(1,1.5,length=20),
    range(1.5,4,length=6)))))
    eD_seq_A0 = sort(unique(Iterators.flatten((range(1,1.5,length=20),
    range(1.5,4,length=6)))))
    eD_seq_A1 = sort(unique(Iterators.flatten((range(1,1.5,length=20),
    range(1.5,4,length=6)))))
    eD_seq_A2 = sort(unique(Iterators.flatten((range(1,1.5,length=20),
    range(1.5,4,length=6)))))
    sig_seq = sort(unique(Iterators.flatten((range(0,10,length=20)))))
	q_star = sort(unique(Iterators.flatten((range(0,10,length=20)))))
	q_star_0 = sort(unique(Iterators.flatten((range(0,10,length=20)))))


""" 
`function AssetGrid(borrow)`


Creates a grid of equally spaced points ranging from the minimum level of assets, i.e. the borrowing constraint level borrow (which has to be specified by the user). 
The maximum number of assets a_max = 24 and the number of assets na=500 have been specified for the sake of the replication exercise.

# Inputs:

- `borrow::Float64` : Borrowing constraint level for the Huggett economy. Notice that borrow<=0 must hold since in the Huggett economy considered the assets are not in positive net supply. 

# Outputs:

- Column vector of dimension na of type Matrix{Float64}
"""
function AssetGrid(borrow)
	
    # create an asset grid based on the borrowing constraint
        
        a_vals = range(borrow,a_max,length=na)
        a_grid = gridmake(a_vals)
    
        return a_grid
end


""" 
`function crra(c,σ)`


Computes the utility associated to a level of consumption c, when the utility function is a CRRA utility with constant relative risk aversion parameter equal to σ. Remember that a CRRA utility function takes the form c^{1-σ}/1-σ

# Inputs:

- `c::Float64`: the level of consumption for which we want to compute the associated utility. Notice that `c>=0` must hold.
- `σ::Float64`: relative risk aversion parameter of the CRRA utility function. 

# Outputs:

- A scalar of type Float64 indicating the utility level.
""" 
function crra(c,σ)
# create the utility function

        if c>= 0
            u = c^(1-σ)/(1-σ)
        else u = vmin
        end
end


""" 
`function initguess(r, borrow)`


Computes an initial guess for the Value Function to start the Value Function Iteration (VFI) procedure used to solve the household problem of the Huggett economy. It is not necessary to run this function in order to solve the household problem, since the Contraction Mapping Theorem guarantees the convergence of the Value Function Iteration procedure, starting from any Value Function Guess. However, running this function makes the code faster. 

# Inputs:

- `r::interest`: the interest rate used to use to solve the household’s problem.
- `borrow::Float64` : the borrowing constraint level to be specified by the user.

# Outputs:

- `ns x na` vector of type Vector{Float64} to be used as a starting guess of the VFI procedure, where na is the number of assets and ns is the number of states for the income process.


""" 
function initguess(r, borrow)

        # using a smarter guess for the value function to make it computationally faster
        
            C = repeat(AssetGrid(borrow),ns,na).*(1+r.r) + kron(y,ones(na,na)) - repeat(AssetGrid(borrow)',ns*na,1)
            U = crra.(C,σ) 
            vtemp =[]
            for j=1:ns
            vtemp = U[(j-1)*na+1:j*na,:]
            V[(ns-1)*na+1:ns*na]=diag(vtemp)./(1-β);
            end
        return V
end


""" 
`function next_guess(V,pol, borrow)`


Implements the Howard’s improvement algorithm, which is used to update the value function after each step of the value function iteration procedure.  The basic idea of Howard’s Improvement Step is that instead of finding a maximizer for each iteration of the VFI procedure, the same maximizer should be used repeatedly for a specified  number of iterations (e.g. 70). The resulting Value Function will be the starting point of the following step of the VFI procedure. 

# Inputs:

- `V::Vector{Float64}` : Vector of dimension ns x na which represents the Value function obtained at the end of any step of the VFI procedure.
- `pol::Vector{Int64}` : Policy Function is a column vector of length ns x na coming from a given iteration of the VFI procedure.
- `borrow::Float64` : Borrowing constraint level to be specified by the user.

# Outputs:

The output of the function is a Value Function of dimension ns x na and type Vector{Float64} to be used in the next iteration of the VFI procedure
""" 
function next_guess(V,pol, borrow)

    # Update the value function guess through Howard's improvement algorithm
    
        global griglia = repeat(AssetGrid(borrow),ns,1).*(1+r.r)
        Vr = reshape(V,na,ns)'
        Vnew = zeros(ns,na)
        for h in 1:70
        Vnew =zeros(ns,na)
        utils = reshape(crra.((griglia+income-AssetGrid(borrow)[pol]),σ),na,ns)'
        for iasset in 1:na
            Vnew[:,iasset] = utils[:,iasset] +β*diag(Π*Vr[:,(reshape(pol,na,ns)')[:,iasset]])
        end
        Vr[:,:] = Vnew
        end
    return Vnew
end


""" 
`function household(v,r, borrow)`


Implements the VFI iteration procedure in order to solve the household problem of the Huggett economy for a given interest rate `r` and given a level of borrowing constraint borrow. The Value function obtained from   the function `initguess` is used as an input, serving as the initial guess in the VFI procedure. Furthermore, the function calls next guess at the end of each step of the VFI procedure: in this way Howard’s improvement algorithm is used in order to obtain the next guess of the VFI procedure. 

# Inputs: 

- `v::Vector{Float64}` : Initial guess for Value Function.
- `r::interest` : Interest rate to be used for solving the household problem.
- `borrow::Float64` : Borrowing constraint level to be specified by the user.

# Output: 

- A tuple of the following type: Tuple{Vector{Float64}, Vector{Int64}, Float64}  where the first element of the Tuple is the Value Function, the second element in the Tuple is the Policy function, i.e. a vector of length ns x na .  The third element of the tuple, instead, is a scalar of Type Float64 indicating the size of the discrepancy between the value functions of two successive iterations. If this parameter is under the tolerance level set then the VFI procedure has converged correctly. 
""" 
function household(v,r, borrow)

    # solving the household problem
        
    C = repeat(AssetGrid(borrow),ns,na).*(1+r.r) + kron(y,ones(na,na)) - repeat(AssetGrid(borrow)',ns*na,1)
    U = crra.(C,σ)
    Vnew = zeros(ns,na)
        
        for i1=1:5000
        res = reshape(v,na,ns)'
        EV = Π*res
        sol = findmax(U+β.* kron(EV, ones(na,1)), dims = 2)
        TV = sol[1]
        pol = getindex.(sol[2],2)
        dV = findmax(abs.(TV-v))[1]
            if dV < tol_v
                return (v,pol,dV)
            else v = TV
            end
        Vnew[:,:] = next_guess(v,pol, borrow)
        v = vcat(Vnew[1,:],Vnew[2,:])
    end
    (v,TV,Vnew)
end

findnearest(A,x) = argmin(abs.(A .- x))

		
""" 
`function lambda(pol_new)`


Taking as input the policy function obtained solving Huggett’s household problem, this function computes the implied stationary distribution of assets in the economy. 

# Inputs:

- `pol_new::Vector{Int64}` : the Policy Function is a column vector of length ns x na coming from a given iteration of the VFI procedure.

# Outputs:

- The output is a vector of length na of type Vector{Float64}, displaying the probability of having any given level of assets. Stated differently, the output is the wealth (or asset) distribution across the agents in the economy. 
""" 
function lambda(pol_new)
# Compute the asset distribution  
            for j=1:ns
                for l=1:ns
                index1 = (j-1)*na+1:j*na
                index2 = (l-1)*na+1:l*na
                PP[index1,index2] = Π[j,l]*a_eye[pol_new[index1],:]
                end
            end
            eigD = eigvals(PP')
            eigV = eigvecs(PP')
            lambda =convert(Array{Float64},eigV[:,findnearest(Diagonal(eigD),1)[1]]/sum(eigV[:,findnearest(Diagonal(eigD),1)[1]]))
end
		
""" 
`function eD(pol_new,lambda, borrow)`


This function retrieves the amount of excess demand in the economy given the Policy Function, the asset’s stationary distribution and the level of borrowing constraints. Stated in an equivalent way, given any level of interest rate `r` this function takes the policy function coming from the `household` function, and the induced stationary distribution of assets coming from the `lambda` function in order to compute the excess demand of the considered Huggett economy (with level of borrowing constraint borrow).

# Inputs:

- `pol_new::Vector{Int64}` : the Policy Function is a column vector of length ns x na coming from the solution of the household problem when the given interest rate is r.
- `lambda::Vector{Int64}` : the stationary distribution of assets in the Huggett economy. 
- `borrow::Float64` : the borrowing constraint level imposed on the Huggett economy.

# Outputs:

- A scalar of type Float64 indicating the excess demand in the considered Huggett economy. Notice that if excess demand = 0 then we the interest rate which delivers pol_new as a policy function is the general equilibrium interest rate!
"""
function eD(pol_new,lambda, borrow)

                # calculating the excess demand
                    
                A = AssetGrid(borrow)[pol_new]' * lambda
                dA = A
end
end
