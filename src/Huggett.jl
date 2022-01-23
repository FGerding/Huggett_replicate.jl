module Huggett


using BenchmarkTools
using LinearAlgebra
using Plots
using GeometryTypes
using QuantEcon

export AssetGrid


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
Just a test
"""
function AssetGrid(borrow)
	
    # create an asset grid based on the borrowing constraint
        
        a_vals = range(borrow,a_max,length=na)
        a_grid = gridmake(a_vals)
    
        return a_grid
end

function crra(c,σ)
# create the utility function

        if c>= 0
            u = c^(1-σ)/(1-σ)
        else u = vmin
        end
end

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

function eD(pol_new,lambda, borrow)

                # calculating the excess demand
                    
                A = AssetGrid(borrow)[pol_new]' * lambda
                dA = A
end
end
