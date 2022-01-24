using BenchmarkTools
using LinearAlgebra
using Plots
using GeometryTypes
using QuantEcon

using Huggett

begin
        for i in 1:length(q)
        r_seq[i] = (1/q[i])-1 			#interest rate based on bond prices
        end
            
        for i in 1:length(r_seq) 		
        #looping over each interest rate and calulate excess demand based on different borrowing contraints
            
        global r.r = r_seq[i] 			
        
        global borrow_temp0=0
        global temp_0 = household(initguess(r,borrow_temp0),r, borrow_temp0)[2]
        eD_seq_A0[i] = eD(temp_0,lambda(temp_0), borrow_temp0)[1]
        
        global borrow_temp1=-1
        global temp_1 = household(initguess(r,borrow_temp1),r, borrow_temp1)[2]
        eD_seq_A1[i] = eD(temp_1,lambda(temp_1), borrow_temp1)[1]
        
        global temp_2 = household(initguess(r, borrow),r, borrow)[2]
        eD_seq_A2[i] = eD(temp_2,lambda(temp_2), borrow)[1]
        
    end
end

begin
            #Replicating the first plot	
            zeroline = zeros(length(r_seq))
            plot(q, eD_seq_A0, color="blue", line=:dash, linewidth=5,label ="Borrowing constraint=0")
            plot!(q, eD_seq_A1,color="red", line=:dash, linewidth=5,label ="Borrowing constraint=-1")
            plot!(q, eD_seq_A2,color="green", line=:dash, linewidth=5,label ="Borrowing constraint=-2")
            plot!(q, zeroline,color="grey", label="")
            plot!([q[findfirst(isequal(0), eD_seq_A0)], 4], [0, 0], color="blue", linewidth=8, label="")
            scatter!([q[findfirst(isequal(0), eD_seq_A0)]],[0.],color="blue", markersize=8, label="q(0)")
            scatter!([q[findfirst(isone, eD_seq_A1.< 0.0)]],[0.],color="red", markersize=8, label="q(-1)")
            scatter!([q[findfirst(isone, eD_seq_A2.< 0.0)]],[0.],color="green", markersize=8, label="q(-2)")
            xlabel!("q")
            ylabel!("excess demand")
end

begin
    for i in 1:length(sig_seq) 		
        # instead of our parameter choices, plug in directly values used in the paper to 			create the plots
            q_star_0[i] = (0.98*( 0.9+(1- 0.9)*(1.08)^sig_seq[i]))*(1.01^(-sig_seq[i]))
            q_star[i] = 0.98*(1.01^(-sig_seq[i])) 
    end
    #Replicating the second plot
    plot(sig_seq, q_star, color="blue", line=:dash, linewidth=5,label ="complete")
    plot!(sig_seq, q_star_0, color="red", line=:dash, linewidth=5,label ="incomplete")
    xlabel!("σ")
    ylabel!("bond price")
end
