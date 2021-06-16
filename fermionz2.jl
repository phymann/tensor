### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 7c78137c-cd7a-11eb-0557-237d0cfc4d25
begin
	using ITensors
	using Plots
	using PlutoUI
	plotly()
end

# ╔═╡ b29f1131-d3a5-4c18-b4c1-09e154c0dcf6
begin
	N = 90
	
	# sites
	sites = [ isodd(mod1(n,3)) ? siteind("Fermion",n,conserve_qns=true) : Index([QN()=>2];tags="S=1/2,Site,n=$n") for n in 1:N]
	
	# construct H
	jv = 0.3
	jh = 1.0
	h = 0.05
	# μ = 0.1
	
	ampo = OpSum()
		# ff horizontal
		for j=1:3:N-5
			ampo += (jh,"Cdag",j,"C",j+3)
			ampo += (jh,"Cdag",j+3,"C",j)
			ampo += (jh,"Cdag",j+2,"C",j+5)
			ampo += (jh,"Cdag",j+5,"C",j+2)
		end
		# fzf vertical
		for j=1:3:N-2
			ampo .+= (jv,"Cdag",j+2,"C",j,"Sz",j+1)
			ampo .+= (jv,"Cdag",j,"C",j+2,"Sz",j+1)
		end
		# sx term
		for j=2:3:N-1
			ampo += (h,"Sx",j)
		end
		# # μ term
		# for j=1:3:3N-2
		# 	ampo += (-μ,"N",j)
		# 	ampo += (-μ,"N",j+2)
		# end
	
	H = MPO(ampo,sites);
	
	# construct initial states
	states = [
	if mod(n,3)==1 && mod(n,6)==1
		"Occ"
	elseif mod(n,3)==2 && iseven(n)
		"Up"
	elseif mod(n,3)==2 && isodd(n)
		"Dn"
	else
		"Emp"
	end
	for n in 1:N];
	
	psi0 = MPS(sites,states);
	
	#DMRG
	sweeps = Sweeps(7)
	setmaxdim!(sweeps,10,20,80,200,300,400)
	setmindim!(sweeps,5)
	setcutoff!(sweeps,1E-5,1E-6,1E-7,1E-8)
	setnoise!(sweeps,1E-5,1E-5,1E-8,1E-9,1E-10)

	energy,psi = dmrg(H,psi0,sweeps)
	
	nothing
end

# ╔═╡ 78953ca8-6a3c-44cc-a7d9-95028d7f916a
# σz
begin
	mz = zeros(div(N,3))
	for j=2:3:N-1
		orthogonalize!(psi,j)
		s = siteind(psi,j)
		mz[div(j+1,3)] = scalar(psi[j]*op(s,"Sz")*dag(prime(psi[j],s)))
	end
	plot([1:div(N,3)],mz,ylims=(-1,1),xlabel="site i",ylabel="〈Sz〉")
end

# ╔═╡ e20298af-0e0a-4696-a86c-088605a24bd2
# σx
begin
	mx = zeros(div(N,3))
	for j=2:3:N-1
		orthogonalize!(psi,j)
		s = siteind(psi,j)
		mx[div(j+1,3)] = scalar(psi[j]*op(s,"Sx")*dag(prime(psi[j],s)))
	end
	plot([1:div(N,3)],mx,ylims=(-1,1),xlabel="site i",ylabel="〈Sx〉")
end

# ╔═╡ 3759fcd2-3065-4e1b-996f-b0f2f18f7a64
# nL
begin
	nl = zeros(div(N,3))
	for j=3:3:N
		orthogonalize!(psi,j)
		s = siteind(psi,j)
		nl[div(j,3)] = scalar(psi[j]*op(s,"N")*dag(prime(psi[j],s)))
	end
	plot([1:div(N,3)],nl,ylims=(-0.1,1.1),xlabel="site i",ylabel="〈n_l〉")
end

# ╔═╡ ac21cc4b-ad0a-4772-999c-0467bd61aa1a
# nL
begin
	nr = zeros(div(N,3))
	for j=1:3:N-2
		orthogonalize!(psi,j)
		s = siteind(psi,j)
		nr[div(j+2,3)] = scalar(psi[j]*op(s,"N")*dag(prime(psi[j],s)))
	end
	plot([1:div(N,3)],nr,ylims=(-0.1,1.1),xlabel="site i",ylabel="〈n_r〉")
end

# ╔═╡ 36994bd4-59a5-481c-bc6a-0555157d766c
# f single-particle correlator for left leg
# a naive way
begin
	corr = zeros(div(N,3),div(N,3))
	for i=1:div(N,3),j=i:div(N,3)
		a = OpSum()
			a += "Cdag",3i,"C",3j
		corr[i,j] = inner(psi,MPO(a,sites),psi)
		corr[j,i] = conj(corr[i,j])
	end
end

# ╔═╡ 0963dfde-c625-41fb-8c57-31cf2c210ac2
plot([1:div(N,3)],corr[div(N,6),:])

# ╔═╡ 48c0c4e4-7776-4d68-8000-1ab7830a979a
# fermion single-particle correlator for left leg
# a faster way
begin
	orthogonalize!(psi,1)
	corr1 = zeros(div(N,3),div(N,3))
	L = ITensor(1.)
	for i=1:N
		if mod(i,3)==0
			Li = L*psi[i]*op("Cdag",sites,i)*dag(prime(psi[i]))
			LiF = Li
			for j=i+3:3:N
				LiF *= psi[j-2]*op("F",sites,j-2)*dag(prime(psi[j-2]))
				LiF *= psi[j-1]*dag(prime(psi[j-1],"Link"))
				lind = commonind(psi[j],LiF)
				LiF *= psi[j]
				cij = LiF*op("C",sites,j)*dag(prime(prime(psi[j],"Site"),lind))
				corr1[div(i,3),div(j,3)] = scalar(cij)
				corr1[div(j,3),div(i,3)] = conj(corr1[div(i,3),div(j,3)])
				LiF *= op("F",sites,j)*dag(prime(psi[j]))
			end
		end
		L *= psi[i]*dag(prime(psi[i],"Link"))
	end
				
	# fill in diagonal elements
	for i=1:div(N,3)
		orthogonalize!(psi,3*i)
		corr1[i,i] = scalar(psi[3i]*op("N",sites,3*i)*dag(prime(psi[3i],"Site")))
	end
end

# ╔═╡ 9807da0b-09c2-4154-90fd-ab96f8f51bba
norm(corr-corr1)

# ╔═╡ 35784c8c-a362-4e54-9471-61b429265022
plot([1:div(N,3)],corr1[div(N,6),:])

# ╔═╡ Cell order:
# ╠═7c78137c-cd7a-11eb-0557-237d0cfc4d25
# ╠═b29f1131-d3a5-4c18-b4c1-09e154c0dcf6
# ╠═78953ca8-6a3c-44cc-a7d9-95028d7f916a
# ╠═e20298af-0e0a-4696-a86c-088605a24bd2
# ╠═3759fcd2-3065-4e1b-996f-b0f2f18f7a64
# ╠═ac21cc4b-ad0a-4772-999c-0467bd61aa1a
# ╠═36994bd4-59a5-481c-bc6a-0555157d766c
# ╠═0963dfde-c625-41fb-8c57-31cf2c210ac2
# ╠═48c0c4e4-7776-4d68-8000-1ab7830a979a
# ╠═9807da0b-09c2-4154-90fd-ab96f8f51bba
# ╠═35784c8c-a362-4e54-9471-61b429265022
