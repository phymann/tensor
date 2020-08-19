using ITensors
using QuadGK

let
    # define basic parameters of the system
    T = 3.0     # temperature
    dMAX = 20     # maximum number of singular values
    topscale = 6    # top-most scale we want to reach with TRG

    # define dimension of indices for ITensor
    dim0 = 2

    # difine an initial Index making up the Ising partition function
    s = Index(dim0, "scale=0")

    # Boltzmann weight tensor "A"
    l = addtags(s, "left")
    r = addtags(s, "right")
    u = addtags(s, "up")
    d = addtags(s, "down")

    A = ITensor(l,r,u,d)

    # convert {1,2} (coming from indices) to {1,-1} (used in Boltzmann weight)
    function Sig(x)
        1.0 - 2.0*(x-1)
    end
    # fill the A tensor with corresponding Boltzmann weights
    for sl in (1:dim0)
        for sd in (1:dim0)
            for sr in (1:dim0)
                for su in (1:dim0)
                    E = Sig(sl)*Sig(sd)+Sig(sd)*Sig(sr)+Sig(sr)*Sig(su)+Sig(su)*Sig(sl);
                    P = exp(-E/T)
                    A[l(sl),d(sd),r(sr),u(su)] = P
                end
            end
        end
    end

    # keep tract of partition function per site, z = Z^(1/N_site)
    z = 1.0

    for scale in (1:topscale)
        println("------- Scale = ",scale," -> ",scale+1,"-------")

        # A = Fl * Fr factorization
        Fl, Fr = factorize(A,(r,d); maxdim = dMAX, tags = "left,scale="*string(scale))
        # grab the new left/right Index
        l_new = commonind(Fl,Fr)
        r_new = replacetags(l_new,"left","right")

        # A = Fu * Fd factorization
        Fu, Fd = factorize(A,(l,d); maxdim = dMAX, tags = "up,scale="*string(scale))
        # grab the new up/down Index
        u_new = commonind(Fu,Fd)
        d_new = replacetags(u_new,"up","down")

        # change Fr's tags from left to right
        Fr *= δ(l_new,r_new)
        # change Fd's tags from up to down
        Fd *= δ(u_new,d_new)

        # set bonds between tensors
        Fl *= delta(r,l)
        Fu *= delta(d,u)
        Fr *= delta(l,r)
        Fd *= delta(u,d)

        # perform tensor contraction
        A = Fl * Fu * Fr * Fd

        println(inds(A))

        # update the indices
        l = l_new
        r = r_new
        u = u_new
        d = d_new

        # normalize the current tensor (which amounts to normalize the partition function)
        # and keep track of the total normalization

        TrA = scalar(A * delta(l,r) * delta(u,d))
        A /= TrA;
        z *= TrA^(1.0/(2^(1+scale)))
    end

    println("log(Z)/N = ",log(z));

    # exact result
    k = sinh(2/T)^(-2)      # a numerical parameter
    exactf = log(2)/2 + 1/(2π) * (quadgk(x -> log(cosh(2/T)^2+sqrt(1 + k^2 - 2*k*cos(2x))/k), 0, π)[1])
    println("exact result = ",exactf)
end