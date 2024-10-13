function ωjk(x,knotvector,j,k)
    return (x-knotvector[j])/(knotvector[j+k-1]-knotvector[j])
end
function Br(x,knotvector,r)
    if knotvector[r] < x < knotvector[r+1]
        return 1.0
    else
        return 0.0
    end
end


#include("./spline.jl")
elem    = 5
no      = elem+1
l,L     = 1.0,6.0
dn      = (L-l)/(n-1)
knotset = collect(l:dn:L)
np      = 100
dp      = (L-l)/(np)
xp      = collect(l:dp:L)

k       = 1
B       = zeros(no-k,np)

no = 1
for (np,x) ∈ enumerate(xp)
    B[no,np] = ωjk(x,knotvector,no,k)
end

