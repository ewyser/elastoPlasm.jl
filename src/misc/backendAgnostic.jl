# https://juliagpu.github.io/KernelAbstractions.jl/stable/
# include("./src/misc/backendAgnostic.jl")
using KernelAbstractions, CUDA, LinearAlgebra


@kernel inbounds=true function mat_kernel!(D)
    i,j = @index(Global, NTuple)
    if 1<i<size(D.A,1) && 1<j<size(D.A,2)
        D.A[i,j] = D.b + D.A[i,j]
    end
end
@kernel inbounds=true function vec_kernel!(D)
    i = @index(Global)
    if 1<i<length(D.a)
        D.a[i] = 2*D.b + D.a[i]
    end
end

@kernel inbounds=true function test!(sig,np)
    k = @index(Global)
    if k<=np
        λ,n = eigen(sig[:,:,k],sortby=nothing)
    end
end

function configDevice(dev)
    blockSize     = KernelAbstractions.isgpu(dev) ? 256 : Threads.nthreads()
    global mat_k! = mat_kernel!(dev, blockSize)
    global vec_k! = vec_kernel!(dev, blockSize)
    return @info "kernel launch parameters & alias(es) done"
end

# initialize
A     = ones(2048,2048)#A     = CuArray(A)
nx,ny = size(A)

# allocate memory, get backend & kernel launch parameter
AT = CUDA.functional() ? CUDA.CuArray : Array
A_D= AT(A)

dev = get_backend(A_D)
configDevice(CPU())

D = (A = A_D,a = Array(ones(100,)), b = 10.0,)

# agnostic backend kernel call
for k in 1:5
    @time mat_k!(D;ndrange=size(D.A))
    KernelAbstractions.synchronize(dev)
end
# print result

for k in 1:1
    vec_k!(D;ndrange=length(D.a))
    KernelAbstractions.synchronize(dev)
end
# print result



sig = ones(3,3,10)
sig_D = sig
global eigN! = test!(CPU(), KernelAbstractions.isgpu(dev) ? 256 : 1024)
np = 10
eigN!(sig_D,np, ndrange=np)
KernelAbstractions.synchronize(dev)