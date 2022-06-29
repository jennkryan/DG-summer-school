
julia --project=~/Downloads/jennifer/MSIAC  (You can also go to the folder and then run julia ] activate . )

julia > include("MSIAC_pkg.jl")

press ]

instantiate

build

note: WriteVTK and VTKDataIO may have some precompilation issues (FINE)

press backspace

---- exit Julia

WARNING: If there was a compilation problem, it should have already occured


export JULIA_NUM_THREADS=8 ( check your actual number of cores)

julia --project=~/Downloads/jennifer/MSIAC

cd ~/Downloads/jennifer/MSIAC/tutorials/

julia> include("tut1.jl")
[ Info: Precompiling MSIAC [89000846-9b51-4315-9022-878a679b62f2]
[ Info: Skipping precompilation since __precompile__(false). Importing MSIAC [89000846-9b51-4315-9022-878a679b62f2].

...
