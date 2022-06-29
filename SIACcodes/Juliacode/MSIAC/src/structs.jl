
mutable struct FIELD
	zwC     ::Vector{gaussZW}
	zw      :: Vector{gaussZW}
	reverse :: Bool
	f       :: matD
	modes   :: matD
	name    :: String
	basis   :: EXPANSION
	function FIELD(type::String, p::Int, collapse::Bool; reverse=false::Bool,zw = nothing::Union{Nothing,Vector{gaussZW}},
		            f = nothing::Union{Nothing,matD}, m = nothing::Union{Nothing,matD},  n = "f"::String)
				   basis = EXPANSION(type,p,collapse)
				   qt = ["legendre","legendre"] ; Q = p + 2
		   		   collapse && (qt[2] = "Radau")
		   		   (zw == nothing) && (zw = [gaussZW(qt[1],p+1),gaussZW(qt[2],p+1)])
				   ff = (f == nothing && m != nothing) ? phys_values(basis,m,zw,reverse) : zeros(1,1)
				   mm = (m == nothing && f != nothing) ? modal_values(basis,f,zw,reverse) : zeros(1,1)
				   return new(zw,zw,reverse,ff,mm,n,basis)
	end
end

mutable struct EXPANSION
	type     :: String
	degree   :: Int
	J        :: Vector{Polynomials.Polynomial}
	mass     :: matD
	collapse :: Bool
	m1       :: vecI
	m2       :: vecI
	function EXPANSION(type::String,p::Int, collapse::Bool)
		if type == "Pk" || type == "pk" || type == "PK"
			m1 = vcat([fill(i,p+1-i) for i = 0:p]...)
			m2 = vcat([[0:p-i;] for i = 0:p]...)
			return new(type,p,Vector{Polynomials.Polynomial}(undef,1),ones(1,1),collapse,m1,m2)
		else
			m1,m2 = tensor_index(p+1)
			if type == "legendre" || type == "LEGENDRE" || type == "Legendre"
				PJ = [SP.basis.(SP.Legendre, i) for i = 0:p]
				J  = convert.(Polynomials.Polynomial, PJ)
			else
				Jac = [SP.Jacobi{1,1}(vcat(zeros(max(i,0)),1)) for i = 0:p-2]
				PJ  = convert.(Polynomials.Polynomial, Jac)
				J   = vcat(Polynomials.Polynomial([0.5,.-0.5]),
				          [Polynomials.Polynomial([0.25,0.0,.-0.25]) * PJ[i] for i = 1:p-1],
					       Polynomials.Polynomial([0.5,0.5]))
		    end
			# Get mass matrix
			q1    = gaussZW("legendre", p + 1)
			q2    = gaussZW("legendre", p + 1)
			t1,t2 = tensor_index(q1.Q , q2.Q)
			pl2   = zeros((p+1)^2,(p+1)^2)
			n1    = q1.nodes[t1] ;	n2    = q2.nodes[t2]
			we    = q1.weights[t1] .* q2.weights[t2]
			for i = 1:length(m1)
				p1 = J[m1[i]].(n1) .* J[m2[i]].(n2)
				for j = 1:length(m1)
					p2 = J[m1[j]].(n1) .* J[m2[j]].(n2)
					pl2[i,j] = sum(p1 .* p2 .* we)
				end
			end
			mij  = inv(pl2)
			return new(type,p,J,mij,collapse,m1,m2)
		end
	end
end

mutable struct DATAFIELD
	fields :: Vector{FIELD}
	fdict  :: Dict{Array{Int},Array{String}}
	mesh   :: MESH
	function DATAFIELD(f::Vector{FIELD}, m::MESH)
		fdict = Dict{Array{Int},Array{String}}()
		[fdict[k] = f[k].name  for k = 1: length(f)]
		return new(f,fdict,m)
	end
end
