# Translating DFM's python version:
include("terms.jl")

mutable struct CeleriteGP{T,TTerm}
  kernel::TTerm
  computed::Bool
  D::Vector{T}
  W::Matrix{T}
  up::Matrix{T}
  phi::Matrix{T}
  x::Vector{T}
  var::Vector{T}
  logdet::T
  n::Int
  J::Int
  function CeleriteGP(kernel::Term{T}) where T
    new{T,typeof(kernel)}(kernel, false, zeros(Float64, 0), zeros(Float64, 0, 0), zeros(Float64, 0, 0), zeros(Float64, 0, 0), zeros(Float64, 0), zeros(Float64, 0), 0.0, 0, 0)
  end
end

function cholesky_ldlt!(a_real::Vector{Float64}, c_real::Vector{Float64},
  a_comp::Vector{Float64}, b_comp::Vector{Float64},
  c_comp::Vector{Float64}, d_comp::Vector{Float64},
  t::Vector{Float64}, var::Vector{Float64}, X::Array{Float64,2},
  phi::Array{Float64,2}, u::Array{Float64,2}, D::Vector{Float64})
  #
  # Fast LDLT Cholesky solver based on low-rank decomposition due to Sivaram, plus
  # real implementation of Celerite term.
  #
  # Compute the dimensions of the problem:
  N = length(t)
  # Number of real components:
  J_real = length(a_real)
  # Number of complex components:
  J_comp = length(a_comp)
  # Rank of semi-separable components:
  J = J_real + 2 * J_comp
  # phi is used to stably compute exponentials between time steps:
  phi = _reshape!(phi, J, N - 1)
  # u, X & D are low-rank matrices and diagonal component:
  u = _reshape!(u, J, N)
  X = _reshape!(X, J, N)
  D = _reshape!(D, N)
  diag = _reshape!(D, N)

  # Sum over the diagonal kernel amplitudes:    
  a_sum = sum(a_real) + sum(a_comp)
  # Compute the first element:
  #    D[1] = sqrt(var[1] + a_sum)
  diag[1] = var[1] + a_sum
  D[1] = diag[1]
  value = 1.0 / D[1]
  for j in 1:J_real
    u[j, 1] = a_real[j]
    X[j, 1] = value
  end
  # We are going to compute cosine & sine recursively - allocate arrays for each complex
  # component:
  cd::Vector{Float64} = zeros(J_comp)
  sd::Vector{Float64} = zeros(J_comp)
  # Initialize the computation of X:
  for j in 1:J_comp
    cd[j] = cos(d_comp[j] * t[1])
    sd[j] = sin(d_comp[j] * t[1])
    u[J_real+2*j-1, 1] = a_comp[j] * cd[j] + b_comp[j] * sd[j]
    u[J_real+2*j, 1] = a_comp[j] * sd[j] - b_comp[j] * cd[j]
    X[J_real+2*j-1, 1] = cd[j] * value
    X[J_real+2*j, 1] = sd[j] * value
  end
  # Allocate array for recursive computation of low-rank matrices:   
  S::Array{Float64,2} = zeros(J, J)
  #    for j in 1:J
  #      for k in 1:j
  #        S[k,j] = X[k,1]*X[j,1]
  #      end
  #    end
  # Allocate temporary variables:
  phij = 0.0
  dx = 0.0
  dcd = 0.0
  dsd = 0.0
  cdtmp = 0.0
  uj = 0.0
  Xj = 0.0
  Dn = 0.0
  Sk = 0.0
  uk = 0.0
  tmp = 0.0
  tn = 0.0
  # Loop over remaining indices:
  @inbounds for n in 2:N
    # Update phi
    tn = t[n]
    # Using time differences stabilizes the exponential component and speeds
    # up cosine/sine computation:
    dx = tn - t[n-1]
    # Compute real components of the low-rank matrices:
    for j in 1:J_real
      phi[j, n-1] = exp(-c_real[j] * dx)
      u[j, n] = a_real[j]
      X[j, n] = 1.0
    end
    # Compute complex components:
    for j in 1:J_comp
      value = exp(-c_comp[j] * dx)
      phi[J_real+2*j-1, n-1] = value
      phi[J_real+2*j, n-1] = value
      cdtmp = cd[j]
      dcd = cos(d_comp[j] * dx)
      dsd = sin(d_comp[j] * dx)
      cd[j] = cdtmp * dcd - sd[j] * dsd
      sd[j] = sd[j] * dcd + cdtmp * dsd
      # Update u and initialize X
      u[J_real+2*j-1, n] = a_comp[j] * cd[j] + b_comp[j] * sd[j]
      u[J_real+2*j, n] = a_comp[j] * sd[j] - b_comp[j] * cd[j]
      X[J_real+2*j-1, n] = cd[j]
      X[J_real+2*j, n] = sd[j]
    end

    # Update S
    Dn = D[n-1]
    for j in 1:J
      phij = phi[j, n-1]
      Xj = X[j, n-1]
      for k in 1:j
        S[k, j] = phij * phi[k, n-1] * (S[k, j] + Dn * Xj * X[k, n-1])
      end
    end

    # Update D and X
    Dn = 0.0
    for j in 1:J
      uj = u[j, n]
      Xj = X[j, n]
      for k in 1:j-1
        Sk = S[k, j]
        tmp = uj * Sk
        uk = u[k, n]
        Dn += uk * tmp
        Xj -= uk * Sk
        X[k, n] -= tmp
      end
      tmp = uj * S[j, j]
      Dn += 0.5 * uj * tmp
      X[j, n] = Xj - tmp
    end
    # Finalize computation of D:
    diag[n] = var[n] + a_sum
    Dn = diag[n] - 2.0 * Dn
    D[n] = Dn
    # Finalize computation of X:
    for j in 1:J
      X[j, n] /= Dn
    end
    # Update S
    #        Xj = 0.0
    #        for j in 1:J
    #            Xj = X[j,n]
    #            for k in 1:j
    #                S[k, j] += Xj*X[k, n]
    #            end
    #        end
  end
  # Finished looping over n.  Now return components to the calling routine
  # so that these may be used in arithmetic:
  return D, X, u, phi
end


function cholesky!(a_real::Vector, c_real::Vector,
  a_comp::Vector, b_comp::Vector,
  c_comp::Vector, d_comp::Vector,
  t::Vector, diag::Vector, X::Matrix,
  phi::Matrix, u::Matrix, D::Vector)

  T = promote_type(
    eltype(a_real),
    eltype(c_real),
    eltype(a_comp),
    eltype(b_comp),
    eltype(c_comp),
    eltype(d_comp),
    eltype(t),
    eltype(diag),
    eltype(X),
    eltype(phi),
    eltype(u),
    eltype(D),
  )
  a_real = convert(Vector{T}, a_real)
  c_real = convert(Vector{T}, c_real)
  a_comp = convert(Vector{T}, a_comp)
  b_comp = convert(Vector{T}, b_comp)
  c_comp = convert(Vector{T}, c_comp)
  d_comp = convert(Vector{T}, d_comp)
  t = convert(Vector{T}, t)
  diag = convert(Vector{T}, diag)
  X = convert(Matrix{T}, X)
  phi = convert(Matrix{T}, phi)
  u = convert(Matrix{T}, u)
  D = convert(Vector{T}, D)
  #
  # Fast Cholesky solver based on low-rank decomposition due to Sivaram, plus
  # real implementation of Celerite term.
  #
  # Compute the dimensions of the problem:
  N = length(t)
  # Number of real components:
  J_real = length(a_real)
  # Number of complex components:
  J_comp = length(a_comp)
  # Rank of semi-separable components:
  J = J_real + 2 * J_comp
  # phi is used to stably compute exponentials between time steps:
  phi = _reshape!(phi, J, N - 1)
  # u, X & D are low-rank matrices and diagonal component:
  u = _reshape!(u, J, N)
  X = _reshape!(X, J, N)
  D = _reshape!(D, N)


  # Sum over the diagonal kernel amplitudes:    
  a_sum = sum(a_real) + sum(a_comp)
  # Compute the first element:
  D[1] = sqrt(diag[1] + a_sum)
  value = 1.0 / D[1]
  for j in 1:J_real
    u[j, 1] = a_real[j]
    X[j, 1] = value
  end
  # We are going to compute cosine & sine recursively - allocate arrays for each complex
  # component:
  cd::Vector{T} = zeros(J_comp)
  sd::Vector{T} = zeros(J_comp)
  # Initialize the computation of X:
  for j in 1:J_comp
    cd[j] = cos(d_comp[j] * t[1])
    sd[j] = sin(d_comp[j] * t[1])
    u[J_real+2*j-1, 1] = a_comp[j] * cd[j] + b_comp[j] * sd[j]
    u[J_real+2*j, 1] = a_comp[j] * sd[j] - b_comp[j] * cd[j]
    X[J_real+2*j-1, 1] = cd[j] * value
    X[J_real+2*j, 1] = sd[j] * value
  end
  # Allocate array for recursive computation of low-rank matrices:   
  S::Array{T,2} = zeros(J, J)
  for j in 1:J
    for k in 1:j
      S[k, j] = X[k, 1] * X[j, 1]
    end
  end
  # Allocate temporary variables:
  phij = 0.0
  dx = 0.0
  dcd = 0.0
  dsd = 0.0
  cdtmp = 0.0
  uj = 0.0
  Xj = 0.0
  Dn = 0.0
  Sk = 0.0
  uk = 0.0
  tmp = 0.0
  tn = 0.0
  # Loop over remaining indices:
  @inbounds for n in 2:N
    # Update phi
    tn = t[n]
    # Using time differences stabilizes the exponential component and speeds
    # up cosine/sine computation:
    dx = tn - t[n-1]
    # Compute real components of the low-rank matrices:
    for j in 1:J_real
      phi[j, n-1] = exp(-c_real[j] * dx)
      u[j, n] = a_real[j]
      X[j, n] = 1.0
    end
    # Compute complex components:
    for j in 1:J_comp
      value = exp(-c_comp[j] * dx)
      phi[J_real+2*j-1, n-1] = value
      phi[J_real+2*j, n-1] = value
      cdtmp = cd[j]
      dcd = cos(d_comp[j] * dx)
      dsd = sin(d_comp[j] * dx)
      cd[j] = cdtmp * dcd - sd[j] * dsd
      sd[j] = sd[j] * dcd + cdtmp * dsd
      # Update u and initialize X
      u[J_real+2*j-1, n] = a_comp[j] * cd[j] + b_comp[j] * sd[j]
      u[J_real+2*j, n] = a_comp[j] * sd[j] - b_comp[j] * cd[j]
      X[J_real+2*j-1, n] = cd[j]
      X[J_real+2*j, n] = sd[j]
    end

    # Update S
    for j in 1:J
      phij = phi[j, n-1]
      for k in 1:j
        S[k, j] = phij * phi[k, n-1] * S[k, j]
      end
    end

    # Update D and X
    Dn = 0.0
    for j in 1:J
      uj = u[j, n]
      Xj = X[j, n]
      for k in 1:j-1
        Sk = S[k, j]
        tmp = uj * Sk
        uk = u[k, n]
        Dn += uk * tmp
        Xj -= uk * Sk
        X[k, n] -= tmp
      end
      tmp = uj * S[j, j]
      Dn += 0.5 * uj * tmp
      X[j, n] = Xj - tmp
    end
    # Finalize computation of D:
    Dn = sqrt(diag[n] + a_sum - 2.0 * Dn)
    D[n] = Dn
    # Finalize computation of X:
    for j in 1:J
      X[j, n] /= Dn
    end
    # Update S
    Xj = 0.0
    for j in 1:J
      Xj = X[j, n]
      for k in 1:j
        S[k, j] += Xj * X[k, n]
      end
    end
  end
  # Finished looping over n.  Now return components to the calling routine
  # so that these may be used in arithmetic:
  return D, X, u, phi
end

function compute_ldlt!(gp::CeleriteGP, x, yerr=0.0)
  # Call the choleksy function to decompose & update
  # the components of gp with X,D,V,U,etc. 
  coeffs = get_all_coefficients(gp.kernel)
  var = yerr .^ 2 + zeros(Float64, length(x))
  gp.n = length(x)
  gp.D, gp.W, gp.up, gp.phi = cholesky_ldlt!(coeffs..., x, var, gp.W, gp.phi, gp.up, gp.D)
  gp.J = size(gp.W)[1]
  # Compute the log determinant (square the determinant of the Cholesky factor):
  #  gp.logdet = sum(log(gp.D))
  logdet = 0.0
  for i = 1:gp.n
    logdet += log(gp.D[i])
  end
  #  gp.logdet = sum(log(gp.D))
  gp.logdet = logdet
  gp.x = x
  gp.computed = true
  return gp.logdet
end

function compute!(gp::CeleriteGP, x, yerr=0.0)
  # Call the choleksy function to decompose & update
  # the components of gp with X,D,V,U,etc. 
  coeffs = get_all_coefficients(gp.kernel)
  if length(gp.var) == 0
    gp.var = yerr .^ 2 + zeros(Float64, length(x))
  else
    gp.var .= yerr .^ 2 + zeros(Float64, length(x))
  end
  gp.n = length(x)
  #  @time gp.D,gp.W,gp.up,gp.phi = cholesky!(coeffs..., x, var, gp.W, gp.phi, gp.up, gp.D)
  gp.D, gp.W, gp.up, gp.phi = cholesky!(coeffs..., x, gp.var, gp.W, gp.phi, gp.up, gp.D)
  gp.J = size(gp.W)[1]
  # Compute the log determinant (square the determinant of the Cholesky factor):
  gp.logdet = 2 * sum(log.(gp.D))
  #  gp.logdet = sum(log(gp.D))
  gp.x = x
  gp.computed = true
  return gp.logdet
end

function invert_lower_ldlt(gp::CeleriteGP, y)
  # Applies just the lower inverse:  L^{-1}.y:
  @assert(gp.computed)
  N = gp.n
  @assert(length(y) == N)
  z = zeros(Float64, N)
  # The following lines solve L.z = y for z:
  z[1] = y[1]
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* z[n-1])
    #    f .= gp.phi[:,n] .* (f .+ gp.W[:,n-1] .* z[n-1])
    z[n] = y[n] - dot(gp.up[:, n], f)
  end
  return z
end

function invert_lower(gp::CeleriteGP, y)
  # Applies just the lower inverse:  L^{-1}.y:
  @assert(gp.computed)
  N = gp.n
  @assert(length(y) == N)
  z = zeros(Float64, N)
  # The following lines solve L.z = y for z:
  z[1] = y[1] / gp.D[1]
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* z[n-1])
    z[n] = (y[n] - sum(gp.up[:, n] .* f)) / gp.D[n]
  end
  return z
end

function apply_inverse_ldlt(gp::CeleriteGP, y)
  # Solves for K.b=y for b with LDLT decomposition.
  @assert(gp.computed)
  N = gp.n
  @assert(length(y) == N)
  z = zeros(Float64, N)
  # The following lines solve L.z = y for z:
  #  z[1] = y[1]/gp.D[1]
  z[1] = y[1]
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* z[n-1])
    z[n] = (y[n] - dot(gp.up[:, n], f))
  end
  # The following solves L^T.z = y for z:
  y = copy(z)
  fill!(z, zero(Float64))
  z[N] = y[N] / gp.D[N]
  fill!(f, zero(Float64))
  for n = N-1:-1:1
    f .= gp.phi[:, n] .* (f .+ gp.up[:, n+1] .* z[n+1])
    z[n] = y[n] / gp.D[n] - dot(gp.W[:, n], f)
  end
  # The result is the solution of L.L^T.z = y for z,
  # or z = {L.L^T}^{-1}.y = L^{T,-1}.L^{-1}.y
  return z
end

function apply_inverse(gp::CeleriteGP, y)
  # Solves for K.b=y for b.
  @assert(gp.computed)
  N = gp.n
  @assert(length(y) == N)
  z = zeros(Float64, N)
  # The following lines solve L.z = y for z:
  z[1] = y[1] / gp.D[1]
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= @views gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* z[n-1])
    z[n] = (y[n] - dot(gp.up[:, n], f)) / gp.D[n]
  end
  # The following solves L^T.z = y for z:
  y = copy(z)
  z = zeros(Float64, N)
  z[N] = y[N] / gp.D[N]
  f = zeros(Float64, gp.J)
  for n = N-1:-1:1
    f .= @views gp.phi[:, n] .* (@views(f .+ gp.up[:, n+1] .* z[n+1]))
    z[n] = @views (y[n] - dot(gp.W[:, n], f)) / gp.D[n]
  end
  # The result is the solution of L.L^T.z = y for z,
  # or z = {L.L^T}^{-1}.y = L^{T,-1}.L^{-1}.y
  return z
end

function simulate_gp_ldlt(gp::CeleriteGP, z)
  # Multiplies lower Cholesky factor times random Gaussian vector (z is N(1,0) ) to simulate
  # a Gaussian process.
  # Check that Cholesky factor has been computed
  @assert(gp.computed)
  # Check for length mismatch:
  N = gp.n
  @assert(length(z) == N)
  # If iid is zeros, then draw from random normal deviates:
  if maximum(abs.(z)) == 0.0
    randn!(z) # In Julia 0.6 this can be replaced with "z .= randn.()"
  end
  y = zeros(Float64, N)
  # Carry out multiplication
  tmp = sqrt(gp.D[1]) * z[1]
  y[1] = tmp
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* tmp)
    tmp = sqrt(gp.D[n]) * z[n]
    y[n] = tmp + dot(gp.up[:, n], f)
  end
  # Return simulated correlated noise vector, y=L.z
  return y
end

function multiply_ldlt(gp::CeleriteGP, x, z, yerr=0.0)
  # Multiplies the full matrix times a vector, z.
  # Need to compute diagonal, so get coefficients:
  a1, c1, a2, b2, c2, d2 = get_all_coefficients(gp.kernel)
  # Number of real components:
  J1 = length(a1)
  # Number of complex components:
  J2 = length(a2)
  # Rank of semi-separable components:
  J = J1 + 2 * J2
  N = length(x)
  @assert(length(z) == N)
  # Carry out multiplication
  asum = sum(a1) + sum(a2)
  phi = zeros(J, N - 1)
  u = zeros(J, N - 1)
  v = zeros(J, N - 1)
  cd = cos.(d2 .* x[1])
  sd = sin.(d2 .* x[1])
  dx = 0.0
  for n = 1:N-1
    dx = x[n+1] - x[n]
    if J1 > 0
      v[1:J1, n] .= 1.0
      u[1:J1, n] .= a1
      phi[1:J1, n] .= exp.(-c1 .* dx)
    end
    if J2 > 0
      v[J1+1:J1+J2, n] .= cd
      v[J1+J2+1:J, n] .= sd
      cd .= cos.(d2 .* x[n+1])
      sd .= sin.(d2 .* x[n+1])
      u[J1+1:J1+J2, n] .= a2 .* cd .+ b2 .* sd
      u[J1+J2+1:J, n] .= a2 .* sd .- b2 .* cd
      phi[J1+1:J1+J2, n] .= exp.(-c2 .* dx)
      phi[J1+J2+1:J, n] .= phi[J1+1:J1+J2, n]
    end
  end
  # This is A_{n,n} z_n from paper:
  y = (yerr .^ 2 .+ asum) .* z
  # sweep upwards in n:
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= phi[:, n-1] .* (f .+ v[:, n-1] .* z[n-1])
    # This is \sum_{j=1}^J \tilde U_{n,j} f^-_{n,j}
    y[n] += dot(u[:, n-1], f)
  end
  # sweep downwards in n:
  fill!(f, zero(Float64))
  for n = N-1:-1:1
    f .= phi[:, n] .* (f .+ u[:, n] .* z[n+1])
    # This is \sum_{j=1}^J \tilde U_{n,j} f^-_{n,j}
    y[n] += dot(v[:, n], f)
  end
  # Return result of multiplication:
  return y
end

function simulate_gp(gp::CeleriteGP, y)
  # Multiplies Cholesky factor times random Gaussian vector (y is N(1,0) ) to simulate
  # a Gaussian process.
  # If iid is zeros, then draw from random normal deviates:
  # Check that Cholesky factor has been computed
  # Carry out multiplication
  # Return simulated correlated noise vector
  N = gp.n
  @assert(length(y) == N)
  z = zeros(Float64, N)
  z[1] = gp.D[1] * y[1]
  f = zeros(Float64, gp.J)
  for n = 2:N
    f .= gp.phi[:, n-1] .* (f .+ gp.W[:, n-1] .* y[n-1])
    z[n] = gp.D[n] * y[n] + dot(gp.up[:, n], f)
  end
  # Returns z=L.y
  return z
end

function log_likelihood_ldlt(gp::CeleriteGP, y)
  # O(N) log likelihood computation once the low-rank Cholesky decomposition is completed
  @assert(gp.computed)
  if size(y, 2) != 1
    error("y must be 1-D")
  end
  alpha = apply_inverse_ldlt(gp, y)
  nll = gp.logdet + gp.n * log(2 * pi)
  for i in 1:gp.n
    nll = nll + alpha[i] * y[i]
  end
  return -0.5 * nll
end

function log_likelihood(gp::CeleriteGP, y)
  # O(N) log likelihood computation once the low-rank Cholesky decomposition is completed
  @assert(gp.computed)
  if size(y, 2) != 1
    error("y must be 1-D")
  end
  alpha = apply_inverse(gp, y)
  nll = gp.logdet + gp.n * log(2 * pi)
  for i in 1:gp.n
    nll = nll + alpha[i] * y[i]
  end
  return -0.5 * nll
end

function full_solve(t::Vector, y0::Vector, aj::Vector, bj::Vector, cj::Vector, dj::Vector, yerr::Vector)
  # This carries out the full GP solver using linear algebra on full covariance matrix.
  # WARNING: do not use this with large datasets.
  N = length(t)
  J = length(aj)
  @assert(length(y0) == length(t))
  u = zeros(Float64, N, J * 2)
  v = zeros(Float64, N, J * 2)
  # Compute the full U/V matrices:
  for i = 1:N
    for j = 1:J
      expcjt = exp(-cj[j] * t[i])
      cosdjt = cos(dj[j] * t[i])
      sindjt = sin(dj[j] * t[i])
      u[i, j*2-1] = aj[j] * expcjt * cosdjt + bj[j] * expcjt * sindjt
      u[i, j*2] = aj[j] * expcjt * sindjt - bj[j] * expcjt * cosdjt
      v[i, j*2-1] = cosdjt / expcjt
      v[i, j*2] = sindjt / expcjt
    end
  end
  # Diagonal components:
  #  diag = fill(yerr.^2 + sum(aj),N)
  diag = yerr .^ 2 .+ sum(aj)

  # Compute the kernel:
  K = zeros(Float64, N, N)
  for i = 1:N
    for j = 1:N
      abs_ti_tj = abs(t[i] - t[j])
      K[i, j] = sum(aj .* exp.(-cj .* abs_ti_tj) .* cos.(dj .* abs_ti_tj) .+ bj .* exp.(-cj .* abs_ti_tj) .* sin.(dj .* abs_ti_tj))
    end
    K[i, i] = diag[i]
  end

  # Check that equation (1) holds:
  K0 = tril(*(u, v'), -1) + triu(*(v, u'), 1)
  for i = 1:N
    K0[i, i] = diag[i]
  end
  println("Semiseparable error: ", maximum(abs.(K .- K0)))
  return logdet(K), K
end

function predict_ldlt!(gp::CeleriteGP, t, y, x)
  # Predict future times, x, based on a 'training set' of values y at times t.
  # Runs in O((M+N)J^2) (variance is not computed, though)
  a_real, c_real, a_comp, b_comp, c_comp, d_comp = get_all_coefficients(gp.kernel)
  N = length(y)
  M = length(x)
  J_real = length(a_real)
  J_comp = length(a_comp)
  J = J_real + 2 * J_comp

  b = apply_inverse_ldlt(gp, y)
  Q = zeros(J)
  X = zeros(J)
  pred = similar(x)
  fill!(pred, 0)

  # Forward pass
  m = 1
  while m < M && x[m] <= t[1]
    m += 1
  end
  for n = 1:N
    if n < N
      tref = t[n+1]
    else
      tref = t[N]
    end
    Q[1:J_real] = (Q[1:J_real] .+ b[n]) .* exp.(-c_real .* (tref - t[n]))
    Q[J_real+1:J_real+J_comp] .= (Q[J_real+1:J_real+J_comp] .+ b[n] .* cos.(d_comp .* t[n])) .*
                                 exp.(-c_comp .* (tref - t[n]))
    Q[J_real+J_comp+1:J] .= (Q[J_real+J_comp+1:J] .+ b[n] .* sin.(d_comp .* t[n])) .*
                            exp.(-c_comp .* (tref - t[n]))

    while m < M + 1 && (n == N || x[m] <= t[n+1])
      X[1:J_real] = a_real .* exp.(-c_real .* (x[m] - tref))
      X[J_real+1:J_real+J_comp] .= a_comp .* exp.(-c_comp .* (x[m] - tref)) .* cos.(d_comp .* x[m]) .+
                                   b_comp .* exp.(-c_comp .* (x[m] - tref)) .* sin.(d_comp .* x[m])
      X[J_real+J_comp+1:J] .= a_comp .* exp.(-c_comp .* (x[m] - tref)) .* sin.(d_comp .* x[m]) .-
                              b_comp .* exp.(-c_comp .* (x[m] - tref)) .* cos.(d_comp .* x[m])

      pred[m] = dot(X, Q)
      m += 1
    end
  end

  # Backward pass
  m = M
  while m >= 1 && x[m] > t[N]
    m -= 1
  end
  fill!(Q, 0.0)
  for n = N:-1:1
    if n > 1
      tref = t[n-1]
    else
      tref = t[1]
    end
    Q[1:J_real] .= (Q[1:J_real] .+ b[n] .* a_real) .*
                   exp.(-c_real .* (t[n] - tref))
    Q[J_real+1:J_real+J_comp] .= (Q[J_real+1:J_real+J_comp] .+ b[n] .* a_comp .* cos.(d_comp .* t[n]) .+
                                  b[n] .* b_comp .* sin.(d_comp .* t[n])) .*
                                 exp.(-c_comp .* (t[n] - tref))
    Q[J_real+J_comp+1:J] .= (Q[J_real+J_comp+1:J] .+ b[n] .* a_comp .* sin.(d_comp .* t[n]) .-
                             b[n] .* b_comp .* cos.(d_comp .* t[n])) .*
                            exp.(-c_comp .* (t[n] - tref))

    while m >= 1 && (n == 1 || x[m] > t[n-1])
      X[1:J_real] .= exp.(-c_real .* (tref - x[m]))
      X[J_real+1:J_real+J_comp] .= exp.(-c_comp .* (tref - x[m])) .* cos.(d_comp .* x[m])
      X[J_real+J_comp+1:J] .= exp.(-c_comp .* (tref - x[m])) .* sin.(d_comp .* x[m])

      pred[m] += dot(X, Q)
      m -= 1
    end
  end
  return pred
end

function predict!(gp::CeleriteGP, t, y, x)
  # Predict future times, x, based on a 'training set' of values y at times t.
  # Runs in O((M+N)J^2) (variance is not computed, though)
  a_real, c_real, a_comp, b_comp, c_comp, d_comp = get_all_coefficients(gp.kernel)
  N = length(y)
  M = length(x)
  J_real = length(a_real)
  J_comp = length(a_comp)
  J = J_real + 2 * J_comp

  b = apply_inverse(gp, y)
  Q = zeros(J)
  X = zeros(J)
  pred = zeros(x)

  # Forward pass
  m = 1
  while m < M && x[m] <= t[1]
    m += 1
  end
  for n = 1:N
    if n < N
      tref = t[n+1]
    else
      tref = t[N]
    end
    Q[1:J_real] .= (Q[1:J_real] .+ b[n]) .* exp.(-c_real .* (tref .- t[n]))
    Q[J_real+1:J_real+J_comp] .= (Q[J_real+1:J_real+J_comp] .+ b[n] .* cos.(d_comp .* t[n])) .*
                                 exp.(-c_comp .* (tref - t[n]))
    Q[J_real+J_comp+1:J] .= (Q[J_real+J_comp+1:J] .+ b[n] .* sin.(d_comp .* t[n])) .*
                            exp.(-c_comp .* (tref - t[n]))

    while m < M + 1 && (n == N || x[m] <= t[n+1])
      X[1:J_real] .= a_real .* exp.(-c_real .* (x[m] - tref))
      X[J_real+1:J_real+J_comp] .= a_comp .* exp.(-c_comp .* (x[m] - tref)) .* cos.(d_comp .* x[m]) .+
                                   b_comp .* exp.(-c_comp .* (x[m] - tref)) .* sin.(d_comp .* x[m])
      X[J_real+J_comp+1:J] .= a_comp .* exp.(-c_comp .* (x[m] - tref)) .* sin.(d_comp .* x[m]) .-
                              b_comp .* exp.(-c_comp .* (x[m] - tref)) .* cos.(d_comp .* x[m])

      pred[m] = dot(X, Q)
      m += 1
    end
  end

  # Backward pass
  m = M
  while m >= 1 && x[m] > t[N]
    m -= 1
  end
  fill!(Q, 0.0)
  for n = N:-1:1
    if n > 1
      tref = t[n-1]
    else
      tref = t[1]
    end
    Q[1:J_real] .= (Q[1:J_real] .+ b[n] .* a_real) .*
                   exp.(-c_real .* (t[n] - tref))
    Q[J_real+1:J_real+J_comp] .= ((Q[J_real+1:J_real+J_comp] .+ b[n] .* a_comp .* cos.(d_comp .* t[n])) .+
                                  b[n] .* b_comp .* sin.(d_comp .* t[n])) .*
                                 exp.(-c_comp .* (t[n] - tref))
    Q[J_real+J_comp+1:J] .= (Q[J_real+J_comp+1:J] .+ b[n] .* a_comp .* sin.(d_comp .* t[n]) .-
                             b[n] .* b_comp .* cos.(d_comp .* t[n])) .*
                            exp.(-c_comp .* (t[n] - tref))

    while m >= 1 && (n == 1 || x[m] > t[n-1])
      X[1:J_real] .= exp.(-c_real .* (tref - x[m]))
      X[J_real+1:J_real+J_comp] .= exp.(-c_comp .* (tref - x[m])) .* cos.(d_comp .* x[m])
      X[J_real+J_comp+1:J] .= exp.(-c_comp .* (tref - x[m])) .* sin.(d_comp .* x[m])

      pred[m] += dot(X, Q)
      m -= 1
    end
  end
  return pred
end

function get_matrix(gp::CeleriteGP, xs...)
  # Gets the full covariance matrix.  Can provide autocorrelation or cross-correlation.
  # WARNING: Do not use with large datasets.
  if length(xs) > 2
    error("At most 2 arguments can be provided")
  end
  local x1::Array
  local x2::Array
  if length(xs) >= 1
    x1 = xs[1]
  else
    if !gp.computed
      error("You must compute the GP first")
    end
    x1 = gp.x
  end
  if length(xs) == 2
    x2 = xs[2]
  else
    x2 = x1
  end

  if size(x1, 2) != 1 || size(x2, 2) != 1
    error("Inputs must be 1D")
  end

  tau = broadcast(-, reshape(x1, length(x1), 1), reshape(x2, 1, length(x2)))
  return get_value(gp.kernel, tau)
end

function predict_full_ldlt(gp::CeleriteGP, y, t; return_cov=true, return_var=false)
  # Prediction with covariance using full covariance matrix
  # WARNING: do not use this with large datasets!
  alpha = apply_inverse_ldlt(gp, y)
  Kxs = get_matrix(gp, t, gp.x)
  mu = Kxs * alpha
  if !return_cov && !return_var
    return mu
  end

  KxsT = transpose(Kxs)
  if return_var
    #        println("Size of y: ",size(y)," size of t: ",size(t)," size of KxsT: ",size(KxsT))
    v = zeros(t)
    for i = 1:length(t)
      v[i] = -sum(KxsT[:, i] .* apply_inverse_ldlt(gp, KxsT[:, i]))
      if mod(i, 100) == 0
        #            println(i," ",v[i])
      end
    end
    v = v + get_value(gp.kernel, [0.0])[1]
    return mu, v[1, :]
  end

  cov = get_matrix(gp, t)
  cov = cov - Kxs * apply_inverse_ldlt(gp, KxsT)
  return mu, cov
end

function predict(gp::CeleriteGP, y, t; return_cov=true, return_var=false)
  # Prediction with covariance using full covariance matrix
  # WARNING: do not use this with large datasets!
  alpha = apply_inverse(gp, y)
  Kxs = get_matrix(gp, t, gp.x)
  mu = Kxs * alpha
  if !return_cov && !return_var
    return mu
  end

  KxsT = transpose(Kxs)
  if return_var
    v = zeros(eltype(t), size(t))
    for i = 1:length(t)
      #        v = -sum(KxsT .* apply_inverse(gp, KxsT), 1)
      v[i] = -sum(KxsT[:, i] .* apply_inverse(gp, KxsT[:, i]))
    end
    v = v .+ get_value(gp.kernel, [0.0])[1]
    #        return mu, v[1, :]
    return mu, v
  end

  cov = get_matrix(gp, t)
  cov = cov - Kxs * apply_inverse(gp, KxsT)
  return mu, cov
end

function predict_full(gp::CeleriteGP, y, t; return_cov=true, return_var=false)
  # Prediction with covariance using full covariance matrix
  # WARNING: do not use this with large datasets!
  alpha = apply_inverse(gp, y)
  Kxs = get_matrix(gp, t, gp.x)
  mu = Kxs * alpha
  if !return_cov && !return_var
    return mu
  end

  KxsT = transpose(Kxs)
  if return_var
    v = -sum(KxsT .* apply_inverse(gp, KxsT), 1)
    v = v + get_value(gp.kernel, [0.0])[1]
    return mu, v[1, :]
  end

  cov = get_matrix(gp, t)
  cov = cov - Kxs * apply_inverse(gp, KxsT)
  return mu, cov
end

function _reshape!(A::Array, dims...)
  # Allocates arrays if size is not correct
  if size(A) != dims
    A = Array{eltype(A)}(undef, dims)
  end
  return A
end
