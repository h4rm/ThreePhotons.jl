include("runs.jl")

#K2 must be at least such that rmax > 35.0 Angstrom (Phasing)
deltar = 1.0 #2 Angstrom resolution maximum
qmax2 = pi/deltar
K2_min = ceil(Int64, 35.0 / deltar)
for K2 = K2_min:3:K2_min+12
  deltaq = qmax2 / K2
  for K3 = 20:2:min(K2,40)
    for L3 = 12:2:floor(Int64,(K2-1)/2)
      run_optimal(K2, K3, L3)
    end
  end
end
