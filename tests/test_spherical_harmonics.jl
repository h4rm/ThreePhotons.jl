using Base.Test
using ThreePhotons

#Also tested in test_structure.jl

function test_Dlms(l, m, mb)
    l1 = l
    l2 = l
    m1 = m
    m2 = m
    m1b = mb
    m2b = mb

    sum = 0
    da = 0.1
    for theta = 0:da:pi
      for phi = 0:da:2*pi
        for gamma = 0:da:2*pi
          sum += sin(theta)* Dlms(l1, m1, m1b, theta, phi, gamma)* conj(Dlms(l2, m2, m2b, theta, phi, gamma)) * da^3
        end
      end
    end

    @test_approx_eq_eps real(sum) 8*pi^2/(2*l+1) 2e-2
end

test_Dlms(4,2,2)
