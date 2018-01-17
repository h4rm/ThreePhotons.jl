using Base.Test
using ThreePhotons

function test_random_rotations()
  p1 = [0,1,0]
  p2 = rand(3)#[0,0,1]#/sqrt(2)
  s1 = 0
  s2 = 0
  for i = 1:10000
      rot = random_rotation(3)
      s1 += rot*p1
      s2 += rot*p2
  end
  println(s1/10000)
  println(s2/10000)
end

function test_random_rotations2()
  val = sum([det(random_rotation(3)) for i = 1:10000])
  println(val)
end

#Tests if random_rotation for D=3 samples all rotations uniformly
function test_complettness_of_rotations(random_rotation_generator=random_rotation)

  angles = [get_euler_angles(random_rotation_generator(3)) for l = 1:100000]
  a = Base.hist(convert(Array{Float64},map((x)->x[1], angles)), 100)
  b = Base.hist(convert(Array{Float64},map((x)->x[2], angles)), 100)
  c = Base.hist(convert(Array{Float64},map((x)->x[3], angles)), 100)

  x1 = collect(a[1])
  x2 = collect(b[1])
  x3 = collect(c[1])
  plot(x1[1:length(x1)-1], a[2]/sum(abs, a[2]), label="phi")
  plot(x2[1:length(x2)-1], b[2]/sum(abs, b[2]), label="theta")
  plot(x3[1:length(x3)-1], c[2]/sum(abs, c[2]), label="gamma")
  xlabel("Angle [rad]")
  ylabel("Distribution of Angles")
  legend()
  # [a[1];]

end

function test_random_rotations3(func=random_rotation)
    global veclist = Float64[];
    init = [1,0,0]
    append!(veclist, init)
    num = int(1e3)
    for i = 1:num
        append!(veclist, func(3)*init)
    end
    # global z = zeros(3,num)
    global v = reshape(veclist, 3, num+1)
    scatter3D( v[1,:], v[2,:], v[3,:])
    xlim(-1, 1)
    ylim(-1,1)
    zlim(-1,1)
end

function test_random_rotations4(func=random_rotation)
    global veclist = Float64[];
    init = [1,0]
    append!(veclist, init)
    num = int(1e3)
    for i = 1:num
        append!(veclist, func(2)*init)
    end
    # global z = zeros(3,num)
    global v = reshape(veclist, 2, num+1)
    scatter( v[1,:], v[2,:])
    xlim(-1, 1)
    ylim(-1,1)
    # zlim(-1,1)
end

function test_random_roation_stepper()
    global veclist = Float64[];
    init = [1,0,0]
    num = int(1e4)
    for i = 1:num
        init = random_rotation_step(3, pi/8)*init
        append!(veclist, init)
    end
    # global z = zeros(3,num)
    global v = reshape(veclist, 3, num)
    scatter3D( v[1,:], v[2,:], v[3,:])
    xlim(-1, 1)
    ylim(-1,1)
    zlim(-1,1)
end

"""Tests SON_parametrized by throwing random angles and plotting the point cloud in 3D. This should give a seamlessly but ununiform sampling of the 3D sphere."""
function test_SON_paramtrized()
  r = 0:0.1:pi
  list = []
  for i = 1:1000
    push!(list, [0 0 1] * SON_parametrized(3, [2*pi*rand(), 2*pi*rand(), 2*pi*rand()]))
  end
  x = [list[i][1] for i=1:length(list)]
  y = [list[i][2] for i=1:length(list)]
  z = [list[i][3] for i=1:length(list)]
  scatter3D(x,y,z)
end
