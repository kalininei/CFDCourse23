x1 = [0.0]
w1 = [2.0]

x2 = [-0.5773502691896257, 0.5773502691896257]
w2 = [1, 1]

x3 = [
    -0.774596669241483377035853079956,
    0.0,
    0.774596669241483377035853079956
]
w3 = [
    0.555555555555555555555555555555,
    0.888888888888888888888888888888,
    0.555555555555555555555555555555
]

x4 = [
  -0.8611363115940526,
  -0.3399810435848563,
  0.3399810435848563,
  0.8611363115940526
]
w4 = [
    0.3478548451374538,
    0.6521451548625461,
    0.6521451548625461,
    0.3478548451374538
]

x = x4
w = w4

# ====== 1D
n = len(x)
# weights
print("{")
for i in range(n):
    e = "," if i < n-1 else ""
    print(f"\t{w[i]}{e}")
print("}")
# points
print("{")
for i in range(n):
    e = "," if i < n-1 else ""
    print("\tPoint{" + f"{x[i]}, 0, 0" + "}" + f"{e}")
print("}")

# ====== 2D
n = len(x)
# weights
print("{")
g = 0
for j in range(n):
    for i in range(n):
        g += 1
        e = "," if g < n*n else ""
        print(f"\t{w[i]*w[j]}{e}")
print("}")
# points
print("{")
g = 0
for j in range(n):
    for i in range(n):
        g += 1
        e = "," if g < n*n else ""
        print("\tPoint{" + f"{x[i]}, {x[j]}, 0" + "}" + f"{e}")
print("}")

# ====== 3D
n = len(x)
# weights
g = 0
print("{")
for k in range(n):
    for j in range(n):
        for i in range(n):
            g += 1
            e = "," if g < n*n*n else ""
            print(f"\t{w[i]*w[j]*w[k]}{e}")
print("}")
# points
g = 0
print("{")
for k in range(n):
    for j in range(n):
        for i in range(n):
            g += 1
            e = "," if g < n*n*n else ""
            print("\tPoint{" + f"{x[i]}, {x[j]}, {x[k]}" + "}" + f"{e}")
print("}")
