using QuantumOptics
using LinearAlgebra

Lx = 5
pos = NLevelBasis(Lx)
H0 = sum([transition(pos,i,i+1)*exp(-im * i * pi / Lx) for i in 1:Lx-1])
H0 = H0 + dagger(H0)

#projection using projector operator:
Ncut = 3
E0s, states0 = eigenstates(dense(H0))
states = states0[1:Ncut]
pos_sub = SubspaceBasis(pos, states)
P1 = projector(pos_sub, pos) 

#projection without using projector operator:
E0s, states0 = eigen(dense(H0).data)
P2 = Operator(pos_sub, pos, states0[:,1:Ncut]')

display(P1.data)
display(P2.data)
P1==P2

#Pt1 = dagger(P1)
#H0_sub1 = P1*H0*Pt1

# The imaginary parts of the matrix elements of P1 and P2 have opposite signs
# @david_pl
# it might be due to taking transpose instead of adjoint
