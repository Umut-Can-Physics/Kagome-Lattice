# Quasihole braiding (09.10.2024)
# E. Kapit et al., PRL 108, 066802 (2012)

if (Lx == 4) && (Ly == 4)
    path = Array{Array{Int,1},1}(undef, 4)
    path[1] = [13, 14, 15, 7]
    path[2] = [7, 6, 5, 15]
    path[3] = [5, 9, 13, 15]
    path[4] = [15, 11, 7, 13]
    
    pathsimultaneous = Array{Array{Int,1},1}(undef, 2)
    pathsimultaneous[1] = [13, 14, 15, 11, 7]
    pathsimultaneous[2] = [7, 6, 5, 9, 13]

elseif (Lx == 6) && (Ly == 6)
    pathsimultaneous = Array{Array{Int,1},1}(undef, 2)
    pathsimultaneous[1] = fill(31, 25)
    pathsimultaneous[2] = [16, 22, 28, 34, 4, 10, 16, 15, 14, 13, 18, 17, 16, 10, 4, 34, 28, 22, 16, 17, 18, 13, 14, 15, 16]
end

nstep = parse(Int, readline(stdin))  # number of steps between two nearest neighbors

simultaneous = parse(Int, readline(stdin))  # simultaneous path? (0 or 1)
if simultaneous == 0
    npathL = 4
    inpathL = length(path[1]) - 2
    pathsize = 4 * (length(path[1]) - 2) * nstep
elseif simultaneous == 1
    npathL = 1
    inpathL = length(pathsimultaneous[1]) - 1
    pathsize = (length(pathsimultaneous[1]) - 1) * nstep
end
println(pathsize)

showfig = parse(Int, readline(stdin))  # show simulation (0 or 1)

Impurity = sparse(zeros(dim, dim))
randImpurity = sparse(zeros(dim, dim))

for si in 1:Ns
    randImpurity += rand() * n_op[si]
end

V_imp = parse(Float64, readline(stdin))  # V_imp
V_randImp = 1e-4
d = 2
UHcurrent = Array{Any,1}(undef, 2)
C = zeros(pathsize - 1)

for npath in 1:npathL
    for inpath in 1:inpathL
        for step in 1:nstep
            if simultaneous == 0
                Impurity = V_imp * ((1 - (step - 1) / nstep) * n_op[path[npath][inpath]] + (step - 1) / nstep * n_op[path[npath][inpath+1]] + n_op[path[npath][end]])
            elseif simultaneous == 1
                Impurity = V_imp * ((1 - (step - 1) / nstep) * (n_op[pathsimultaneous[npath][inpath]] +
                n_op[pathsimultaneous[npath+1][inpath]]) + 
                (step - 1) / nstep * (n_op[pathsimultaneous[npath][inpath+1]] + 
                n_op[pathsimultaneous[npath+1][inpath+1]]))
            end
            
            H = H0 + Interaction + Impurity + V_randImp * randImpurity
            HNp = H[projectNp, projectNp]
            if (npath == 1) && (inpath == 1) && (step == 1)
                EHplot = eigen(SparseMatrixCSC(HNp))
                plot(diag(EHplot), ".")
                degeneracy = parse(Int, readline(stdin))
                closeall()
                EHinitial = eigen(SparseMatrixCSC(HNp))[1:degeneracy]
                if showfig == 1
                    figure()
                end
                UHcurrent[1] = EHinitial
            else
                EH = eigen(SparseMatrixCSC(HNp))[1:degeneracy]
                UHcurrent[d] = EH
                
                A = UHcurrent[d-1]' * UHcurrent[d]
                C[d-1] = abs(det(A))
                
                # SVD
                U_svd, S_svd, V_svd = svd(A)
                M_svd = U_svd * V_svd'
                UHcurrent[d] = UHcurrent[d] * M_svd'
                UHcurrent[d-1] = nothing
                
                if showfig == 1
                    siteoccupation = zeros(Ns)
                    for si in 1:Ns
                        occ = sum([basis(si, projectNp) * abs(UH[:, ii])^2 for ii in 1:degeneracy])
                        siteoccupation[si] = occ / degeneracy
                    end
                    SiteOcc = reshape(siteoccupation, Ly, Lx)
                    imagesc(SiteOcc')
                    gca().ydir = "normal"
                    axis("square")
                    pause(0.02)
                elseif showfig == 0
                    if d % 10 == 0
                        println(d)
                    end
                end
                d += 1
            end
        end
    end
end

figure()
plot(C)

M = UHcurrent[d-1]' * UHinitial
eigenvalues = eigen(M)
Abs_eigenvalues = abs(eigenvalues)
Phases = angle(eigenvalues) / π